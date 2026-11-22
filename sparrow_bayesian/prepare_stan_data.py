"""Build the Stan data JSON for the Bayesian SPARROW member.

Inputs (the Ontario project's built products, in
``~/workspace/fei_lake_ontario/sparrow/``):

* ``reach_attributes_v1.csv`` -- per-reach sources, land cover and
  climate/soil attributes for the 23,963-reach extended network,
  extracted from the mid-continental SPARROW data1 table (Vouk et al.
  deposit, base year 2002).
* the same data1 table itself, for ``rhload`` (areal hydraulic load),
  which the attribute extraction did not carry over.
* ``station_year_loads_v4.csv`` -- annual TP loads 1990-2021 at the
  Balasmeh A/B stations (build_loads.py).  Each station's mean annual
  load over its record is one observation of the model.

Output: a Stan JSON matching ``sparrow.stan``'s data block, with
reaches in topological (upstream-to-downstream) order and the
downstream-index arrays (dummy sink = position N + 1) described in
the Stan file's header.

Station inclusion rule: only stations snapped ON the network enter
the fit.  The 15 off-network Central Lake Ontario CA stations are
excluded, because the source network contains no reaches in that
shore strip (see stations_snapped_v4.csv and the proposal package);
their station-year rows in the loads file (on_network == False) are
dropped here and the drop is reported when the script runs.

Observation standard errors (log space), per station, combine the two
components named in the design:

* the load-estimator error of the mean: sqrt(sum(se_t^2)) / n, from
  the per-year Beale / regression-flow standard errors, expressed
  relative to the mean load;
* the interannual variance of the station's annual loads: the
  standard deviation of log annual loads divided by sqrt(n).

Both are computed here and combined in quadrature into ``obs_se``;
the Stan model adds the fitted residual SD ``sigma`` on top.

Units: loads are converted to kg/yr inside the JSON (the loads file
is in t/yr); source columns of data1 are already kg/yr; land areas
are km^2 (incremental area x land-cover percent / 100).

Environment: mapvenv python (numpy / pandas).  No cmdstanpy.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from classical_sparrow import topological_order  # noqa: E402

SPARROW_DIR = os.path.expanduser("~/workspace/fei_lake_ontario/sparrow")
DATA1 = os.path.expanduser(
    "~/workspace/lake_of_the_woods/sparrow_midcontinental_nrc/work/"
    "midcont_SPARROW_input/data1_midcont.txt"
)
DEFAULT_OUT = os.path.join(SPARROW_DIR, "sparrow_stan_data.json")

VEL_FLOOR_MPS = 0.02  # guard against zero mean velocity in data1


def load_network() -> pd.DataFrame:
    """Extended-network reach table with rhload merged from data1."""
    a = pd.read_csv(os.path.join(SPARROW_DIR, "reach_attributes_v1.csv"))
    rh = pd.read_csv(DATA1, usecols=["sparrowid", "rhload"]).rename(
        columns={"sparrowid": "reach_sparrowid"}
    )
    a = a.merge(rh, on="reach_sparrowid", how="left")
    # In data1, rhload is populated exactly on reservoir reaches
    # (rchtype == 1) and is NaN on stream reaches; that is the source
    # table's own convention, not missing data.
    if a.loc[a["rchtype"] == 1, "rhload"].isna().any():
        raise ValueError("rhload missing for a reservoir reach after merge")
    a["rhload"] = a["rhload"].fillna(0.0)
    return a


def rhload_scale(df: pd.DataFrame) -> float:
    """Median rhload over reservoir reaches; the Stan input is rhload
    divided by this, so 1.0 is a typical reservoir reach and the
    reservoir coefficient's prior is on an interpretable scale."""
    return float(df.loc[df["rchtype"] == 1, "rhload"].median())


def _standardise(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    med = np.nanmedian(x)
    x = np.where(np.isnan(x), med, x)
    sd = x.std()
    return (x - x.mean()) / (sd if sd > 0 else 1.0)


def build_arrays(df: pd.DataFrame, scale: float | None = None) -> dict:
    """All per-reach Stan arrays for a reach table, in topological order.

    ``df`` needs the reach_attributes_v1 columns plus ``rhload``.
    Returns a dict of numpy arrays keyed by Stan data names (without
    the observation block), plus bookkeeping (reach ids in order,
    position lookup).  Used for the full network here and for the
    synthetic-test sub-basin in run_synthetic_test.py, so both go
    through exactly the same construction.
    """
    df = df.reset_index(drop=True).copy()
    order = topological_order(df["fnode"].to_numpy(), df["tnode"].to_numpy())
    d = df.iloc[order].reset_index(drop=True)
    n = len(d)
    pos_of_reach = {int(rid): i + 1 for i, rid in enumerate(d["reach_sparrowid"])}
    # reaches leaving each node, in topological position
    leaving: dict[int, list[int]] = {}
    for i, fn in enumerate(d["fnode"]):
        leaving.setdefault(int(fn), []).append(i + 1)
    down1 = np.full(n, n + 1, dtype=int)
    down2 = np.full(n, n + 1, dtype=int)
    for i, tn in enumerate(d["tnode"]):
        outs = leaving.get(int(tn), [])
        if len(outs) >= 1:
            down1[i] = outs[0]
        if len(outs) >= 2:
            down2[i] = outs[1]
        if len(outs) > 2:
            raise ValueError("network node with more than two leaving reaches")
    is_res = (d["rchtype"].to_numpy() == 1).astype(int)
    vel = np.clip(d["meanv_mps"].to_numpy(dtype=float), VEL_FLOOR_MPS, None)
    ttime = d["rch_len_m"].to_numpy(dtype=float) / (vel * 86400.0)
    ttime = np.where(is_res == 1, 0.0, ttime)
    if scale is None:
        scale = rhload_scale(df)
    rhload = np.where(is_res == 1,
                      d["rhload"].to_numpy(dtype=float) / scale, 0.0)
    area = d["demiarea"].to_numpy(dtype=float)
    # 639 reaches in data1 have zero incremental area and blank
    # land-cover / manure / climate fields; blank percent is zero area
    # in that class, so NaN percents are filled with 0 (the delivery
    # modifiers are median-filled inside _standardise instead, since a
    # zero-area reach still transports upstream load).
    def _pct(col: str) -> np.ndarray:
        return np.nan_to_num(d[col].to_numpy(dtype=float))
    arrays = dict(
        N=n,
        src_point=np.nan_to_num(d["pcs_p02_s"].to_numpy(dtype=float)),
        src_fert=np.nan_to_num(d["fert_p"].to_numpy(dtype=float)),
        src_man=np.nan_to_num(d["man_p_us_kg"].to_numpy(dtype=float)),
        area_urban=area * _pct("lulc_urb_perc") / 100.0,
        area_bg=area
        * (_pct("lulc_for_perc") + _pct("lulc_grsh_perc")
           + _pct("lulc_wet_perc") + _pct("lulc_past_perc"))
        / 100.0,
        area_water=area * _pct("lulc_wat_perc") / 100.0,
        z_runoff=_standardise(d["runoff_mm"].to_numpy(dtype=float)),
        z_perm=_standardise(d["soil_perm"].to_numpy(dtype=float)),
        z_tile=_standardise(d["tiles90_02_perc1"].to_numpy(dtype=float)),
        ttime=ttime,
        is_res=is_res,
        rhload=np.nan_to_num(rhload),
        frac=d["frac"].to_numpy(dtype=float),
        downstream1=down1,
        downstream2=down2,
    )
    return dict(arrays=arrays, reach_ids=d["reach_sparrowid"].to_numpy(),
                pos_of_reach=pos_of_reach, frame=d, rhload_scale=scale)


def station_observations(arrays_info: dict) -> pd.DataFrame:
    """Mean annual loads and log-space SEs for the on-network stations."""
    loads = pd.read_csv(
        os.path.join(SPARROW_DIR, "station_year_loads_v4.csv"),
        dtype={"station_id": str},
    )
    n_off = int((~loads["on_network"]).sum())
    loads = loads[loads["on_network"]].copy()
    print(f"station-year rows dropped as off-network (CLOCA): {n_off}")
    rows = []
    for sid, g in loads.groupby("station_id"):
        n = len(g)
        mean_t = float(g["load_t"].mean())
        se_mean_rel = float(np.sqrt((g["se_t"] ** 2).sum()) / n / mean_t)
        pos = g.loc[g["load_t"] > 0, "load_t"]
        if len(pos) >= 3:
            # Zero-load years exist in the loads file (one 0.1 km^2
            # station rounds to 0.0 t in dry years); the interannual
            # component is computed on the positive years only.
            inter = float(np.std(np.log(pos), ddof=1) / np.sqrt(len(pos)))
        else:
            inter = float(se_mean_rel)
        obs_se = float(np.hypot(se_mean_rel, inter))
        rows.append(
            dict(
                station_id=sid,
                reach_sparrowid=int(g["reach_sparrowid"].iloc[0]),
                n_years=n,
                mean_load_kg=mean_t * 1000.0,
                log_obs=float(np.log(mean_t * 1000.0)),
                obs_se=obs_se,
            )
        )
    obs = pd.DataFrame(rows)
    obs["obs_reach"] = obs["reach_sparrowid"].map(arrays_info["pos_of_reach"])
    if obs["obs_reach"].isna().any():
        raise ValueError("a station reach is missing from the network table")
    return obs


def to_stan_json(arrays: dict, obs: pd.DataFrame | None,
                 outlet_idx: list[int]) -> dict:
    out: dict = {}
    for k, v in arrays.items():
        out[k] = v.tolist() if isinstance(v, np.ndarray) else v
    if obs is None:
        out.update(n_obs=0, obs_reach=[], log_obs=[], obs_se=[])
    else:
        out.update(
            n_obs=int(len(obs)),
            obs_reach=[int(x) for x in obs["obs_reach"]],
            log_obs=[float(x) for x in obs["log_obs"]],
            obs_se=[float(x) for x in obs["obs_se"]],
        )
    out.update(n_outlet=len(outlet_idx),
               outlet_idx=[int(x) for x in outlet_idx])
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=DEFAULT_OUT)
    args = ap.parse_args()
    df = load_network()
    info = build_arrays(df)
    obs = station_observations(info)
    basins = pd.read_csv(os.path.join(SPARROW_DIR, "extended_basins_v1.csv"))
    outlet_idx = [info["pos_of_reach"][int(r)] for r in basins["outlet_sparrowid"]]
    data = to_stan_json(info["arrays"], obs, outlet_idx)
    with open(args.out, "w") as f:
        json.dump(data, f)
    print(f"reaches: {data['N']} | stations fitted: {data['n_obs']} "
          f"| basin outlets: {data['n_outlet']} "
          f"| rhload scale (reservoir median): {info['rhload_scale']:.6g}")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
