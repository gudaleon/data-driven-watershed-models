"""Python port of the GREEN model as implemented in the GREENeR R package.

GREENeR (GPL-3) is the JRC authors' implementation of GREEN (Geospatial
Regression Equation for European Nutrient losses; Grizzetti et al. 2005,
2008, 2012, 2021). This module ports the *core load computation* of
GREENeR 1.0.2 (file ``R/rgreen.R``: ``launch_green``, ``aggregate_loop``,
``green_shares``, ``region_nut_balance``, ``read_NSdata``, ``shreve``)
to numpy/pandas, faithfully to the R equations.

Deliberately NOT ported (see PORT_NOTES.md): the Latin-Hypercube
calibration (``calib_green``), goodness-of-fit metrics, and all plotting
/ mapping functions. Calibration of this port is intended to be
Bayesian (Stan) later.

Core equations (per catchment i, per year), exactly as in ``launch_green``:

    AA_i          = exp(-alpha_p * InvNrmRain_i)          # 1 - basin retention
    SumDS_i       = Min + Man + Fix + Soil + Atm * (1 - ForestFraction)
    SumPSCmp1_i   = Ps + atm_coeff * Atm * ForestFraction + sd_coef * Sd
    CatchToRiver_i = SumDS_i * AA_i + SumPSCmp1_i
    BB_i          = exp(-alpha_l * NrmLengthKm_i)         # 1 - river retention
    CatchLoad_i   = (1 - LakeFrRet_i) * (CatchToRiver_i + U_i) * BB_i

where U_i is the sum of CatchLoad over the catchments draining directly
into i, evaluated in topological order (GREENeR groups by its "Shreve"
levels, which are topological generations: headwaters = 1).

For the TP scenario GREENeR renames the background source ``Bg`` to
``Atm`` and sets ``Fix = Soil = 0`` (``append_empty_cols``), so the TP
background term passes through the same ``atm_coeff`` split as TN
atmospheric deposition. This port replicates the code, including that
behaviour (the vignette's TP equation prints the forest-fraction
background term without the 0.38 factor; the code applies it).

Input tables follow the GREENeR schemas:

catch_data columns:
    HydroID, To_catch (downstream catchment id, -1 for the outlet),
    Shreve (optional; recomputed if absent), LakeFrRet, NrmLengthKm

annual_data columns (TN):
    BasinID, YearValue, HydroID, NextDownID, Atm, Min, Man, Fix, Soil,
    Sd, Ps, YearlyMass, ForestFraction, InvNrmRain
annual_data columns (TP):
    BasinID, YearValue, HydroID, NextDownID, Bg, Min, Man, Sd, Ps,
    YearlyMass, ForestFraction, InvNrmRain

Units are those of the inputs (GREENeR examples use ton/yr); the model
is unit-agnostic because all retention terms are dimensionless.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

__all__ = [
    "compute_shreve_levels",
    "launch_green",
    "green_shares",
    "region_nut_balance",
    "read_nsdata",
    "TN_SOURCES",
    "TP_SOURCES",
]

TN_SOURCES = ["Atm", "Min", "Man", "Fix", "Soil", "Sd", "Ps"]
TP_SOURCES = ["Bg", "Min", "Man", "Sd", "Ps"]

_CATCH_COLS = ["HydroID", "To_catch", "LakeFrRet", "NrmLengthKm"]
_TN_COLS = ["BasinID", "YearValue", "HydroID", "NextDownID", "Atm", "Min",
            "Man", "Fix", "Soil", "Sd", "Ps", "YearlyMass",
            "ForestFraction", "InvNrmRain"]
_TP_COLS = ["BasinID", "YearValue", "HydroID", "NextDownID", "Bg", "Min",
            "Man", "Sd", "Ps", "YearlyMass", "ForestFraction", "InvNrmRain"]


# --------------------------------------------------------------------------
# Topology helpers
# --------------------------------------------------------------------------

def compute_shreve_levels(catch_data: pd.DataFrame) -> pd.DataFrame:
    """Replicates GREENeR's ``shreve()``: topological generation levels.

    Headwater catchments (nobody drains into them... precisely, catchments
    that are not the downstream target of any other catchment) get level 1;
    every other catchment gets ``1 + max(level of its direct upstreams)``.
    Note this is *not* the classic Shreve magnitude (which sums upstream
    orders); it is the level definition GREENeR's ``aggregate_loop`` relies
    on for sequential evaluation.
    """
    df = catch_data.copy()
    hydro = df["HydroID"].to_numpy()
    to = df["To_catch"].to_numpy()
    ids = set(hydro.tolist())
    upstreams: dict = {h: [] for h in ids}
    for h, t in zip(hydro, to):
        if t in upstreams:
            upstreams[t].append(h)
    level: dict = {}
    remaining = set(ids)
    while remaining:
        progressed = False
        for h in list(remaining):
            ups = upstreams[h]
            if not ups:
                level[h] = 1
            elif all(u in level for u in ups):
                # R loop: max Shreve of the feeding catchments, + 1.
                level[h] = max(level[u] for u in ups) + 1
            else:
                continue
            remaining.discard(h)
            progressed = True
        if not progressed:
            raise ValueError("Topology contains a cycle or unreachable "
                             "catchments; cannot assign Shreve levels.")
    df["Shreve"] = df["HydroID"].map(level).astype(int)
    return df


def _topological_order(catch_data: pd.DataFrame) -> tuple[np.ndarray, dict]:
    """Return catchment indices in an order where every catchment follows
    all of its direct upstreams (Kahn's algorithm), plus the upstream
    index lists. Equivalent to GREENeR's Shreve-grouped loop for valid
    inputs, and robust to whatever Shreve convention the caller used.
    """
    hydro = catch_data["HydroID"].to_numpy()
    to = catch_data["To_catch"].to_numpy()
    pos = {h: i for i, h in enumerate(hydro)}
    upstreams: dict[int, list[int]] = {i: [] for i in range(len(hydro))}
    for i, t in enumerate(to):
        if t in pos:
            upstreams[pos[t]].append(i)
    indeg = {i: len(upstreams[i]) for i in range(len(hydro))}
    # downstream pointer per node (index or -1)
    down = np.array([pos.get(t, -1) for t in to], dtype=int)
    queue = [i for i in range(len(hydro)) if indeg[i] == 0]
    order: list[int] = []
    remaining = dict(indeg)
    while queue:
        i = queue.pop()
        order.append(i)
        d = down[i]
        if d >= 0:
            remaining[d] -= 1
            if remaining[d] == 0:
                queue.append(d)
    if len(order) != len(hydro):
        raise ValueError("Topology contains a cycle; routing impossible.")
    return np.array(order, dtype=int), upstreams


# --------------------------------------------------------------------------
# Input preparation (mirrors check_colnames_* / append_empty_cols)
# --------------------------------------------------------------------------

def _prepare(catch_data: pd.DataFrame, annual_data: pd.DataFrame):
    missing = [c for c in _CATCH_COLS if c not in catch_data.columns]
    if missing:
        raise ValueError(f"catch_data is missing columns: {missing}")
    catch = catch_data.copy()
    if "Shreve" not in catch.columns:
        catch = compute_shreve_levels(catch)
    catch = catch[["HydroID", "To_catch", "Shreve", "LakeFrRet",
                   "NrmLengthKm"]].copy()

    annual = annual_data.copy()
    if all(c in annual.columns for c in _TP_COLS) and "Atm" not in annual.columns:
        # TP scenario: GREENeR renames Bg -> Atm and appends Fix = Soil = 0.
        annual = annual[_TP_COLS].rename(columns={"Bg": "Atm"}).copy()
        annual["Fix"] = 0.0
        annual["Soil"] = 0.0
        nutrient = "TP"
    elif all(c in annual.columns for c in _TN_COLS):
        annual = annual[_TN_COLS].copy()
        nutrient = "TN"
    else:
        raise ValueError("annual_data columns match neither the GREENeR "
                         "TN nor the TP schema.")
    return catch, annual, nutrient


def _check_years(loc_years, available):
    if loc_years is None:
        raise ValueError("Years cannot be null (GREENeR check_years).")
    years = [y for y in loc_years if y in set(available)]
    if len(years) != len(list(loc_years)):
        import warnings
        warnings.warn("Some years are not in the scenario.")
    return years


# --------------------------------------------------------------------------
# Core routing (mirrors aggregate_loop)
# --------------------------------------------------------------------------

def _route(catch: pd.DataFrame, per_catch_year: pd.DataFrame,
           years: list[int], balance: bool = False) -> pd.DataFrame:
    """Accumulate CatchToRiver downstream.

    per_catch_year must carry HydroID, Year, CatchToRiver (and, for the
    balance variant, the extra flux columns of region_nut_balance).
    """
    order, upstreams = _topological_order(catch)
    pos = {h: i for i, h in enumerate(catch["HydroID"].to_numpy())}
    lake = catch["LakeFrRet"].to_numpy(dtype=float)
    bb = catch["BB"].to_numpy(dtype=float)
    out_frames = []
    for yr in years:
        sub = per_catch_year[per_catch_year["Year"] == yr]
        ctr = np.zeros(len(catch))
        ctr[sub["HydroID"].map(pos).to_numpy()] = \
            sub["CatchToRiver"].to_numpy(dtype=float)
        load = np.zeros(len(catch))
        riv_ret = np.zeros(len(catch))
        lake_ret = np.zeros(len(catch))
        for i in order:
            up = load[upstreams[i]].sum() if upstreams[i] else 0.0
            load[i] = (1.0 - lake[i]) * (ctr[i] + up) * bb[i]
            if balance:
                # Faithful to the R sequential := updates in
                # aggregate_loop(type="N"): CatchRivRet and CatchLakeRet
                # are computed from the ALREADY UPDATED CatchLoad.
                riv_ret[i] = (ctr[i] + load[i]) * (1.0 - bb[i])
                lake_ret[i] = lake[i] * (ctr[i] + load[i]) * bb[i]
        frame = catch.copy()
        frame["Year"] = yr
        frame["CatchToRiver"] = ctr
        frame["CatchLoad"] = load
        if balance:
            frame["CatchRivRet"] = riv_ret
            frame["CatchLakeRet"] = lake_ret
            extra = sub.drop(columns=["CatchToRiver"], errors="ignore")
            frame = frame.merge(extra, left_on=["HydroID", "Year"],
                                right_on=["HydroID", "Year"], how="left")
        else:
            obs = sub.set_index("HydroID")["ObsLoad"]
            frame["ObsLoad"] = frame["HydroID"].map(obs)
        out_frames.append(frame)
    return pd.concat(out_frames, ignore_index=True)


def _local_fluxes(annual: pd.DataFrame, alpha_p: float, sd_coef: float,
                  atm_coeff: float) -> pd.DataFrame:
    """The per-catchment, per-year local terms of launch_green."""
    a = annual.copy()
    a["AA"] = np.exp(-alpha_p * a["InvNrmRain"].to_numpy(dtype=float))
    a["SumDS"] = (a["Min"] + a["Man"] + a["Fix"] + a["Soil"]
                  + a["Atm"] * (1.0 - a["ForestFraction"]))
    a["SumPSCmp1"] = (a["Ps"] + atm_coeff * a["Atm"] * a["ForestFraction"]
                      + sd_coef * a["Sd"])
    a["CatchToRiver"] = a["SumDS"] * a["AA"] + a["SumPSCmp1"]
    return a


# --------------------------------------------------------------------------
# Public model functions
# --------------------------------------------------------------------------

def launch_green(catch_data: pd.DataFrame, annual_data: pd.DataFrame,
                 alpha_p: float, alpha_l: float, sd_coef: float,
                 loc_years, atm_coeff: float = 0.38) -> pd.DataFrame:
    """Port of GREENeR's ``launch_green``.

    Returns one row per catchment and year with the GREENeR output
    columns (HydroID, To_catch, Shreve, LakeFrRet, NrmLengthKm, BB,
    Year, ObsLoad, CatchToRiver, CatchLoad). BB and AA are included for
    transparency; GREENeR carries BB in its output table as well.
    """
    catch, annual, _ = _prepare(catch_data, annual_data)
    years = _check_years(loc_years, annual["YearValue"].unique())
    catch = catch.copy()
    catch["BB"] = np.exp(-alpha_l * catch["NrmLengthKm"].to_numpy(dtype=float))
    a = _local_fluxes(annual, alpha_p, sd_coef, atm_coeff)
    per = a.rename(columns={"YearValue": "Year",
                            "YearlyMass": "ObsLoad"})[
        ["HydroID", "Year", "ObsLoad", "CatchToRiver"]]
    result = _route(catch, per, years)
    aa = a.set_index(["HydroID", "YearValue"])["AA"]
    result["AA"] = [aa.get((h, y), np.nan)
                    for h, y in zip(result["HydroID"], result["Year"])]
    cols = ["HydroID", "To_catch", "Shreve", "LakeFrRet", "NrmLengthKm",
            "BB", "AA", "Year", "ObsLoad", "CatchToRiver", "CatchLoad"]
    return result[cols].sort_values(["Year", "HydroID"]).reset_index(drop=True)


def green_shares(catch_data: pd.DataFrame, annual_data: pd.DataFrame,
                 alpha_p: float, alpha_l: float, sd_coef: float,
                 loc_years, atm_coeff: float = 0.38) -> pd.DataFrame:
    """Port of GREENeR's ``green_shares``: source apportionment.

    Runs the model once per source with all other source inputs zeroed
    (valid because the GREEN equations are linear in the sources for
    fixed retention factors) and returns, per catchment and year, the
    routed load attributable to each source plus their sum, CatchLoad.
    """
    nutrient = "TP" if ("Bg" in annual_data.columns
                        and "Atm" not in annual_data.columns) else "TN"
    sources = TN_SOURCES if nutrient == "TN" else TP_SOURCES
    frames = []
    for src in sources:
        mod = annual_data.copy()
        for other in sources:
            if other != src:
                mod[other] = 0.0
        res = launch_green(catch_data, mod, alpha_p, alpha_l, sd_coef,
                           loc_years, atm_coeff=atm_coeff)
        res = res[["HydroID", "To_catch", "Year", "CatchLoad"]].rename(
            columns={"CatchLoad": src})
        frames.append(res.set_index(["HydroID", "To_catch", "Year"]))
    out = pd.concat(frames, axis=1).reset_index()
    out["CatchLoad"] = out[sources].sum(axis=1)
    return out.sort_values(["Year", "HydroID"]).reset_index(drop=True)


def region_nut_balance(catch_data: pd.DataFrame, annual_data: pd.DataFrame,
                       alpha_p: float, alpha_l: float, sd_coef: float,
                       loc_years, atm_coeff: float = 0.38) -> pd.DataFrame:
    """Port of GREENeR's ``region_nut_balance``.

    Adds the diffuse / land-retention decomposition columns and the
    routed CatchRivRet / CatchLakeRet terms. NOTE: following the R
    code's sequential data.table updates, CatchRivRet and CatchLakeRet
    are computed from the already-updated CatchLoad (see PORT_NOTES.md).
    """
    catch, annual, _ = _prepare(catch_data, annual_data)
    years = _check_years(loc_years, annual["YearValue"].unique())
    catch = catch.copy()
    catch["BB"] = np.exp(-alpha_l * catch["NrmLengthKm"].to_numpy(dtype=float))
    a = _local_fluxes(annual, alpha_p, sd_coef, atm_coeff)
    a["Atm2DiffAgri"] = a["Atm"] * (1.0 - a["ForestFraction"])
    a["Atm2Diff"] = atm_coeff * a["Atm"] * a["ForestFraction"]
    a["SD2Diff"] = sd_coef * a["Sd"]
    a["Diffuse"] = a["SumDS"] * a["AA"] + a["Atm2Diff"] + a["SD2Diff"]
    a["Atm2LandRet"] = (1.0 - atm_coeff) * a["Atm"] * a["ForestFraction"]
    a["SD2LandRet"] = (1.0 - sd_coef) * a["Sd"]
    a["DiffAgri2LandRet"] = a["SumDS"] * (1.0 - a["AA"])
    a["LandRet"] = (a["DiffAgri2LandRet"] + a["Atm2LandRet"]
                    + a["SD2LandRet"])
    per = a.rename(columns={"YearValue": "Year"})[[
        "HydroID", "Year", "Atm", "Sd", "Min", "Man", "Fix", "Soil", "Ps",
        "SumDS", "Atm2DiffAgri", "Atm2Diff", "SD2Diff", "Diffuse",
        "DiffAgri2LandRet", "Atm2LandRet", "SD2LandRet", "LandRet",
        "CatchToRiver"]].rename(columns={"SumDS": "DiffuseAgri",
                                         "Ps": "PS"})
    result = _route(catch, per, years, balance=True)
    return result.sort_values(["Year", "HydroID"]).reset_index(drop=True)


# --------------------------------------------------------------------------
# GREENeR_data CSV reader (port of read_NSdata)
# --------------------------------------------------------------------------

def read_nsdata(path: str):
    """Port of GREENeR's ``read_NSdata`` for a folder of CSV files in
    the GREENeR_data layout (TS_nutrients.csv, Obs_monitoring.csv,
    ForestFr.csv, Precipitation.csv, Topology.csv, LakeProperties.csv,
    Length.csv). Returns (catch_data, annual_data).
    """
    import os
    def rd(name):
        return pd.read_csv(os.path.join(path, name))
    tsn = rd("TS_nutrients.csv")
    nutrient = "TN" if tsn.shape[1] == 9 else "TP"
    topo = rd("Topology.csv").rename(columns={"Next_HydroID": "To_catch"})
    lak = rd("LakeProperties.csv")
    length = rd("Length.csv")
    lake_r_par = 7.3 if nutrient == "TN" else 26.0
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = lake_r_par / lak["AvgDepth_m"].to_numpy(dtype=float) \
            * lak["ResTime"].to_numpy(dtype=float)
        fr = 1.0 - 1.0 / (1.0 + ratio)
    fr = np.where(np.isnan(fr), 0.0, fr)
    fr = np.where(fr > 0.1, 0.1, fr)
    lak = lak.assign(LakeFrRet=fr)[["HydroID", "LakeFrRet"]]
    length = length.assign(
        NrmLengthKm=length["LengthKm"] / length["LengthKm"].max())[
        ["HydroID", "NrmLengthKm"]]
    topo = compute_shreve_levels(topo.rename(
        columns={"To_catch": "To_catch"}))
    catch = topo.merge(lak, on="HydroID").merge(length, on="HydroID")
    catch = catch[["HydroID", "To_catch", "Shreve", "LakeFrRet",
                   "NrmLengthKm"]]
    obs = rd("Obs_monitoring.csv")
    ff = rd("ForestFr.csv")
    pre = rd("Precipitation.csv")
    min_rain = 50.0
    pre = pre.assign(InvNrmRain=min_rain / np.maximum(min_rain, pre["Rain"]))
    scen = tsn.merge(obs, on=["YearValue", "HydroID"], how="left")
    scen = scen.merge(ff, on=["YearValue", "HydroID"], how="left")
    scen = scen.merge(pre[["YearValue", "HydroID", "InvNrmRain"]],
                      on=["YearValue", "HydroID"], how="left")
    scen = scen.merge(topo[["HydroID", "To_catch"]], on="HydroID")
    scen = scen.rename(columns={"To_catch": "NextDownID"})
    scen["BasinID"] = 1234
    return catch, scen
