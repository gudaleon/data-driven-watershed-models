"""Synthetic recovery test for the Bayesian SPARROW (sparrow.stan).

This is the build's verification gate, run BEFORE any real-data fit:

1. Take one real sub-basin of the extended network (all reaches whose
   ``outlet_sparrowid`` is the chosen basin's outlet).
2. Fix known parameter values (TRUE below), generate the reach loads
   with a numpy mirror of the Stan accumulation, and cross-check that
   mirror against the verified classical port (ClassicalSparrow with
   an equivalent spec, unconditioned prediction) -- the mirror must
   match the port to 1e-8, so the data-generating code is the verified
   code, not a second unchecked implementation.
3. Simulate station observations (the basin's real load stations plus
   additional synthetic sites) in log space with the model's own
   observation error.
4. Fit with the compiled Stan model (CmdStan CLI, no cmdstanpy) and
   check that the true values fall inside the posterior 90% intervals.
5. Fit the same synthetic data with a single-block random-walk
   Metropolis sampler (the Gibbs-style reference; see BUILD_NOTES.md
   for why this, and not a BUGS run, is the honest comparison) and
   report the posterior SD of ks (the in-stream attenuation
   coefficient) from both samplers and their ratio.  Alex's earlier
   MultiBUGS-to-Stan port found ks "way narrow" in Stan; the reported
   quantity for the real fit is this same Stan/Metropolis SD ratio.

Outputs: synthetic/synthetic_results.json and the Stan draw CSVs
beside it.  Nothing here touches the real loads.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import prepare_stan_data as prep  # noqa: E402
from classical_sparrow import (  # noqa: E402
    ClassicalSparrow,
    Parameter,
    SparrowSpec,
)

OUT = os.path.join(HERE, "synthetic")
os.makedirs(OUT, exist_ok=True)
SEED = int(os.environ.get("SYNTH_SEED", "20261121"))
TAG = "" if SEED == 20261121 else f"_s{SEED}"
RNG = np.random.default_rng(SEED)

TRUE = dict(
    alpha_fert=0.35, alpha_man=0.25, rate_urban=280.0, rate_bg=22.0,
    rate_dep=480.0, beta_runoff=0.25, beta_perm=-0.20, beta_tile=0.15,
    ks=0.18, beta_res=0.45, sigma=0.15,
)
OBS_SE = 0.08
PARAMS = list(TRUE.keys())


# ----------------------------------------------------------------------
# numpy mirror of the Stan accumulation
# ----------------------------------------------------------------------
def simulate(arrays: dict, p: dict) -> np.ndarray:
    n = arrays["N"]
    deliv = np.exp(
        p["beta_runoff"] * arrays["z_runoff"]
        + p["beta_perm"] * arrays["z_perm"]
        + p["beta_tile"] * arrays["z_tile"]
    )
    G = (
        arrays["src_point"]
        + deliv * (p["alpha_fert"] * arrays["src_fert"]
                   + p["alpha_man"] * arrays["src_man"])
        + p["rate_urban"] * arrays["area_urban"]
        + p["rate_bg"] * arrays["area_bg"]
        + p["rate_dep"] * arrays["area_water"]
    )
    is_res = arrays["is_res"] == 1
    A = np.where(is_res, 1.0, np.exp(-p["ks"] * arrays["ttime"]))
    B = np.where(is_res, 1.0 / (1.0 + p["beta_res"] * arrays["rhload"]), 1.0)
    inc = np.where(is_res, B, np.sqrt(A)) * G
    carry = arrays["frac"] * np.where(is_res, B, A)
    inflow = np.zeros(n + 1)
    mu = np.empty(n)
    d1, d2 = arrays["downstream1"], arrays["downstream2"]
    for r in range(n):
        tot = inc[r] + carry[r] * inflow[r]
        mu[r] = tot
        inflow[d1[r] - 1] += tot
        if d2[r] != d1[r]:
            inflow[d2[r] - 1] += tot
    return mu


def cross_check_with_classical(frame: pd.DataFrame, arrays: dict,
                               mu_mirror: np.ndarray) -> float:
    """Reproduce the mirror's loads with the verified classical port."""
    d = frame.copy()
    deliv = np.exp(
        TRUE["beta_runoff"] * arrays["z_runoff"]
        + TRUE["beta_perm"] * arrays["z_perm"]
        + TRUE["beta_tile"] * arrays["z_tile"]
    )
    # fold the shared delivery factor and the source coefficients into
    # an equivalent classical spec: sources carry coefficient 1 and the
    # Stan source terms are pre-multiplied.
    d["s_point"] = arrays["src_point"]
    d["s_fert"] = deliv * arrays["src_fert"]
    d["s_man"] = deliv * arrays["src_man"]
    d["s_urb"] = arrays["area_urban"]
    d["s_bg"] = arrays["area_bg"]
    d["s_wat"] = arrays["area_water"]
    d["dec_t"] = np.where(arrays["is_res"] == 1, 0.0, arrays["ttime"])
    d["res_r"] = np.where(arrays["is_res"] == 1, arrays["rhload"], 0.0)
    d["depvar"] = 0.0
    d["iftran"] = 1.0  # data1 carries no iftran column; transfer is whole
    # classical_sparrow works on the frame in hydrologic order internally
    # but expects the ORIGINAL (unordered) frame; our frame is already in
    # topological order, and its own topological sort is a valid order
    # too -- predictions are returned in input row order either way.
    spec = SparrowSpec(
        sources=["s_point", "s_fert", "s_man", "s_urb", "s_bg", "s_wat"],
        stream_decay_vars=["dec_t"],
        reservoir_decay_vars=["res_r"],
        mean_adjust_delivery=False,
    )
    model = ClassicalSparrow(d, spec)
    beta = np.array([1.0, TRUE["alpha_fert"], TRUE["alpha_man"],
                     TRUE["rate_urban"], TRUE["rate_bg"], TRUE["rate_dep"],
                     TRUE["ks"], TRUE["beta_res"]])
    pred = model.predict_loads(beta, condition=False)
    # `pred` is in the frame's row order, which is topological order --
    # the same order as mu_mirror.
    denom = np.maximum(np.abs(mu_mirror), 1e-12)
    return float(np.max(np.abs(pred - mu_mirror) / denom))


# ----------------------------------------------------------------------
# Metropolis reference sampler (same posterior as the Stan model)
# ----------------------------------------------------------------------
LOGP_PRIOR = dict(  # (median, log-sd) for lognormal priors
    alpha_fert=(0.3, 0.5), alpha_man=(0.3, 0.5), rate_urban=(300.0, 0.5),
    rate_bg=(25.0, 0.5), rate_dep=(500.0, 0.25), ks=(0.15, 0.5),
    beta_res=(0.5, 0.75),
)
POSITIVE = ["alpha_fert", "alpha_man", "rate_urban", "rate_bg",
            "rate_dep", "ks", "beta_res", "sigma"]


def _lists(arrays: dict) -> dict:
    """Plain-python-list copy of the arrays for the fast MH simulator."""
    return {k: (v.tolist() if isinstance(v, np.ndarray) else v)
            for k, v in arrays.items()}


def simulate_fast(L: dict, p: dict) -> list:
    """Pure-python twin of simulate() (same arithmetic, no numpy
    scalars) so the Metropolis reference can afford ~10^4 evaluations."""
    import math

    n = L["N"]
    inflow = [0.0] * (n + 1)
    mu = [0.0] * n
    for r in range(n):
        deliv = math.exp(p["beta_runoff"] * L["z_runoff"][r]
                         + p["beta_perm"] * L["z_perm"][r]
                         + p["beta_tile"] * L["z_tile"][r])
        G = (L["src_point"][r]
             + deliv * (p["alpha_fert"] * L["src_fert"][r]
                        + p["alpha_man"] * L["src_man"][r])
             + p["rate_urban"] * L["area_urban"][r]
             + p["rate_bg"] * L["area_bg"][r]
             + p["rate_dep"] * L["area_water"][r])
        if L["is_res"][r] == 1:
            B = 1.0 / (1.0 + p["beta_res"] * L["rhload"][r])
            tot = B * G + L["frac"][r] * B * inflow[r]
        else:
            A = math.exp(-p["ks"] * L["ttime"][r])
            tot = math.sqrt(A) * G + L["frac"][r] * A * inflow[r]
        mu[r] = tot
        inflow[L["downstream1"][r] - 1] += tot
        if L["downstream2"][r] != L["downstream1"][r]:
            inflow[L["downstream2"][r] - 1] += tot
    return mu


def log_posterior(theta: np.ndarray, L: dict, obs_reach: np.ndarray,
                  log_obs: np.ndarray, obs_se: np.ndarray) -> float:
    p = {}
    for i, name in enumerate(PARAMS):
        p[name] = np.exp(theta[i]) if name in POSITIVE else theta[i]
    lp = 0.0
    for name, (med, sd) in LOGP_PRIOR.items():
        lp += -0.5 * ((np.log(p[name]) - np.log(med)) / sd) ** 2  # jacobian cancels -log x
    for name in ["beta_runoff", "beta_perm", "beta_tile"]:
        lp += -0.5 * (p[name] / 0.5) ** 2
    lp += -0.5 * (p["sigma"] / 0.5) ** 2 + np.log(p["sigma"])  # half-normal + jacobian
    mu = simulate_fast(L, p)
    sd_o = np.sqrt(obs_se ** 2 + p["sigma"] ** 2)
    ll = 0.0
    for o, r in enumerate(obs_reach):
        if mu[r] <= 0:
            return -np.inf
        ll += -0.5 * ((log_obs[o] - np.log(mu[r])) / sd_o[o]) ** 2 - np.log(sd_o[o])
    return lp + float(ll)


def metropolis(arrays, obs_reach, log_obs, obs_se,
               n_iter=12000, burn=4000, seed=7):
    """Single-block random-walk Metropolis on the transformed
    parameters (log for positive parameters), same posterior as Stan."""
    rng = np.random.default_rng(seed)
    L = _lists(arrays)
    theta = np.array([np.log(TRUE[k]) if k in POSITIVE else TRUE[k]
                      for k in PARAMS], dtype=float)
    theta += rng.normal(0, 0.3, size=len(PARAMS))  # start off-truth
    step = np.full(len(PARAMS), 0.045)
    cur = log_posterior(theta, L, obs_reach, log_obs, obs_se)
    draws = []
    acc = 0
    for it in range(n_iter):
        prop = theta + rng.normal(0, step)
        val = log_posterior(prop, L, obs_reach, log_obs, obs_se)
        if np.log(rng.random()) < val - cur:
            theta, cur = prop, val
            acc += 1
        if it == burn - 1:  # one-step adaptation, then fixed
            rate = acc / burn
            step = step * (1.5 if rate > 0.35 else 0.55 if rate < 0.15 else 1.0)
            acc = 0
        if it >= burn:
            draws.append(theta.copy())
    draws = np.array(draws)
    out = {}
    for j, name in enumerate(PARAMS):
        v = np.exp(draws[:, j]) if name in POSITIVE else draws[:, j]
        out[name] = v
    out["_accept_rate_post_burn"] = acc / max(n_iter - burn, 1)
    return out


# ----------------------------------------------------------------------
def main() -> None:
    full = prep.load_network()
    scale = prep.rhload_scale(full)
    basins = pd.read_csv(os.path.join(prep.SPARROW_DIR, "extended_basins_v1.csv"))
    loads = pd.read_csv(os.path.join(prep.SPARROW_DIR, "station_year_loads_v4.csv"),
                        dtype={"station_id": str})
    on_reaches = set(loads.loc[loads["on_network"], "reach_sparrowid"])
    cand = basins[(basins["n_reaches"].between(150, 600))].copy()
    cand["n_stations"] = cand["outlet_sparrowid"].map(
        full.groupby("outlet_sparrowid")["reach_sparrowid"]
        .apply(lambda s: int(s.isin(on_reaches).sum())))
    cand = cand[cand["n_stations"] >= 6].sort_values("n_reaches")
    if cand.empty:
        raise RuntimeError("no sub-basin met the size/station screen")
    basin = cand.iloc[0]
    outlet_id = int(basin["outlet_sparrowid"])
    print(f"synthetic basin: outlet {outlet_id}, {int(basin['n_reaches'])} "
          f"reaches, {int(basin['n_stations'])} real load stations")
    sub = full[full["outlet_sparrowid"] == outlet_id].copy()
    info = prep.build_arrays(sub, scale=scale)
    arrays, frame = info["arrays"], info["frame"]
    mu = simulate(arrays, TRUE)
    rel_err = cross_check_with_classical(frame, arrays, mu)
    print(f"mirror vs classical port: max relative difference {rel_err:.3g}")
    if rel_err > 1e-8:
        raise RuntimeError("numpy mirror disagrees with the classical port")
    fast = np.array(simulate_fast(_lists(arrays), TRUE))
    fast_err = float(np.max(np.abs(fast - mu) / np.maximum(np.abs(mu), 1e-12)))
    print(f"fast MH simulator vs mirror: max relative difference {fast_err:.3g}")
    if fast_err > 1e-10:
        raise RuntimeError("fast Metropolis simulator disagrees with the mirror")

    # observation sites: the basin's real stations + synthetic sites on
    # larger reaches, so the recovery test is identified
    site_pos = sorted({info["pos_of_reach"][r] for r in
                       set(frame.loc[frame["reach_sparrowid"].isin(on_reaches),
                                     "reach_sparrowid"])})
    rng = RNG
    big = np.where((frame["demtarea"].to_numpy() > 15) & (mu > 0))[0] + 1
    extra = rng.choice(big, size=min(25, len(big)), replace=False).tolist()
    site_pos = sorted(set(site_pos) | set(extra))
    obs_reach = np.array(site_pos, dtype=int) - 1  # 0-based for numpy
    sd_o = np.sqrt(OBS_SE ** 2 + TRUE["sigma"] ** 2)
    log_obs = np.log(mu[obs_reach]) + rng.normal(0, sd_o, size=len(obs_reach))
    obs = pd.DataFrame(dict(obs_reach=site_pos, log_obs=log_obs,
                            obs_se=np.full(len(site_pos), OBS_SE)))
    data = prep.to_stan_json(arrays, obs, [info["pos_of_reach"][outlet_id]])
    data_path = os.path.join(OUT, f"synthetic_data{TAG}.json")
    with open(data_path, "w") as f:
        json.dump(data, f)
    print(f"synthetic observations: {len(obs)} sites")

    # ---- Stan fit ------------------------------------------------------
    exe = os.path.join(HERE, "sparrow")
    out_base = os.path.join(OUT, f"draws{TAG}")
    subprocess.run(
        [exe, "sample", "num_chains=4", "num_warmup=600", "num_samples=600",
         "adapt", "delta=0.9", "data", f"file={data_path}",
         "output", f"file={out_base}.csv", "random", f"seed={SEED}"],
        check=True, capture_output=True, text=True)
    draws = pd.concat([pd.read_csv(f"{out_base}_{c}.csv", comment="#")
                       for c in range(1, 5)], ignore_index=True)
    stan_summary = {}
    for name in PARAMS:
        v = draws[name].to_numpy()
        stan_summary[name] = dict(
            true=TRUE[name], mean=float(v.mean()), sd=float(v.std(ddof=1)),
            q05=float(np.quantile(v, 0.05)), q95=float(np.quantile(v, 0.95)),
            covered=bool(np.quantile(v, 0.05) <= TRUE[name]
                         <= np.quantile(v, 0.95)))
    n_div = int((draws["divergent__"] == 1).sum())
    max_rhat_proxy = None

    # ---- Metropolis reference ------------------------------------------
    mh = metropolis(arrays, obs_reach, log_obs,
                    np.full(len(obs_reach), OBS_SE), seed=SEED % 1000)
    mh_acc = mh.pop("_accept_rate_post_burn")
    mh_summary = {name: dict(mean=float(v.mean()), sd=float(v.std(ddof=1)))
                  for name, v in mh.items()}

    ks_ratio = stan_summary["ks"]["sd"] / mh_summary["ks"]["sd"]
    results = dict(
        basin_outlet=outlet_id, n_reaches=int(arrays["N"]),
        n_obs=int(len(obs)), mirror_vs_classical_max_rel=rel_err,
        stan=stan_summary, metropolis=mh_summary,
        metropolis_accept_rate_post_burn=float(mh_acc),
        ks_sd_stan=stan_summary["ks"]["sd"], ks_sd_metropolis=mh_summary["ks"]["sd"],
        ks_sd_ratio_stan_over_metropolis=float(ks_ratio),
        stan_divergences=n_div,
        n_covered=int(sum(s["covered"] for s in stan_summary.values())),
        n_params=len(PARAMS),
    )
    with open(os.path.join(OUT, f"synthetic_results{TAG}.json"), "w") as f:
        json.dump(results, f, indent=2)
    print(json.dumps({k: results[k] for k in
                      ["n_covered", "n_params", "stan_divergences",
                       "ks_sd_stan", "ks_sd_metropolis",
                       "ks_sd_ratio_stan_over_metropolis"]}, indent=2))
    for name in PARAMS:
        s = stan_summary[name]
        print(f"  {name:11s} true {s['true']:10.4g}  post mean {s['mean']:10.4g}"
              f"  sd {s['sd']:9.4g}  90% CI [{s['q05']:.4g}, {s['q95']:.4g}]"
              f"  covered={s['covered']}")


if __name__ == "__main__":
    main()
