"""Verification gate against RSPARROW's published tutorial outputs.

Uses the bundled UserTutorial data (MRB3 total nitrogen) and tutorial
Model8 -- the fullest tutorial model: 13 estimated parameters
(5 SOURCE, 5 DELIVF, 2 STRM, 1 RESV), monitoring-adjusted SSE
115.6632 on 708 calibration sites (Model8_summary.txt).

Checks, all against files RSPARROW itself wrote:
  1. unconditioned per-reach loads vs predict/Model8_predicts_load.csv
     (`pload_total`), using the published NLLS estimates;
  2. conditioned per-reach loads vs `mpload_total`;
  3. per-station log residuals vs estimate/Model8_residuals.csv and
     the SSE vs the summary's 115.6632;
  4. an NLLS refit from the parameters.csv starting values, compared
     with the published estimates.

Run:  ~/workspace/mapvenv/bin/python verify_classical_sparrow.py
"""

import sys

import numpy as np
import pandas as pd

sys.path.insert(0, ".")
from classical_sparrow import ClassicalSparrow, Parameter, SparrowSpec  # noqa: E402

TUT = "/home/hatch/workspace/fei_lake_ontario/codebase/RSPARROW/UserTutorial"
M8 = f"{TUT}/results/Model8"

# Published NLLS estimates, Model8_summary.txt (PARAMETER ESTIMATES).
PUBLISHED = {
    "point": 0.789478293, "ndep": 0.513053881, "MANC_N": 0.290962380,
    "FARM_N": 0.130656653, "Fixation": 6.251067012,
    "ldrainden": 0.134403345, "PPT30MEAN": 0.001592016,
    "meanTemp": -0.041235668, "tiles_perc": 1.132953147,
    "soil_CLAYAVE": 0.013925452,
    "rchdecay1": 0.424093222, "rchdecay2": 0.233019423,
    "iresload": 6.710379414,
}
PUBLISHED_SSE = 115.6632  # monitoring-adjusted, Model8_summary.txt

SOURCES = ["point", "ndep", "MANC_N", "FARM_N", "Fixation"]
DELIV = ["ldrainden", "PPT30MEAN", "meanTemp", "tiles_perc", "soil_CLAYAVE"]
STRM = ["rchdecay1", "rchdecay2"]
RESV = ["iresload"]
BETA_ORDER = SOURCES + DELIV + STRM + RESV

# ----------------------------------------------------------------------
# Build the reach table the way RSPARROW's data prep does:
# depvar = observed load at calibration sites only
# (sitedata = subdata[depvar > 0 & calsites == 1], calsites = Tagsite
# per Model8_userModifyData.R), rows in hydseq order (handled by the
# model class).
# ----------------------------------------------------------------------
d1 = pd.read_csv(f"{TUT}/data/data1.csv")
n_sta_load = int(((d1["staid"] > 0) & (d1["LOAD_A_00600"] > 0)).sum())
d1["depvar"] = np.where(
    (d1["staid"] > 0) & (d1["Tagsite"] == 1), d1["LOAD_A_00600"], 0.0
)
d1["depvar"] = d1["depvar"].fillna(0.0)
n_cal = int((d1["depvar"] > 0).sum())
print(f"reaches: {len(d1)}; reaches with staid>0 and load>0: {n_sta_load}; "
      f"calibration sites (Tagsite==1): {n_cal}")

design_df = pd.read_csv(f"{M8}/Model8_design_matrix.csv", index_col=0)
design = design_df.loc[SOURCES, DELIV].to_numpy(dtype=float)

spec = SparrowSpec(
    sources=SOURCES,
    delivery_vars=DELIV,
    design=design,
    stream_decay_vars=STRM,
    reservoir_decay_vars=RESV,
    mean_adjust_delivery=True,
)
model = ClassicalSparrow(d1, spec)
beta_pub = np.array([PUBLISHED[n] for n in BETA_ORDER])

# ----------------------------------------------------------------------
# 1-2. Per-reach loads vs the published prediction file
# ----------------------------------------------------------------------
pred_pub = pd.read_csv(f"{M8}/predict/Model8_predicts_load.csv")
mine = pd.DataFrame(
    {
        "mrb_id": d1["mrb_id"].to_numpy(),
        "pload_py": model.predict_loads(beta_pub, condition=False),
        "mpload_py": model.predict_loads(beta_pub, condition=True),
    }
).merge(pred_pub[["mrb_id", "pload_total", "mpload_total"]], on="mrb_id")
for col_py, col_pub in [("pload_py", "pload_total"), ("mpload_py", "mpload_total")]:
    diff = (mine[col_py] - mine[col_pub]).abs()
    rel = diff / mine[col_pub].abs().clip(lower=1.0)
    print(f"{col_pub}: n={len(mine)}, max abs diff={diff.max():.6g}, "
          f"max rel diff (loads>=1)={rel.max():.3e}, "
          f"corr={np.corrcoef(mine[col_py], mine[col_pub])[0, 1]:.12f}")

# RSPARROW's prediction files multiply every load by a single
# retransformation factor, bootcorrection = mean_exp_weighted_error
# (the summary's "Mean Exp Weighted Error"; predict.R line 413,
# controlFileTasksModel.R line 276, estimateNLLSmetrics.R line 313):
#     mean(exp(boot_resid)),  boot_resid = Resids / sqrt(1 - leverage)
# with leverage from the NLLS Jacobian.  Verify that the published
# file is exactly the port's accumulation times that factor, with the
# factor recomputed from the port's own residuals and Jacobian.
from scipy.optimize._numdiff import approx_derivative  # noqa: E402

ratio = mine["pload_total"] / mine["pload_py"]
resid_pub_beta = model.residuals(beta_pub)
jac = approx_derivative(model.residuals, beta_pub, method="2-point")
hat = jac @ np.linalg.solve(jac.T @ jac, jac.T)
leverage = np.diag(hat)
boot_resid = resid_pub_beta / np.sqrt(1.0 - leverage)
bootcorr_py = float(np.mean(np.exp(boot_resid)))
print(f"published/raw ratio: min={ratio.min():.8f} max={ratio.max():.8f}")
print(f"bootcorrection from Python residuals+Jacobian: {bootcorr_py:.8f}")
for col_py, col_pub in [("pload_py", "pload_total"), ("mpload_py", "mpload_total")]:
    corrected = mine[col_py] * bootcorr_py
    rel = (corrected - mine[col_pub]).abs() / mine[col_pub].abs().clip(lower=1.0)
    print(f"{col_pub} after bootcorrection: max rel diff={rel.max():.3e}")

# ----------------------------------------------------------------------
# 3. Residuals and SSE vs the published estimate file
# ----------------------------------------------------------------------
resid_pub = pd.read_csv(f"{M8}/estimate/Model8_residuals.csv")
site_rows = model.data.loc[model.is_site, "mrb_id"].to_numpy()
mine_resid = pd.DataFrame({"mrb_id": site_rows, "resid_py": model.residuals(beta_pub)})
both = mine_resid.merge(
    resid_pub[["waterid", "Resids"]], left_on="mrb_id", right_on="waterid"
)
rdiff = (both["resid_py"] - both["Resids"]).abs()
print(f"residuals: n={len(both)}, max abs diff vs published={rdiff.max():.3e}")
sse_py = model.sse(beta_pub)
print(f"SSE at published estimates: Python={sse_py:.4f}  "
      f"RSPARROW={PUBLISHED_SSE}  diff={abs(sse_py - PUBLISHED_SSE):.4f}")

# ----------------------------------------------------------------------
# 4. NLLS refit from the parameters.csv starting values
# ----------------------------------------------------------------------
parms = pd.read_csv(f"{M8}/Model8_parameters.csv").set_index("sparrowNames")
plist = [
    Parameter(n, parms.loc[n, "parmType"], float(parms.loc[n, "parmInit"]),
              float(parms.loc[n, "parmMin"]), float(parms.loc[n, "parmMax"]))
    for n in BETA_ORDER
]
fit = model.fit(plist)
print(f"\nrefit: SSE={fit['sse']:.4f} (published {PUBLISHED_SSE}), "
      f"n_parm={fit['n_parm']}, sites={fit['n_site']}")
print(f"{'parameter':14s} {'published':>14s} {'refit':>14s} {'SE(refit)':>12s}")
for i, n in enumerate(BETA_ORDER):
    print(f"{n:14s} {PUBLISHED[n]:14.6f} {fit['beta'][i]:14.6f} "
          f"{fit['se'][i]:12.6f}")
