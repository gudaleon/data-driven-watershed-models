"""
verify_against_guide.py -- verification gate for the Python EEMMA port.

Reference values are the R script's own printed outputs on the demo
data, as published in Kirchner's *eemma user's guide* (EnviDat
doi:10.16904/envidat.410), which reproduces the vignette analyses of
eemma_vignette_script.r.  (R itself is not installed on this machine --
see PORT_NOTES.md -- so the guide's printed tables are the reference.)

Two independent brute-force cross-checks are included for the OLS
machinery:
  * complete case: constrained least squares sum(f) = 1 solved directly
    via the KKT system (no difference transformation);
  * incomplete case: plain unconstrained least squares with intercept,
    with f_unk = 1 - sum(f).

The EIV path has no published reference numbers in the guide; it is
exercised as a smoke test and its output printed for the record.

Usage:  python verify_against_guide.py [path/to/vignette_data.txt]
Exit code 0 if every published value matches within 5e-6 (the guide
prints 6-9 significant digits), else 1.
"""

import sys

import numpy as np
import pandas as pd

from eemma import eemma

TOL = 5e-6

# (title, y col, X cols, kwargs, guide f, guide se)
GUIDE = [
    ("xylem1 complete", "xylem1", ["precip", "soil", "gw"], {},
     [0.01458161, 0.56922839, 0.41619000],
     [0.01107822, 0.04771557, 0.04883837]),
    ("xylem1 incomplete, gw missing", "xylem1", ["precip", "soil"],
     {"complete": False},
     [0.01628097, 0.52489645, -4.57869373, 0.45882257],
     [0.01085917, 0.06343070, 0.62528608, 0.06443678]),
    ("xylem1 incomplete, all EMs", "xylem1", ["precip", "soil", "gw"],
     {"complete": False},
     [0.01502141, 0.53569228, 0.18118903, -2.67960750, 0.26809728],
     [0.01067974, 0.06267287, 0.13817945, 1.57246935, 0.15855543]),
    ("xylem1 complete, gw omitted", "xylem1", ["precip", "soil"], {},
     [0.03474711, 0.96525289],
     [0.02363558, 0.02363558]),
    ("xylem1 memory", "xylem1", ["precip", "soil", "gw"],
     {"memory": True},
     [0.02102953, 0.57081759, 0.40815288, 0.13393360],
     [0.01485511, 0.05857590, 0.06163491, 0.13135832]),
    ("xylem2 memory", "xylem2", ["precip", "soil", "gw"],
     {"memory": True},
     [-0.002256761, 0.564638607, 0.437618153, 0.770389817],
     [0.05352275, 0.26108935, 0.26765053, 0.11561457]),
    ("xylem2 no memory", "xylem2", ["precip", "soil", "gw"], {},
     [-0.00887984, 0.40622847, 0.60265137],
     [0.02256012, 0.09716985, 0.09945635]),
    ("stream complete", "stream", ["precip", "soil", "gw"], {},
     [0.2068911, 0.2101869, 0.5829221],
     [0.01141577, 0.04916945, 0.05032645]),
    ("stream incomplete, all EMs", "stream", ["precip", "soil", "gw"],
     {"complete": False},
     [0.2070744, 0.1713226, 0.4222108, -1.9784636, 0.1993922],
     [0.01163847, 0.06829907, 0.15058395, 1.71363139, 0.17278910]),
    ("stream memory", "stream", ["precip", "soil", "gw"],
     {"memory": True},
     [0.2229728, 0.1774278, 0.5995994, 0.0552703],
     [0.02364249, 0.05755744, 0.05431148, 0.06634318]),
    ("stream incomplete, precip only", "stream", ["precip"],
     {"complete": False},
     [0.2098147, -7.8219550, 0.7901853],
     [0.01459053, 0.14846348, 0.01459053]),
]

# guide, stream complete case, $f_cov (printed in the user's guide)
GUIDE_FCOV_STREAM = np.array([
    [1.303198e-04, -7.601319e-06, -0.0001227185],
    [-7.601319e-06, 2.417635e-03, -0.0024100334],
    [-1.227185e-04, -2.410033e-03, 0.0025327519],
])


def main(path="vignette_data.txt"):
    dat = pd.read_csv(path, sep="\t")
    worst = 0.0
    ok = True
    for title, ycol, xcols, kw, gf, gse in GUIDE:
        res = eemma(dat[ycol], dat[xcols], **kw)
        df_f = np.max(np.abs(res["f"].to_numpy() - np.array(gf)))
        df_se = np.max(np.abs(res["se"].to_numpy() - np.array(gse)))
        worst = max(worst, df_f, df_se)
        status = "PASS" if max(df_f, df_se) <= TOL else "FAIL"
        ok &= status == "PASS"
        print(f"[{status}] {title:36s} max|df|={df_f:.3e}  "
              f"max|dse|={df_se:.3e}")

    res = eemma(dat["stream"], dat[["precip", "soil", "gw"]])
    d = np.max(np.abs(res["f_cov"].to_numpy() - GUIDE_FCOV_STREAM))
    status = "PASS" if d <= 2e-9 else "FAIL"
    ok &= status == "PASS"
    print(f"[{status}] {'stream complete f_cov':36s} max|dCov|={d:.3e}")

    # ---- independent cross-check 1: KKT constrained least squares ----
    y = dat["stream"].to_numpy()
    X = dat[["precip", "soil", "gw"]].to_numpy()
    n, k = X.shape
    KKT = np.block([[X.T @ X, np.ones((k, 1))],
                    [np.ones((1, k)), np.zeros((1, 1))]])
    rhs = np.concatenate([X.T @ y, [1.0]])
    f_kkt = np.linalg.solve(KKT, rhs)[:k]
    d = np.max(np.abs(f_kkt - res["f"].to_numpy()))
    status = "PASS" if d <= 1e-10 else "FAIL"
    ok &= status == "PASS"
    print(f"[{status}] {'KKT cross-check (stream complete)':36s} "
          f"max|df|={d:.3e}")

    # ---- independent cross-check 2: plain lstsq, incomplete case ----
    y = dat["xylem1"].to_numpy()
    X = dat[["precip", "soil"]].to_numpy()
    A = np.hstack([X, np.ones((n, 1))])
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    f_ref = np.array([beta[0], beta[1], beta[2], 1.0 - beta[0] - beta[1]])
    res2 = eemma(dat["xylem1"], dat[["precip", "soil"]], complete=False)
    d = np.max(np.abs(f_ref - res2["f"].to_numpy()))
    status = "PASS" if d <= 1e-10 else "FAIL"
    ok &= status == "PASS"
    print(f"[{status}] {'lstsq cross-check (xylem1 incomplete)':36s} "
          f"max|df|={d:.3e}")

    # ---- EIV smoke test (no published reference values) ----
    print("\nEIV smoke test, stream complete (defaults: equal relative "
          "errors, alpha=1):")
    res_eiv = eemma(dat["stream"], dat[["precip", "soil", "gw"]],
                    method="EIV")
    print(res_eiv["coeffs"].to_string(float_format=lambda v: f"{v:.8g}"))
    print(f"lambda = {res_eiv['lambda']:.8g}, "
          f"RMSE = {res_eiv['RMSE']:.8g}")
    print("\nEIV smoke test, xylem1 incomplete (precip+soil):")
    res_eiv2 = eemma(dat["xylem1"], dat[["precip", "soil"]],
                     complete=False, method="EIV")
    print(res_eiv2["coeffs"].to_string(float_format=lambda v: f"{v:.8g}"))
    print(f"lambda = {res_eiv2['lambda']:.8g}, "
          f"RMSE = {res_eiv2['RMSE']:.8g}")

    print(f"\nWorst published-value discrepancy: {worst:.3e} "
          f"(tolerance {TOL:.0e})")
    print("OVERALL:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "vignette_data.txt"))
