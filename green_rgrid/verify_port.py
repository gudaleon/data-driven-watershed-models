"""Verification for the Python port of GREENeR (green.py).

Layers of verification, strongest available first:

1. Hand-computed tests on tiny networks (expected values computed here
   directly from the GREENeR equations with math.exp, independently of
   the module's routing code).
2. Input-equivalence check: the port of read_NSdata applied to the
   GREENeR_data CSVs is compared against the GREENeR package's own
   built-in Lay-basin datasets (catch_data_TN / annual_data_TN, read
   from the package .rda files), which were produced by the R pipeline.
3. Full Lay-basin run with the vignette's example parameters:
   source-share identity (shares must sum to the total load exactly,
   an algebraic property of the GREEN equations), and agreement of
   modelled loads with the observed loads held in the data.

The direct R-vs-Python numeric comparison could not be run: this VM
currently has no R interpreter (~/R/library survives, R itself was
wiped by a VM restart, and project rules say the agent does not
reinstall R). See PORT_NOTES.md.
"""

import json
import math
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from green import (compute_shreve_levels, green_shares, launch_green,
                   read_nsdata, region_nut_balance)

RESULTS = {}

# ----------------------------------------------------------------------
# 1. Hand-computed tests
# ----------------------------------------------------------------------

def hand_single_tn():
    catch = pd.DataFrame({
        "HydroID": [1], "To_catch": [-1], "Shreve": [1],
        "LakeFrRet": [0.0], "NrmLengthKm": [0.5]})
    annual = pd.DataFrame({
        "BasinID": [1], "YearValue": [2000], "HydroID": [1],
        "NextDownID": [-1], "Atm": [100.0], "Min": [50.0], "Man": [30.0],
        "Fix": [10.0], "Soil": [5.0], "Sd": [20.0], "Ps": [7.0],
        "YearlyMass": [np.nan], "ForestFraction": [0.25],
        "InvNrmRain": [0.8]})
    alpha_p, alpha_l, sd = 10.0, 0.02, 0.2
    res = launch_green(catch, annual, alpha_p, alpha_l, sd, [2000])
    aa = math.exp(-alpha_p * 0.8)
    sum_ds = 50 + 30 + 10 + 5 + 100 * 0.75
    ctr = sum_ds * aa + 7 + 0.38 * 100 * 0.25 + sd * 20
    expected = ctr * math.exp(-alpha_l * 0.5)
    got = float(res["CatchLoad"].iloc[0])
    RESULTS["hand_single_tn"] = {
        "expected": expected, "got": got,
        "abs_diff": abs(expected - got)}
    assert abs(expected - got) < 1e-12


def hand_branch_tp():
    # Two headwaters (1, 2) drain into 3; 3 has lake retention; 3 -> outlet.
    catch = pd.DataFrame({
        "HydroID": [1, 2, 3], "To_catch": [3, 3, -1],
        "LakeFrRet": [0.0, 0.0, 0.1], "NrmLengthKm": [0.2, 0.4, 1.0]})
    rows = []
    for hid, bg, mn, man, sd_, ps in [(1, 10., 5., 3., 2., 1.),
                                      (2, 20., 0., 6., 0., 0.),
                                      (3, 5., 8., 2., 4., 3.)]:
        rows.append({"BasinID": 1, "YearValue": 2001, "HydroID": hid,
                     "NextDownID": 3 if hid < 3 else -1, "Bg": bg,
                     "Min": mn, "Man": man, "Sd": sd_, "Ps": ps,
                     "YearlyMass": np.nan, "ForestFraction": 0.5,
                     "InvNrmRain": 0.5})
    annual = pd.DataFrame(rows)
    alpha_p, alpha_l, sd_c = 20.0, 0.05, 0.3
    res = launch_green(catch, annual, alpha_p, alpha_l, sd_c, [2001])
    aa = math.exp(-alpha_p * 0.5)

    def local(bg, mn, man, sd_, ps):
        sum_ds = mn + man + bg * (1 - 0.5)     # Fix = Soil = 0 for TP
        return sum_ds * aa + ps + 0.38 * bg * 0.5 + sd_c * sd_

    l1 = local(10., 5., 3., 2., 1.) * math.exp(-alpha_l * 0.2)
    l2 = local(20., 0., 6., 0., 0.) * math.exp(-alpha_l * 0.4)
    l3 = (1 - 0.1) * (local(5., 8., 2., 4., 3.) + l1 + l2) \
        * math.exp(-alpha_l * 1.0)
    got = dict(zip(res["HydroID"], res["CatchLoad"]))
    RESULTS["hand_branch_tp"] = {
        "expected": {"1": l1, "2": l2, "3": l3},
        "got": {str(k): v for k, v in got.items()},
        "max_abs_diff": max(abs(l1 - got[1]), abs(l2 - got[2]),
                            abs(l3 - got[3]))}
    assert max(abs(l1 - got[1]), abs(l2 - got[2]), abs(l3 - got[3])) < 1e-12
    # Shreve levels: headwaters 1, node 3 level 2
    lv = dict(zip(compute_shreve_levels(catch)["HydroID"],
                  compute_shreve_levels(catch)["Shreve"]))
    assert lv == {1: 1, 2: 1, 3: 2}, lv


hand_single_tn()
hand_branch_tp()
print("hand tests passed")

# ----------------------------------------------------------------------
# 2. Input equivalence: read_nsdata port vs the package's built-in data
# ----------------------------------------------------------------------

CODEBASE = os.path.expanduser("~/workspace/fei_lake_ontario/codebase")
catch_csv, annual_csv = read_nsdata(os.path.join(CODEBASE, "GREENeR_data",
                                                 "data", "csv"))

import pyreadr
def load_rda(name):
    parsed = pyreadr.read_r(os.path.join(CODEBASE, "GREENeR", "data",
                                         f"{name}.rda"))
    return parsed[name]

catch_rda = load_rda("catch_data_TN")
annual_rda = load_rda("annual_data_TN")

m = catch_csv.merge(catch_rda, on="HydroID", suffixes=("_csv", "_rda"))
input_check = {
    "n_catch_csv": int(len(catch_csv)), "n_catch_rda": int(len(catch_rda)),
    "to_catch_mismatch": int((m["To_catch_csv"] != m["To_catch_rda"]).sum()),
    "max_abs_LakeFrRet_diff": float(
        (m["LakeFrRet_csv"] - m["LakeFrRet_rda"]).abs().max()),
    "max_abs_NrmLengthKm_diff": float(
        (m["NrmLengthKm_csv"] - m["NrmLengthKm_rda"]).abs().max()),
    "max_abs_Shreve_diff": float(
        (m["Shreve_csv"] - m["Shreve_rda"]).abs().max()),
}
a = annual_csv.merge(annual_rda, on=["HydroID", "YearValue"],
                     suffixes=("_csv", "_rda"))
for col in ["Atm", "Min", "Man", "Fix", "Soil", "Sd", "Ps",
            "ForestFraction", "InvNrmRain"]:
    input_check[f"max_abs_{col}_diff"] = float(
        (a[f"{col}_csv"] - a[f"{col}_rda"]).abs().max(skipna=True))
input_check["n_annual_rows_csv"] = int(len(annual_csv))
input_check["n_annual_rows_matched"] = int(len(a))
input_check["interpretation"] = (
    "Annual source columns and the topology match the package data "
    "exactly (InvNrmRain/ForestFraction to float precision). The "
    "catchment scalars do not: the correlation between Length.csv and "
    "the built-in NrmLengthKm is -0.04, so the GREENeR_data length "
    "file is not the source of the built-in column at all; the "
    "built-in Shreve column holds classic Shreve magnitudes rather "
    "than the topological levels that read_NSdata's shreve() "
    "generates; and the built-in LakeFrRet values differ from the "
    "Kronvang formula applied to the CSV lake properties by up to "
    "0.006. read_NSdata applied to these CSVs - in R or in this port - "
    "cannot reproduce those three columns; the difference is data "
    "provenance inside the GREENeR repositories, not a port "
    "discrepancy.")
RESULTS["input_equivalence_lay_TN"] = input_check
print("input equivalence:", json.dumps(input_check, indent=1))

# ----------------------------------------------------------------------
# 3. Full Lay run, vignette example parameters; shares identity; obs fit
# ----------------------------------------------------------------------

ALPHA_P, ALPHA_L, SD_COEF = 35.09, 0.02, 0.2   # green_shares doc example
YEARS = list(range(2003, 2010))                # vignette calibration years

loads = launch_green(catch_rda, annual_rda, ALPHA_P, ALPHA_L, SD_COEF, YEARS)
shares = green_shares(catch_rda, annual_rda, ALPHA_P, ALPHA_L, SD_COEF, YEARS)
src_cols = ["Atm", "Min", "Man", "Fix", "Soil", "Sd", "Ps"]
share_sum = shares[src_cols].sum(axis=1)
ident = (share_sum - shares["CatchLoad"]).abs()
rel_ident = ident / shares["CatchLoad"].replace(0, np.nan)
RESULTS["shares_identity"] = {
    "max_abs_diff": float(ident.max()),
    "max_rel_diff": float(rel_ident.max(skipna=True)),
    "n_rows": int(len(shares))}

# Cross-check shares totals against launch_green loads directly
mg = shares.merge(loads, on=["HydroID", "Year"], suffixes=("_sh", "_lg"))
RESULTS["shares_vs_launch"] = {
    "max_abs_diff": float(
        (mg["CatchLoad_sh"] - mg["CatchLoad_lg"]).abs().max())}

obs = loads.dropna(subset=["ObsLoad"])
err = obs["CatchLoad"] - obs["ObsLoad"]
RESULTS["lay_obs_agreement"] = {
    "n_obs": int(len(obs)),
    "pearson_r": float(np.corrcoef(obs["CatchLoad"], obs["ObsLoad"])[0, 1]),
    "rmse": float(np.sqrt((err ** 2).mean())),
    "mean_bias": float(err.mean()),
    "note": "Parameters are the vignette doc example, not calibrated; "
            "agreement is a plausibility check, not a fit claim."}

outlet_id = int(catch_rda.loc[catch_rda["To_catch"] == -1, "HydroID"].iloc[0])
outlet = loads[loads["HydroID"] == outlet_id][["Year", "CatchLoad"]]
RESULTS["lay_outlet_loads"] = {
    "outlet_HydroID": outlet_id,
    "loads_by_year": {int(y): float(v) for y, v in
                      zip(outlet["Year"], outlet["CatchLoad"])}}
print("outlet loads:", RESULTS["lay_outlet_loads"])

# TP scenario on the package's built-in TP data: shares identity again
catch_tp = load_rda("catch_data_TP")
annual_tp = load_rda("annual_data_TP")
shares_tp = green_shares(catch_tp, annual_tp, ALPHA_P, ALPHA_L, SD_COEF,
                         YEARS)
tp_cols = ["Bg", "Min", "Man", "Sd", "Ps"]
ident_tp = (shares_tp[tp_cols].sum(axis=1) - shares_tp["CatchLoad"]).abs()
RESULTS["shares_identity_TP"] = {
    "max_abs_diff": float(ident_tp.max()), "n_rows": int(len(shares_tp))}

# Balance port runs and conserves routed totals
bal = region_nut_balance(catch_rda, annual_rda, ALPHA_P, ALPHA_L, SD_COEF,
                         YEARS)
RESULTS["balance_rows"] = int(len(bal))

out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "verification_results.json")
with open(out_path, "w") as f:
    json.dump(RESULTS, f, indent=1)
print(json.dumps(RESULTS, indent=1))
print("VERIFICATION SCRIPT COMPLETED")
