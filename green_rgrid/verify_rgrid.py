"""Verification for the GREEN-Rgrid adaptation (green_rgrid.py).

Gates:
1. Grid topology builder: D8-from-DEM and D8-from-codes cases whose
   built To_catch maps are asserted exactly against the intended
   networks; the builder's NrmLengthKm values are asserted too.
2. Degeneracy: a grid whose cells coincide one-to-one with a small
   catchment network reproduces green.launch_green's loads through
   the Rgrid routing with forms="catchment" (same topological
   accumulation, Rgrid form differences switched off). The builder
   supplies the topology; per-cell attribute columns (NrmLengthKm,
   LakeFrRet) are then set from the catchment table under the cell
   mapping, because the gate is about routing, not about the
   builder's length defaults.
3. Hand-computed cell tests of the Malagó et al. (2019) forms:
   soil factor with irrigation (TN) and precipitation-only (TP),
   runoff-based river retention (Eq. 3), lake retention from
   hydraulic load (Eq. 4), the full Eq. (1) on one cell, and the
   placement of the upstream load UL in Eq. (1) on a two-cell chain.
   Expected values are computed here with math functions directly
   from the paper's equations. Paper Table 3 TP optima (alpha=1.89,
   A=0.84, B=0.76, C=0.11) are used as test inputs where a published
   value exists; n is not printed in the paper, so the tests pass
   the explicit test value n=0.05.
4. Source-apportionment identity under the Rgrid forms: per-source
   routed loads sum to the total load (an exact algebraic property
   of Eq. 1 at fixed factors).
5. Catchment path unchanged: green.py's own hand case from
   verify_port.py still reproduces its expected value.
"""

import hashlib
import json
import math
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import green  # noqa: E402
import green_rgrid as gr  # noqa: E402

RESULTS = {}
TOL = 1e-12

# ----------------------------------------------------------------------
# Gate 1: topology builder
# ----------------------------------------------------------------------

# 1a. DEM, 2x2 with one inactive cell: (0,0) and (0,1) drain to
# (1,0); (1,0) is a pit/outlet.
dem = np.array([[30.0, 30.0], [20.0, np.nan]])
down, active = gr.d8_from_dem(dem)
catch_a, cells_a = gr.build_grid_catchments(down, active,
                                            cellsize_m=1000.0)
to_a = dict(zip(catch_a["HydroID"], catch_a["To_catch"]))
assert to_a == {0: 2, 1: 2, 2: -1}, to_a
# Flow lengths: 0->2 cardinal (1 km), 1->2 diagonal (sqrt2 km),
# 2 outlet (1 km); normalised by sqrt2.
len_a = dict(zip(catch_a["HydroID"], catch_a["NrmLengthKm"]))
exp_len = {0: 1 / math.sqrt(2), 1: 1.0, 2: 1 / math.sqrt(2)}
assert all(abs(len_a[k] - v) < 1e-15 for k, v in exp_len.items()), len_a
shreve_a = dict(zip(catch_a["HydroID"], catch_a["Shreve"]))
assert shreve_a == {0: 1, 1: 1, 2: 2}, shreve_a
RESULTS["gate1a_dem_topology"] = {"to_catch": {str(k): int(v) for k, v
                                               in to_a.items()},
                                  "nrm_length": len_a, "passed": True}

# 1b. ESRI codes, 3x3, everything drains to cell 7, the outlet.
codes = np.array([[2, 4, 8], [1, 4, 16], [1, 0, 16]])
down_b, active_b = gr.d8_from_codes(codes)
catch_b, _ = gr.build_grid_catchments(down_b, active_b)
to_b = dict(zip(catch_b["HydroID"], catch_b["To_catch"]))
exp_b = {0: 4, 1: 4, 2: 4, 3: 4, 4: 7, 5: 4, 6: 7, 7: -1, 8: 7}
assert to_b == exp_b, to_b
RESULTS["gate1b_codes_topology"] = {"to_catch": {str(k): int(v) for k, v
                                                 in to_b.items()},
                                    "passed": True}
print("gate 1 (topology builder) passed")

# ----------------------------------------------------------------------
# Gate 2: degeneracy against green.launch_green
# ----------------------------------------------------------------------

SD = 0.3  # sd_coef shared by both engines in case (a)


def grid_from_network(network, nrow, ncol, active_cells):
    """Build grid catchments from an explicit downstream map."""
    down = np.full((nrow, ncol), -1, dtype=int)
    act = np.zeros((nrow, ncol), dtype=bool)
    for (r, c), d in network.items():
        act[r, c] = True
        down[r, c] = d
    return gr.build_grid_catchments(down, act)


# (a) The branching TP case of verify_port.py: catchments 1, 2 -> 3,
# catchment 3 has lake retention 0.1.
catch_tp = pd.DataFrame({
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
annual_tp = pd.DataFrame(rows)
alpha_p, alpha_l = 20.0, 0.05
ref = green.launch_green(catch_tp, annual_tp, alpha_p, alpha_l, SD,
                         [2001])
# Grid: cells (0,0)=catchment 1, (0,1)=catchment 2, (1,0)=catchment 3.
gcatch, gcells = grid_from_network({(0, 0): 2, (0, 1): 2, (1, 0): -1},
                                   2, 2, None)
idmap = {0: 1, 1: 2, 2: 3}  # grid HydroID -> catchment HydroID
assert dict(zip(gcatch["HydroID"], gcatch["To_catch"])) == \
    {0: 2, 1: 2, 2: -1}
# Attribute columns from the catchment table under the mapping.
attr = catch_tp.set_index("HydroID")
gcatch["NrmLengthKm"] = [attr.loc[idmap[h], "NrmLengthKm"]
                         for h in gcatch["HydroID"]]
gcatch["LakeFrRet"] = [attr.loc[idmap[h], "LakeFrRet"]
                       for h in gcatch["HydroID"]]
grows = []
for _, r in annual_tp.iterrows():
    gid = {v: k for k, v in idmap.items()}[r["HydroID"]]
    sum_ds = r["Min"] + r["Man"] + r["Bg"] * (1 - 0.5)  # Fix=Soil=0
    sum_ps = r["Ps"] + 0.38 * r["Bg"] * 0.5 + SD * r["Sd"]
    grows.append({"HydroID": gid, "YearValue": 2001,
                  "Surplus": sum_ds, "Ps": sum_ps,
                  "InvNrmRain": r["InvNrmRain"]})
got = gr.launch_green_rgrid(gcatch, pd.DataFrame(grows), [2001],
                            forms="catchment", alpha_p=alpha_p,
                            alpha_l=alpha_l)
diff_a = max(abs(float(ref.loc[ref["HydroID"] == idmap[g],
                               "CatchLoad"].iloc[0])
                 - float(got.loc[got["HydroID"] == g,
                                 "CellLoad"].iloc[0]))
             for g in idmap)
RESULTS["gate2a_degeneracy_tp"] = {"max_abs_diff": diff_a,
                                   "passed": bool(diff_a < TOL)}
assert diff_a < TOL, diff_a

# (b) TN, four catchments (1, 2 -> 3 -> 4), two years.
catch_tn = pd.DataFrame({
    "HydroID": [1, 2, 3, 4], "To_catch": [3, 3, 4, -1],
    "LakeFrRet": [0.0, 0.05, 0.0, 0.1],
    "NrmLengthKm": [0.3, 0.6, 0.8, 1.0]})
src = {  # HydroID: (Atm, Min, Man, Fix, Soil, Sd, Ps, FF, InvNrmRain)
    2005: {1: (80., 40., 25., 8., 4., 15., 5., 0.3, 0.7),
           2: (60., 10., 30., 0., 6., 0., 0., 0.5, 0.9),
           3: (90., 55., 12., 5., 3., 20., 9., 0.2, 0.6),
           4: (70., 20., 18., 2., 5., 8., 4., 0.4, 0.8)},
    2006: {1: (85., 42., 22., 8., 4., 14., 6., 0.3, 0.65),
           2: (58., 12., 28., 0., 6., 1., 0., 0.5, 0.85),
           3: (95., 50., 15., 5., 3., 22., 8., 0.2, 0.55),
           4: (66., 22., 16., 2., 5., 9., 5., 0.4, 0.75)}}
rows = []
for yr, per in src.items():
    for hid, (atm, mn, man, fix, soil, sd_, ps, ff, inv) in per.items():
        rows.append({"BasinID": 1, "YearValue": yr, "HydroID": hid,
                     "NextDownID": {1: 3, 2: 3, 3: 4, 4: -1}[hid],
                     "Atm": atm, "Min": mn, "Man": man, "Fix": fix,
                     "Soil": soil, "Sd": sd_, "Ps": ps,
                     "YearlyMass": np.nan, "ForestFraction": ff,
                     "InvNrmRain": inv})
annual_tn = pd.DataFrame(rows)
ap2, al2, sd2 = 15.0, 0.04, 0.25
ref_tn = green.launch_green(catch_tn, annual_tn, ap2, al2, sd2,
                            [2005, 2006])
gcatch2, _ = grid_from_network({(0, 0): 2, (0, 1): 2, (1, 0): 3,
                                (1, 1): -1}, 2, 2, None)
idmap2 = {0: 1, 1: 2, 2: 3, 3: 4}
attr2 = catch_tn.set_index("HydroID")
gcatch2["NrmLengthKm"] = [attr2.loc[idmap2[h], "NrmLengthKm"]
                          for h in gcatch2["HydroID"]]
gcatch2["LakeFrRet"] = [attr2.loc[idmap2[h], "LakeFrRet"]
                        for h in gcatch2["HydroID"]]
grows = []
for _, r in annual_tn.iterrows():
    gid = {v: k for k, v in idmap2.items()}[r["HydroID"]]
    sum_ds = (r["Min"] + r["Man"] + r["Fix"] + r["Soil"]
              + r["Atm"] * (1 - r["ForestFraction"]))
    sum_ps = (r["Ps"] + 0.38 * r["Atm"] * r["ForestFraction"]
              + sd2 * r["Sd"])
    grows.append({"HydroID": gid, "YearValue": r["YearValue"],
                  "Surplus": sum_ds, "Ps": sum_ps,
                  "InvNrmRain": r["InvNrmRain"]})
got_tn = gr.launch_green_rgrid(gcatch2, pd.DataFrame(grows),
                               [2005, 2006], forms="catchment",
                               alpha_p=ap2, alpha_l=al2)
diff_b = 0.0
for _, r in ref_tn.iterrows():
    gid = {v: k for k, v in idmap2.items()}[r["HydroID"]]
    g = float(got_tn.loc[(got_tn["HydroID"] == gid)
                         & (got_tn["Year"] == r["Year"]),
                         "CellLoad"].iloc[0])
    diff_b = max(diff_b, abs(g - float(r["CatchLoad"])))
RESULTS["gate2b_degeneracy_tn"] = {"max_abs_diff": diff_b,
                                   "passed": bool(diff_b < TOL)}
assert diff_b < TOL, diff_b
print("gate 2 (degeneracy) passed:", diff_a, diff_b)

# ----------------------------------------------------------------------
# Gate 3: hand-computed Malagó form tests
# ----------------------------------------------------------------------

# Soil factor with irrigation (TN): P=[500,250], I=[300,0], alpha=2.
# water=[800,250]; inv=[1/800,1/250]; x=[0.3125,1]; S=[e^-0.625, e^-2].
s_tn = gr.soil_reduction_factor([500., 250.], [300., 0.], 2.0,
                                nutrient="TN")
exp_tn = [math.exp(-0.625), math.exp(-2.0)]
d = max(abs(a - b) for a, b in zip(s_tn, exp_tn))
assert d < TOL, d
# Soil factor for TP: precipitation only. water=[500,250]; x=[0.5,1].
s_tp = gr.soil_reduction_factor([500., 250.], [300., 0.], 2.0,
                                nutrient="TP")
exp_tp = [math.exp(-1.0), math.exp(-2.0)]
d2 = max(abs(a - b) for a, b in zip(s_tp, exp_tp))
assert d2 < TOL, d2
RESULTS["gate3_soil_factor"] = {"tn_max_abs_diff": float(d),
                                "tp_max_abs_diff": float(d2),
                                "passed": True}

# River retention, Eq. (3), paper TP optima A=0.84, B=0.76, qs=10.
r = float(gr.river_reduction_factor([10.0], 0.84, 0.76)[0])
exp_r = 1.0 / (1.0 + 0.84 * 10.0 ** 0.76)
assert abs(r - exp_r) < TOL
RESULTS["gate3_river_factor"] = {"expected": exp_r, "got": r,
                                 "abs_diff": abs(r - exp_r),
                                 "passed": True}

# Lake retention, Eq. (4), paper TP C=0.11, test n=0.05, HL=10;
# NaN (no lake) must give 0.
res = gr.lake_retention_factor([10.0, np.nan], 0.11, 0.05)
exp_res = 1.0 / (1.0 + 0.05 * 10.0 ** 0.11)
assert abs(float(res[0]) - exp_res) < TOL and float(res[1]) == 0.0
RESULTS["gate3_lake_factor"] = {"expected": exp_res,
                                "got": float(res[0]),
                                "abs_diff": abs(float(res[0]) - exp_res),
                                "no_lake_value": float(res[1]),
                                "passed": True}

# Full Eq. (1), one cell, TN: SUR=120, PS=30, P=400, I=100 (single
# cell -> x=1), alpha=1.89 (paper TP alpha, used here as a value),
# qs=5 with A=0.84, B=0.76; HL=8, C=0.11, n=0.05.
one_catch, one_cells = gr.build_grid_catchments(
    np.array([[-1]]), np.array([[True]]))
one_annual = gr.assemble_grid_annual(
    one_cells, 2010, {"Surplus": np.array([[120.0]])},
    np.array([[30.0]]), np.array([[400.0]]), np.array([[100.0]]),
    np.array([[5.0]]), np.array([[8.0]]))
one = gr.launch_green_rgrid(one_catch, one_annual, [2010],
                            nutrient="TN", alpha=1.89, A=0.84,
                            B=0.76, C=0.11, n=0.05)
S = math.exp(-1.89)
R = 1.0 / (1.0 + 0.84 * 5.0 ** 0.76)
RES = 1.0 / (1.0 + 0.05 * 8.0 ** 0.11)
exp_one = (120.0 * S + 30.0) * R * (1.0 - RES)
got_one = float(one["CellLoad"].iloc[0])
assert abs(got_one - exp_one) < TOL, (got_one, exp_one)
assert one["Mode"].iloc[0] == gr.MODE_MALAGO
RESULTS["gate3_full_eq1"] = {"expected": exp_one, "got": got_one,
                             "abs_diff": abs(got_one - exp_one),
                             "passed": True}

# Upstream-load placement, two-cell chain: alpha=1, equal precip ->
# x=1 in both cells, S=e^-1; qs=0 -> R=1; no lakes.
# L_up = 100 e^-1; L_down = (50 e^-1 + 10 + 100 e^-1) = 150 e^-1 + 10.
chain_catch, chain_cells = gr.build_grid_catchments(
    np.array([[1, -1]]), np.array([[True, True]]))
chain_annual = pd.DataFrame([
    {"HydroID": 0, "YearValue": 2010, "Surplus": 100.0, "Ps": 0.0,
     "Precip": 600.0, "Irrig": 0.0, "SpecRunoff": 0.0,
     "HydLoad": np.nan},
    {"HydroID": 1, "YearValue": 2010, "Surplus": 50.0, "Ps": 10.0,
     "Precip": 600.0, "Irrig": 0.0, "SpecRunoff": 0.0,
     "HydLoad": np.nan}])
chain = gr.launch_green_rgrid(chain_catch, chain_annual, [2010],
                              nutrient="TN", alpha=1.0, A=0.84,
                              B=0.76, C=0.11, n=0.05)
loads = dict(zip(chain["HydroID"], chain["CellLoad"]))
exp_up = 100.0 * math.exp(-1.0)
exp_down = 150.0 * math.exp(-1.0) + 10.0
assert abs(loads[0] - exp_up) < TOL and abs(loads[1] - exp_down) < TOL
wrong_down = 50.0 * math.exp(-1.0) + 10.0 + 100.0 * math.exp(-2.0)
assert abs(loads[1] - wrong_down) > 1.0  # test discriminates
RESULTS["gate3_upstream_placement"] = {
    "expected_down": exp_down, "got_down": float(loads[1]),
    "abs_diff": abs(float(loads[1]) - exp_down), "passed": True}
print("gate 3 (hand-computed Malagó forms) passed")

# ----------------------------------------------------------------------
# Gate 4: source-apportionment identity under the Rgrid forms
# ----------------------------------------------------------------------

# 2x3 grid: 0->3, 1->4, 2->5, 3->4, 4->5, 5->outlet; cell 5 has a lake.
grid_catch, grid_cells = gr.build_grid_catchments(
    np.array([[3, 4, 5], [4, 5, -1]]),
    np.array([[True] * 3, [True] * 3]))
rng = np.random.default_rng(20261121)
frames = []
for yr in (2010, 2011):
    frames.append(gr.assemble_grid_annual(
        grid_cells, yr,
        {"SurplusAgri": rng.uniform(0, 200, (2, 3)),
         "SurplusBack": rng.uniform(0, 40, (2, 3))},
        rng.uniform(0, 30, (2, 3)), rng.uniform(300, 900, (2, 3)),
        rng.uniform(0, 150, (2, 3)), rng.uniform(1, 25, (2, 3)),
        np.array([[np.nan, np.nan, np.nan],
                  [np.nan, np.nan, 12.0]])))
grid_annual = pd.concat(frames, ignore_index=True)
pars = dict(nutrient="TP", alpha=1.89, A=0.84, B=0.76, C=0.11, n=0.05,
            surplus_components=("SurplusAgri", "SurplusBack"))
total = gr.launch_green_rgrid(grid_catch, grid_annual, [2010, 2011],
                              **pars)
shares = gr.green_rgrid_shares(grid_catch, grid_annual, [2010, 2011],
                               **pars)
src_cols = ["SurplusAgri", "SurplusBack", "Ps"]
ident = (shares[src_cols].sum(axis=1) - shares["CellLoad"]).abs()
mg = shares.merge(total, on=["HydroID", "Year"],
                  suffixes=("_sh", "_tot"))
vs_launch = (mg["CellLoad_sh"] - mg["CellLoad_tot"]).abs()
RESULTS["gate4_shares_identity"] = {
    "max_abs_diff_shares_sum": float(ident.max()),
    "max_abs_diff_vs_launch": float(vs_launch.max()),
    "n_rows": int(len(shares)),
    "passed": bool(ident.max() < TOL and vs_launch.max() < TOL)}
assert ident.max() < TOL and vs_launch.max() < TOL
print("gate 4 (shares identity) passed:", float(ident.max()),
      float(vs_launch.max()))

# ----------------------------------------------------------------------
# Gate 5: catchment path unchanged
# ----------------------------------------------------------------------

catch1 = pd.DataFrame({
    "HydroID": [1], "To_catch": [-1], "Shreve": [1],
    "LakeFrRet": [0.0], "NrmLengthKm": [0.5]})
annual1 = pd.DataFrame({
    "BasinID": [1], "YearValue": [2000], "HydroID": [1],
    "NextDownID": [-1], "Atm": [100.0], "Min": [50.0], "Man": [30.0],
    "Fix": [10.0], "Soil": [5.0], "Sd": [20.0], "Ps": [7.0],
    "YearlyMass": [np.nan], "ForestFraction": [0.25],
    "InvNrmRain": [0.8]})
res1 = green.launch_green(catch1, annual1, 10.0, 0.02, 0.2, [2000])
aa = math.exp(-10.0 * 0.8)
ctr = (50 + 30 + 10 + 5 + 100 * 0.75) * aa + 7 + 0.38 * 100 * 0.25 \
    + 0.2 * 20
exp1 = ctr * math.exp(-0.02 * 0.5)
got1 = float(res1["CatchLoad"].iloc[0])
assert abs(got1 - exp1) < TOL
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "green.py"), "rb") as f:
    green_md5 = hashlib.md5(f.read()).hexdigest()
RESULTS["gate5_catchment_path"] = {"expected": exp1, "got": got1,
                                   "abs_diff": abs(got1 - exp1),
                                   "green_py_md5": green_md5,
                                   "passed": True}
print("gate 5 (catchment path unchanged) passed")

out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "adaptation_verification_results.json")
with open(out_path, "w") as f:
    json.dump(RESULTS, f, indent=1)
print(json.dumps(RESULTS, indent=1))
print("RGRID VERIFICATION COMPLETED")
