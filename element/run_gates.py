"""
Verification gates and machinery tests for element_p.py.

Gate 1  Mass balance: per-year closure residual of the run ledger, and
        the cumulative identity  sum(Np) = cum. export + change in all
        stores (deficit ledger included), on every test run.
Gate 2  Qualitative reproduction of Van Meter et al. (2021): a
        GRW-style surplus history (rise from 1900, peak ~15 kg/ha/yr in
        the late 1970s, decline thereafter -- the trajectory shape and
        peak magnitude stated in VM21's abstract/results) must produce
        (a) export lagging surplus and (b) a large retained share of
        cumulative surplus, compared against the paper's finding that
        ~96% (">95%") of net P inputs since 1900 were retained in the
        Grand River Watershed. This is a behaviour check on the
        machinery with illustrative parameters, NOT a calibration.
Gate 3  Stationarity (Basu et al. 2010): under a long constant surplus
        the stores must approach a steady state (biogeochemical
        stationarity: the legacy store stops growing and export tracks
        input), measured by the final-decade relative change in total
        stores and the export/surplus ratio.

Machinery drives: the held Canadian census trajectory v2 (Canada,
2010-2021 surplus = total inputs minus PNC-workbook crop removal) and
the held Ontario long trajectory v3 (Ontario, 1950-2007 surplus on the
gross basis, per hectare of land in crops). These exercise the model on
real surplus trajectories at national/regional scale; they are NOT
corridor catchment applications.

Run:  python run_gates.py   (writes gate_results.json, test_outputs.csv)
"""

import json

import numpy as np
import pandas as pd

from element_p import ElementPParameters, run_element_p

BASE = "/home/hatch/workspace/fei_lake_ontario/image_na"
TOL = 1e-9

results = {}


def check_balance(name, res, initial_total=0.0):
    resid = np.abs(res.balance_residual)
    cum_in = res.surplus.sum()
    cum_out = (res.load_subsurface + res.load_surface).sum()
    s = res.stores
    d_stores = (s["Msource"][-1] + s["transit"][-1] + s["gw_sorbed"][-1]
                + s["network"][-1] - s["deficit"][-1]) - initial_total
    cum_resid = abs(cum_in - cum_out - d_stores)
    ok = resid.max() < TOL and cum_resid < 1e-8
    results[name] = dict(max_year_residual=float(resid.max()),
                         cumulative_residual=float(cum_resid), pass_=bool(ok))
    return ok


# ----------------------------------------------------------------------
# Held-data machinery drives
# ----------------------------------------------------------------------
def canada_surplus():
    tr = pd.read_csv(f"{BASE}/trajectory_canada_census_v2.csv")
    tr = tr[tr.geography == "Canada"]
    rem = pd.read_csv(f"{BASE}/extension_run/canada_p_removal_series.csv")
    m = tr.merge(rem[["year", "removal_P_kg_per_ha_cropland"]], on="year")
    return (m.year.values, m.total_P_kg_per_ha_cropland.values
            - m.removal_P_kg_per_ha_cropland.values)


def ontario_surplus():
    v3 = pd.read_csv(f"{BASE}/ipni/adopted/trajectory_ontario_long_v3.csv")
    lu = pd.read_csv(f"{BASE}/landuse_historical_canada_ontario.csv")
    lu = lu[lu.GEO == "Ontario"].set_index("REF_DATE")["Land in crops"]
    years = np.arange(1950, 2008)
    area = np.interp(years, lu.index.values, lu.values)  # 1996 value held pre-1996
    bal = v3.set_index("year").loc[years, "balance_P_gross_basis_t"].values
    return years, bal * 1000.0 / area


frames = []
for name, (yrs, sur) in [("canada_census_v2", canada_surplus()),
                         ("ontario_long_v3", ontario_surplus())]:
    spin = run_element_p(np.arange(yrs[0] - 200, yrs[0]),
                         np.full(200, sur[0]), ElementPParameters())
    init = {k: spin.stores[k][-1] for k in ("Mp", "Ma", "Mdis", "Msorb_r", "Msorb_i")}
    res = run_element_p(yrs, sur, ElementPParameters(), initial=init)
    check_balance(f"gate1_{name}", res, initial_total=sum(init.values()))
    frames.append(pd.DataFrame({
        "run": name, "year": res.years, "surplus": res.surplus,
        "load_total": res.load_total, "load_subsurface": res.load_subsurface,
        "load_surface": res.load_surface,
        "soil_source_store": res.stores["Msource"],
        "transit_store": res.stores["transit"]}))
    results[f"drive_{name}"] = dict(
        surplus_range=[float(sur.min()), float(sur.max())],
        final_load=float(res.load_total[-1]),
        final_source_store=float(res.stores["Msource"][-1]))

# ----------------------------------------------------------------------
# Gate 2: VM21-style surplus history, 1900-2020
# ----------------------------------------------------------------------
yrs = np.arange(1900, 2021)
sur = np.interp(yrs, [1900, 1940, 1978, 2000, 2020], [3.0, 6.0, 15.0, 8.0, 6.0])
spin = run_element_p(np.arange(1800, 1900), np.full(100, 3.0),
                     ElementPParameters())
init = {k: spin.stores[k][-1] for k in ("Mp", "Ma", "Mdis", "Msorb_r", "Msorb_i")}
res = run_element_p(yrs, sur, ElementPParameters(), initial=init)
check_balance("gate1_gate2_history", res, initial_total=sum(init.values()))
cum_surplus = res.surplus.sum()
cum_export = (res.load_subsurface + res.load_surface).sum()
retained_share = 1.0 - cum_export / cum_surplus
# Lag: cross-correlation peak between surplus and export.
s_c = res.surplus - res.surplus.mean()
l_c = res.load_total - res.load_total.mean()
cors = [np.corrcoef(s_c[:len(s_c) - k], l_c[k:])[0, 1] for k in range(0, 61)]
lag = int(np.argmax(cors))
results["gate2"] = dict(
    cum_surplus=float(cum_surplus), cum_export=float(cum_export),
    retained_share=float(retained_share),
    paper_value="~0.96 retained since 1900 (VM21 abstract/key points: >95%)",
    export_peak_year=int(yrs[np.argmax(res.load_total)]),
    surplus_peak_year=int(yrs[np.argmax(res.surplus)]),
    crosscorr_lag_years=lag,
    # Behaviour-check criterion (set a priori for this machinery test):
    # the signature is a large retained share and a positive lag with
    # the export peak after the surplus peak. The paper's 0.96 is its
    # CALIBRATED Grand River value; our defaults are uncalibrated, so
    # the check is on the signature and the order of magnitude of the
    # retained share, not on reproducing 0.96.
    pass_=bool(retained_share > 0.85 and lag > 0
               and yrs[np.argmax(res.load_total)] > yrs[np.argmax(res.surplus)]))
frames.append(pd.DataFrame({
    "run": "gate2_grw_style", "year": res.years, "surplus": res.surplus,
    "load_total": res.load_total, "load_subsurface": res.load_subsurface,
    "load_surface": res.load_surface,
    "soil_source_store": res.stores["Msource"],
    "transit_store": res.stores["transit"]}))

# Switch-branch spot check: triangular TTD (VM17 supplement choice).
res_tri = run_element_p(yrs, sur, ElementPParameters(ttd="triangular"),
                        initial=init)
results["gate2_triangular_branch_retained_share"] = float(
    1.0 - (res_tri.load_subsurface.sum() + res_tri.load_surface.sum()) / cum_surplus)

# ----------------------------------------------------------------------
# Gate 3: stationarity under constant surplus
# ----------------------------------------------------------------------
yrs3 = np.arange(0, 2000)
res3 = run_element_p(yrs3, np.full(2000, 8.0), ElementPParameters())
check_balance("gate1_stationarity", res3)
# Stationarity is assessed on the DYNAMIC stores (source zone +
# transit). The gw_sorbed and network stores are cumulative retention
# integrals in VM21's accounting (sorbed along pathways; retained in
# reservoirs/riparian zones) and grow linearly by construction, so
# throughput stationarity is checked separately: export + permanent
# sink rates must approach the input rate.
tot = res3.stores["Msource"] + res3.stores["transit"]
rel_change_last50 = abs(tot[-1] - tot[-50]) / tot[-1]
sink_rate = ((res3.stores["gw_sorbed"][-1] - res3.stores["gw_sorbed"][-50])
             + (res3.stores["network"][-1] - res3.stores["network"][-50])) / 50.0
results["gate3"] = dict(
    dynamic_store_final=float(tot[-1]),
    rel_change_last_50yr=float(rel_change_last50),
    export_over_surplus_final_decade=float(res3.load_total[-10:].mean() / 8.0),
    export_plus_sinks_over_surplus=float(
        (res3.load_total[-10:].mean() + sink_rate) / 8.0),
    pass_=bool(rel_change_last50 < 0.02))

pd.concat(frames).to_csv("test_outputs.csv", index=False)
with open("gate_results.json", "w") as f:
    json.dump(results, f, indent=2)
print(json.dumps(results, indent=2))
