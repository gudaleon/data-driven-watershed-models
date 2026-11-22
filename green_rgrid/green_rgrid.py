"""GREEN-Rgrid adaptation of the GREEN member (Malagó et al. 2019).

This module adds the grid-cell variant of GREEN published by Malagó,
Bouraoui, Grizzetti & De Roo (2019, J. Hydrology: Regional Studies 22,
100592) to the verified GREENeR catchment port in ``green.py``. The
catchment path in ``green.py`` is NOT modified by this module and
remains reachable unchanged (``green.launch_green`` etc.); the Rgrid
machinery below reuses only its topological ordering
(``green._topological_order`` / ``green.compute_shreve_levels``), so
both modes accumulate downstream through the same ordering code.

The mode switch is explicit: ``launch_green_rgrid(..., forms=...)``
with ``forms="malago2019"`` (the published Rgrid forms) or
``forms="catchment"`` (the Rgrid routing arrangement fed with the
GREENeR catchment factors, used by the degeneracy gate in
``verify_rgrid.py``). Every output row carries a ``Mode`` column
naming the form set that produced it.

Malagó et al. (2019) model equations (paper numbering), per grid
cell i and year:

    L_i = [SUR_i * S_i * R_i + (PS_i + UL_i) * R_i] * (1 - RES_i)   (1)

where L is the annual nutrient load (ton/y), SUR the nutrient surplus
in the cell (ton/y), PS the point sources (ton/y) routed as direct
inputs, UL the upstream load (ton/y) -- added to the receiving cell
"as an additional point source" (paper section 2.2), so UL is
reduced by R_i and RES_i but NOT by the soil factor S_i -- and S, R,
RES are the soil reduction factor, river reduction factor and
lake/reservoir retention (all dimensionless):

    S_i   = exp(-alpha * x_i)                                       (2)
    R_i   = 1 / (1 + A * qs_i ** B)                                 (3)
    RES_i = 1 / (1 + n * HL_i ** C)                                 (4)

Eq. (2): x_i is the inverse of the cell's water input normalised by
the maximum inverse over the whole study area (Frank & Todeschini
1994 maximum scaling, paper section 2.2). The water input is
precipitation PLUS irrigation for nitrogen; the paper's input
section (2.3.2) states the model "uses annual precipitation and
irrigation as proxies for the basin retention of nitrogen, and only
precipitation for phosphorus", so for phosphorus the water input is
precipitation alone. alpha is a calibration parameter.

Eq. (3): qs_i is the average area specific runoff (l/s/km2), in the
paper computed from LISFLOOD annual discharge divided by the cell's
upstream area. A and B are calibration parameters. This replaces
GREENeR's reach-length river retention BB = exp(-alpha_L *
NrmLengthKm) -- the second of the two Rgrid form differences.

Eq. (4): HL_i is the hydraulic load (m), in the paper the lake's
streamflow output divided by its area; C is a calibration parameter
and n is "a fixed coefficient retrieved from Behrendt et al. (2003)".
Malagó et al. do not print the value of n, so this module takes n as
a REQUIRED named parameter and never supplies a default for it.
The same applies to alpha, A, B and C: they are calibration
parameters of the paper (its Table 3 lists Mediterranean optima per
determinand) and are required arguments here, not built-in defaults.

SUR construction (paper section 2.3.3.2), for reference: cropland
and fodder surplus = mineral + manure fertiliser + N fixation (N
only) + atmospheric deposition (N only) - crop uptake; grassland
surplus = inputs - 60% of inputs (Bouwman et al. 2005); the diffuse
pool further includes atmospheric deposition on natural areas
(natural background), erosion-borne nutrients, urban washoff,
phosphorus weathering and scattered dwellings. This module takes
SUR as given (optionally split into named components for source
apportionment) and does not rebuild the source layers.

Grid schemas:

catch_data (as in green.py): HydroID, To_catch (-1 = outlet),
LakeFrRet, NrmLengthKm. In Rgrid mode LakeFrRet is unused (lake
retention comes from HydLoad via Eq. 4); it is carried so one table
serves both modes, and ``build_grid_catchments`` fills it.

annual_data, forms="malago2019":
    HydroID, YearValue, Ps, Precip (mm/y), Irrig (mm/y),
    SpecRunoff (l/s/km2), HydLoad (m; NaN where the cell holds no
    lake), one column per surplus component (default "Surplus"),
    optional ObsLoad / YearlyMass (carried through, unused).
annual_data, forms="catchment":
    HydroID, YearValue, Surplus (the diffuse total that the soil
    factor multiplies, i.e. GREENeR's SumDS), Ps (the direct total,
    i.e. GREENeR's SumPSCmp1), InvNrmRain.

Lake rule: Eq. (4) applies only to cells holding a lake or reservoir
(the paper positions at most one lake per cell, from HydroLAKES).
A cell with HydLoad = NaN is a non-lake cell and gets RES = 0; at
HL = 0 Eq. (4) would give RES = 1, so the formula cannot be applied
to non-lake cells. HydLoad <= 0 on a lake cell follows the formula's
limit behaviour and is the caller's responsibility.

Units follow the inputs (the paper uses ton/y for loads); all
reduction factors are dimensionless.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import green as _catchment  # noqa: E402  (sibling module, same folder)

__all__ = [
    "d8_from_dem",
    "d8_from_codes",
    "build_grid_catchments",
    "assemble_grid_annual",
    "soil_water_input",
    "soil_reduction_x",
    "soil_reduction_factor",
    "river_reduction_factor",
    "lake_retention_factor",
    "launch_green_rgrid",
    "green_rgrid_shares",
    "MODE_MALAGO",
    "MODE_CATCHMENT_FORMS",
]

MODE_MALAGO = "GREEN-Rgrid (Malagó et al. 2019 forms)"
MODE_CATCHMENT_FORMS = "GREEN-Rgrid routing, GREENeR catchment forms"


# --------------------------------------------------------------------------
# Grid topology builder
# --------------------------------------------------------------------------

# Neighbour scan order for the DEM steepest-descent rule: row by row
# from north-west to south-east. On an exact tie of slope the FIRST
# neighbour in this order wins (documented determinism rule).
_NEIGHBOURS = [(-1, -1), (-1, 0), (-1, 1),
               (0, -1),           (0, 1),
               (1, -1),  (1, 0),  (1, 1)]

# ESRI D8 flow-direction codes -> (d_row, d_col).
_D8_ESRI = {1: (0, 1), 2: (1, 1), 4: (1, 0), 8: (1, -1),
            16: (0, -1), 32: (-1, -1), 64: (-1, 0), 128: (-1, 1)}


def d8_from_dem(dem, nodata=None):
    """Steepest-descent D8 downstream indices from an elevation array.

    dem: 2D array of elevations; cells equal to ``nodata`` (or NaN)
    are inactive. Returns (downstream, active): downstream is an int
    array of flat (row-major) indices of the receiving cell, -1 for
    outlets; active is a boolean mask.

    Rules (stated limits of the builder): single flow direction per
    cell (D8, no flow dispersion); the receiving cell is the active
    neighbour with the largest drop per unit distance (cardinal
    distance 1, diagonal sqrt(2), in cell units); a strictly positive
    drop is required, so pits and flats with no lower neighbour
    become outlets -- there is NO depression filling and NO flat
    resolution; ties are broken by the fixed neighbour scan order
    _NEIGHBOURS; cells draining off the grid edge are outlets.
    """
    z = np.asarray(dem, dtype=float)
    if z.ndim != 2:
        raise ValueError("dem must be a 2D array.")
    active = ~np.isnan(z)
    if nodata is not None:
        active &= z != nodata
    nrow, ncol = z.shape
    downstream = np.full(z.shape, -1, dtype=int)
    for r in range(nrow):
        for c in range(ncol):
            if not active[r, c]:
                continue
            best_slope = 0.0
            best = -1
            for dr, dc in _NEIGHBOURS:
                rr, cc = r + dr, c + dc
                if not (0 <= rr < nrow and 0 <= cc < ncol):
                    continue
                if not active[rr, cc]:
                    continue
                dist = np.hypot(dr, dc)
                slope = (z[r, c] - z[rr, cc]) / dist
                if slope > best_slope:
                    best_slope = slope
                    best = rr * ncol + cc
            downstream[r, c] = best
    return downstream, active


def d8_from_codes(codes, encoding="esri"):
    """Downstream indices from a D8 flow-direction code array.

    Only the ESRI encoding is implemented (1=E, 2=SE, 4=S, 8=SW,
    16=W, 32=NW, 64=N, 128=NE). Code 0 (ESRI "undefined") and any
    unrecognised code make the cell an outlet; a direction pointing
    off the grid or into an inactive cell also resolves to outlet.
    Every cell with a recognised code is active. Returns
    (downstream, active) as in d8_from_dem.
    """
    if encoding != "esri":
        raise ValueError(f"Unsupported D8 encoding: {encoding!r} "
                         "(only 'esri' is implemented).")
    arr = np.asarray(codes)
    if arr.ndim != 2:
        raise ValueError("codes must be a 2D array.")
    nrow, ncol = arr.shape
    active = np.isin(arr, list(_D8_ESRI) + [0])
    downstream = np.full(arr.shape, -1, dtype=int)
    for r in range(nrow):
        for c in range(ncol):
            step = _D8_ESRI.get(int(arr[r, c])) if active[r, c] else None
            if step is None:
                continue
            rr, cc = r + step[0], c + step[1]
            if 0 <= rr < nrow and 0 <= cc < ncol and active[rr, cc]:
                downstream[r, c] = rr * ncol + cc
    return downstream, active


def build_grid_catchments(downstream, active, cellsize_m=1.0,
                          lake_fr=None):
    """Express a grid as a green.py catchment table.

    downstream / active: the arrays returned by d8_from_dem or
    d8_from_codes. HydroID is the flat row-major cell index; To_catch
    is the downstream cell's HydroID, or -1 where the cell drains off
    the grid, into an inactive cell, or is itself a pit/flat outlet.
    NrmLengthKm is the cell's flow length -- the distance to its
    downstream neighbour (cellsize, or sqrt(2) * cellsize for a
    diagonal step; cellsize for outlets, whose reach still crosses
    the cell) -- normalised by the grid maximum, GREENeR's
    definition. cellsize_m is the cell size in metres. lake_fr, if
    given, is a 2D array of lake fractions of retention for the
    catchment mode (LakeFrRet); default 0.

    Returns (catch_data, grid_cells): catch_data follows the green.py
    schema (HydroID, To_catch, Shreve, LakeFrRet, NrmLengthKm) with
    Shreve levels from green.compute_shreve_levels; grid_cells maps
    HydroID -> (row, col). Raises if the topology contains a cycle.
    """
    downstream = np.asarray(downstream)
    active = np.asarray(active, dtype=bool)
    if downstream.shape != active.shape or downstream.ndim != 2:
        raise ValueError("downstream and active must be 2D arrays of "
                         "the same shape.")
    nrow, ncol = downstream.shape
    rows = []
    grid_cells = {}
    for r in range(nrow):
        for c in range(ncol):
            if not active[r, c]:
                continue
            hid = r * ncol + c
            grid_cells[hid] = (r, c)
            d = int(downstream[r, c])
            if d >= 0 and active.flat[d]:
                dr = abs(d // ncol - r)
                dc = abs(d % ncol - c)
                step_km = cellsize_m / 1000.0 * np.hypot(dr, dc)
                to = d
            else:
                step_km = cellsize_m / 1000.0
                to = -1
            rows.append({"HydroID": hid, "To_catch": to,
                         "LakeFrRet": 0.0 if lake_fr is None
                         else float(np.asarray(lake_fr)[r, c]),
                         "_len_km": step_km})
    catch = pd.DataFrame(rows)
    catch["NrmLengthKm"] = catch["_len_km"] / catch["_len_km"].max()
    catch = catch.drop(columns=["_len_km"])
    catch = _catchment.compute_shreve_levels(catch)
    return catch[["HydroID", "To_catch", "Shreve", "LakeFrRet",
                  "NrmLengthKm"]], grid_cells


def assemble_grid_annual(grid_cells, year, surplus, ps, precip, irrig,
                         spec_runoff, hyd_load):
    """Assemble one year of Rgrid annual_data from 2D field arrays.

    grid_cells: the HydroID -> (row, col) map from
    build_grid_catchments. surplus: dict {component_name: 2D array}
    of surplus masses per cell (ton/y in the paper's units); ps,
    precip, irrig, spec_runoff, hyd_load: 2D arrays in the units of
    the module docstring (hyd_load NaN where a cell holds no lake).
    This is the assembly point the documented real-grid path feeds.
    """
    rows = []
    for hid, (r, c) in sorted(grid_cells.items()):
        row = {"HydroID": hid, "YearValue": year,
               "Ps": float(ps[r, c]), "Precip": float(precip[r, c]),
               "Irrig": float(irrig[r, c]),
               "SpecRunoff": float(spec_runoff[r, c]),
               "HydLoad": float(hyd_load[r, c])}
        for name, arr in surplus.items():
            row[name] = float(arr[r, c])
        rows.append(row)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# Rgrid reduction factors (Malagó et al. 2019, Eqs. 2-4)
# --------------------------------------------------------------------------

def soil_water_input(precip, irrig, nutrient="TN", use_irrigation=None):
    """The water input driving Eq. (2)'s x, as published.

    Precipitation plus irrigation for nitrogen; precipitation only
    for phosphorus (paper section 2.3.2). ``use_irrigation`` may
    override the published per-nutrient choice explicitly.
    """
    precip = np.asarray(precip, dtype=float)
    irrig = np.asarray(irrig, dtype=float)
    if use_irrigation is None:
        use_irrigation = nutrient.upper().startswith("TN") \
            or nutrient.upper() in ("N", "N-NO3", "NO3")
    return precip + irrig if use_irrigation else precip.copy()


def soil_reduction_x(water):
    """x of Eq. (2): inverse water input, maximum-scaled.

    x_i = (1 / water_i) / max_j(1 / water_j) over the values given
    (the caller passes one year's cells, or all cells and years for
    the global scope). Cells with water <= 0 would have an unbounded
    inverse; they are assigned x = 1 (the maximum-reduction end) and
    excluded from the maximum, a rule the paper does not discuss
    because Mediterranean precipitation is positive everywhere.
    """
    water = np.asarray(water, dtype=float)
    x = np.ones_like(water)
    pos = water > 0
    if pos.any():
        inv = np.where(pos, 1.0 / np.where(pos, water, 1.0), 0.0)
        x = np.where(pos, inv / inv.max(), 1.0)
    return x


def soil_reduction_factor(precip, irrig, alpha, nutrient="TN",
                          use_irrigation=None):
    """S of Eq. (2): exp(-alpha * x) for one set of cells (one year,
    or whatever scope the caller passes -- the normalisation maximum
    is taken over exactly the values given)."""
    water = soil_water_input(precip, irrig, nutrient, use_irrigation)
    return np.exp(-alpha * soil_reduction_x(water))


def river_reduction_factor(spec_runoff, A, B):
    """R of Eq. (3): 1 / (1 + A * qs**B), qs in l/s/km2."""
    qs = np.asarray(spec_runoff, dtype=float)
    if np.isnan(qs).any():
        raise ValueError("SpecRunoff contains NaN; the paper derives "
                         "qs for every cell from modelled discharge, "
                         "so a missing value is a data gap, not a "
                         "default.")
    return 1.0 / (1.0 + A * qs ** B)


def lake_retention_factor(hyd_load, C, n):
    """RES of Eq. (4): 1 / (1 + n * HL**C) on lake cells.

    hyd_load NaN marks a non-lake cell and returns RES = 0 there
    (Eq. 4 at HL = 0 would return 1, so the formula is meaningful
    only where a lake exists). n is the fixed coefficient the paper
    attributes to Behrendt et al. (2003); its value is not printed
    in Malagó et al. (2019) and is a required argument here.
    """
    hl = np.asarray(hyd_load, dtype=float)
    res = np.zeros_like(hl)
    lake = ~np.isnan(hl)
    res[lake] = 1.0 / (1.0 + n * hl[lake] ** C)
    return res


# --------------------------------------------------------------------------
# Core Rgrid computation
# --------------------------------------------------------------------------

def _prepare_rgrid(catch_data, annual_data, surplus_components):
    missing = [c for c in ("HydroID", "To_catch", "LakeFrRet",
                           "NrmLengthKm") if c not in catch_data.columns]
    if missing:
        raise ValueError(f"catch_data is missing columns: {missing}")
    catch = catch_data.copy()
    if "Shreve" not in catch.columns:
        catch = _catchment.compute_shreve_levels(catch)
    need = ["HydroID", "YearValue", "Ps"] + list(surplus_components)
    missing = [c for c in need if c not in annual_data.columns]
    if missing:
        raise ValueError(f"annual_data is missing columns: {missing}")
    return catch, annual_data.copy()


def launch_green_rgrid(catch_data, annual_data, loc_years, *,
                       forms="malago2019", nutrient="TN",
                       surplus_components=("Surplus",),
                       alpha=None, A=None, B=None, C=None, n=None,
                       alpha_p=None, alpha_l=None,
                       use_irrigation=None, x_scope="per_year"):
    """Run GREEN-Rgrid (or the Rgrid routing with catchment forms).

    forms="malago2019": factors from Eqs. (2)-(4); alpha, A, B, C
    and n are then REQUIRED (they are the paper's calibration
    parameters / Behrendt coefficient; this module sets no defaults
    for them). x_scope="per_year" normalises x by each year's own
    spatial maximum (the paper normalises over the study area);
    x_scope="global" uses one maximum over all years supplied.

    forms="catchment": S = exp(-alpha_p * InvNrmRain) and
    R = exp(-alpha_l * NrmLengthKm), RES = LakeFrRet -- the GREENeR
    factors -- fed through Eq. (1)'s arrangement. With "Surplus"
    carrying GREENeR's SumDS and "Ps" its SumPSCmp1 this reproduces
    green.launch_green exactly (the degeneracy gate), which is its
    purpose; production Rgrid runs use forms="malago2019".

    Returns one row per cell and year: HydroID, To_catch, Year,
    Mode, S, R, RES, X (Eq. 2's x; NaN in catchment mode), Surplus
    (sum of components), Ps, UpstreamLoad, CellLoad, and ObsLoad
    where the input carried YearlyMass or ObsLoad.
    """
    if forms not in ("malago2019", "catchment"):
        raise ValueError(f"Unknown forms: {forms!r}")
    if forms == "malago2019":
        missing_pars = [name for name, val in
                        (("alpha", alpha), ("A", A), ("B", B),
                         ("C", C), ("n", n)) if val is None]
        if missing_pars:
            raise ValueError(
                "forms='malago2019' requires the paper's calibration "
                f"parameters {missing_pars}: Malagó et al. (2019) "
                "leave alpha, A, B, C to calibration (their Table 3 "
                "lists Mediterranean optima) and n to Behrendt et "
                "al. (2003); this module chooses none of them.")
        need = ["Precip", "Irrig", "SpecRunoff", "HydLoad"]
    else:
        if alpha_p is None or alpha_l is None:
            raise ValueError("forms='catchment' requires alpha_p and "
                             "alpha_l (the GREENeR parameters).")
        need = ["InvNrmRain"]
    catch, annual = _prepare_rgrid(catch_data, annual_data,
                                   surplus_components)
    missing = [c for c in need if c not in annual.columns]
    if missing:
        raise ValueError(f"annual_data is missing columns: {missing}")
    years = _catchment._check_years(
        loc_years, annual["YearValue"].unique())

    pos = {h: i for i, h in enumerate(catch["HydroID"].to_numpy())}
    order, upstreams = _catchment._topological_order(catch)
    lake_fr = catch["LakeFrRet"].to_numpy(dtype=float)
    nrm_len = catch["NrmLengthKm"].to_numpy(dtype=float)

    def cell_frame(sub, yr):
        """Per-cell factor and source arrays for one year."""
        idx = sub["HydroID"].map(pos).to_numpy()
        sur = np.zeros(len(catch))
        ps = np.zeros(len(catch))
        sur[idx] = sub[list(surplus_components)].sum(axis=1).to_numpy(
            dtype=float)
        ps[idx] = sub["Ps"].to_numpy(dtype=float)
        if forms == "malago2019":
            water = np.zeros(len(catch))
            water[idx] = soil_water_input(
                sub["Precip"].to_numpy(dtype=float),
                sub["Irrig"].to_numpy(dtype=float),
                nutrient, use_irrigation)
            qs = np.zeros(len(catch))
            qs[idx] = sub["SpecRunoff"].to_numpy(dtype=float)
            hl = np.full(len(catch), np.nan)
            hl[idx] = sub["HydLoad"].to_numpy(dtype=float)
            rr = river_reduction_factor(qs, A, B)
            res = lake_retention_factor(hl, C, n)
            if x_scope == "global":
                # x filled in by the caller, who owns the global max.
                xx = np.full(len(catch), np.nan)
            else:
                xx = np.zeros(len(catch))
                xx[idx] = soil_reduction_x(water[idx])
            ss = np.exp(-alpha * xx) if x_scope != "global" else xx
        else:
            inv_rain = np.zeros(len(catch))
            inv_rain[idx] = sub["InvNrmRain"].to_numpy(dtype=float)
            ss = np.exp(-alpha_p * inv_rain)
            rr = np.exp(-alpha_l * nrm_len)
            res = lake_fr.copy()
            xx = np.full(len(catch), np.nan)
        return sur, ps, ss, rr, res, xx

    frames = []
    if forms == "malago2019" and x_scope == "global":
        waters = {}
        for yr in years:
            sub = annual[annual["YearValue"] == yr]
            idx = sub["HydroID"].map(pos).to_numpy()
            w = np.zeros(len(catch))
            w[idx] = soil_water_input(
                sub["Precip"].to_numpy(dtype=float),
                sub["Irrig"].to_numpy(dtype=float),
                nutrient, use_irrigation)
            waters[yr] = (idx, w)
        allw = np.concatenate([w[idx] for idx, w in waters.values()])
        gmax_inv = (1.0 / allw[allw > 0]).max()
    for yr in years:
        sub = annual[annual["YearValue"] == yr]
        sur, ps, ss, rr, res, xx = cell_frame(sub, yr)
        if forms == "malago2019" and x_scope == "global":
            idx, w = waters[yr]
            xx = np.where(w > 0, (1.0 / np.where(w > 0, w, 1.0))
                          / gmax_inv, 1.0)
            ss = np.exp(-alpha * xx)
        load = np.zeros(len(catch))
        up_load = np.zeros(len(catch))
        for i in order:
            ul = load[upstreams[i]].sum() if upstreams[i] else 0.0
            up_load[i] = ul
            # Eq. (1): the upstream load joins the point-source term:
            # reduced by R and RES, not by S.
            load[i] = (sur[i] * ss[i] + ps[i] + ul) * rr[i] \
                * (1.0 - res[i])
        frame = catch[["HydroID", "To_catch"]].copy()
        frame["Year"] = yr
        frame["Mode"] = MODE_MALAGO if forms == "malago2019" \
            else MODE_CATCHMENT_FORMS
        frame["S"] = ss
        frame["R"] = rr
        frame["RES"] = res
        frame["X"] = xx
        frame["Surplus"] = sur
        frame["Ps"] = ps
        frame["UpstreamLoad"] = up_load
        frame["CellLoad"] = load
        obs_col = "YearlyMass" if "YearlyMass" in sub.columns else \
            ("ObsLoad" if "ObsLoad" in sub.columns else None)
        if obs_col is not None:
            obs = np.full(len(catch), np.nan)
            idx = sub["HydroID"].map(pos).to_numpy()
            obs[idx] = sub[obs_col].to_numpy(dtype=float)
            frame["ObsLoad"] = obs
        frames.append(frame)
    out = pd.concat(frames, ignore_index=True)
    return out.sort_values(["Year", "HydroID"]).reset_index(drop=True)


def green_rgrid_shares(catch_data, annual_data, loc_years, **kwargs):
    """Source apportionment for GREEN-Rgrid, mirroring green_shares.

    Runs launch_green_rgrid once per source (each surplus component
    and Ps) with all other sources zeroed -- valid because Eq. (1)
    is linear in (SUR, PS) at fixed S, R, RES -- and returns, per
    cell and year, the routed load attributable to each source plus
    their sum, CellLoad. Accepts the same keyword arguments as
    launch_green_rgrid; surplus_components defaults to ("Surplus",).
    """
    components = list(kwargs.get("surplus_components", ("Surplus",)))
    sources = components + ["Ps"]
    frames = []
    for src in sources:
        mod = annual_data.copy()
        for other in sources:
            if other != src:
                mod[other] = 0.0
        res = launch_green_rgrid(catch_data, mod, loc_years, **kwargs)
        res = res[["HydroID", "To_catch", "Year", "Mode",
                   "CellLoad"]].rename(columns={"CellLoad": src})
        frames.append(res.set_index(["HydroID", "To_catch", "Year",
                                     "Mode"]))
    out = pd.concat(frames, axis=1).reset_index()
    out["CellLoad"] = out[sources].sum(axis=1)
    return out.sort_values(["Year", "HydroID"]).reset_index(drop=True)
