# GREEN-Rgrid adaptation notes — Malagó et al. (2019) forms

**Paper of record:** Malagó, Bouraoui, Grizzetti & De Roo (2019),
"Modelling nutrient fluxes into the Mediterranean Sea with
GREEN-Rgrid", *J. Hydrology: Regional Studies* 22, 100592
(PMC6472532 is the same paper). **Ruling implemented:** Alex,
2026-11-21 — adapt GREEN-Rgrid with Malagó's published forms as
written. **Code:** `green_rgrid.py` in this folder, beside the
verified GREENeR catchment port `green.py`, which is **unmodified**
(md5 `b10fc7266981d40f2032d4d6ac08b68d`; its `verify_port.py` was
re-run after this work and reproduced `verification_results.json`
byte-identically). **Verification:** `verify_rgrid.py`, numbers in
`adaptation_verification_results.json`.

The adaptation is an additional mode of the GREEN member, not a
replacement. The switch is the `forms` argument of
`launch_green_rgrid`: `forms="malago2019"` for production Rgrid
runs, `forms="catchment"` for the degeneracy configuration. Every
output row carries a `Mode` column naming the form set used, and
the Rgrid routing reuses `green.py`'s topological ordering, so both
modes accumulate downstream through the same ordering code.

## The model as published (Eqs. 1–4)

Per grid cell *i* and year, with loads in ton/y:

```
L_i = [SUR_i · S_i · R_i + (PS_i + UL_i) · R_i] · (1 − RES_i)     (1)
S_i   = exp(−α · x_i)                                             (2)
R_i   = 1 / (1 + A · qs_i^B)                                      (3)
RES_i = 1 / (1 + n · HL_i^C)                                      (4)
```

`SUR` is the cell's nutrient surplus (the diffuse source term),
`PS` the point sources routed as direct inputs, and `UL` the
upstream load: the paper (§2.2) adds the upstream load to the
receiving cell "as an additional point source", so `UL` is reduced
by `R_i` and `RES_i` but **not** by the soil factor `S_i`. Gate 3
includes a two-cell chain test that discriminates exactly this
placement.

## The two form differences from GREENeR's catchment GREEN

**1. Soil reduction factor — precipitation plus irrigation.**
GREENeR: `AA = exp(−alpha_P · InvNrmRain)` with
`InvNrmRain = 50 / max(50, Rain)` (rainfall alone, floored
normalisation). Malagó Eq. (2): `S = exp(−α · x)` where `x` is the
inverse of the water input normalised by the **maximum inverse
over the whole study area** (Frank & Todeschini 1994 maximum
scaling). One refinement the task summary did not carry: the
paper's input section (§2.3.2) states the model "uses annual
precipitation and irrigation as proxies for the basin retention of
**nitrogen**, and **only precipitation for phosphorus**". The
implementation follows the paper: water input = precipitation +
irrigation for TN (and N-NO3), precipitation alone for TP (and
P-PO4), with an explicit `use_irrigation` override for callers who
need the other reading.

**2. River retention — specific runoff, not reach length.**
GREENeR: `BB = exp(−alpha_L · NrmLengthKm)`. Malagó Eq. (3):
`R = 1 / (1 + A · qs^B)`, a rational function of the average area
specific runoff `qs` (l/s/km²; in the paper, LISFLOOD annual
discharge divided by the cell's upstream area), with two
calibration parameters `A` and `B`.

**Lake retention** also changes form, though it was not one of the
two flagged differences: GREENeR carries a per-catchment
`LakeFrRet` fraction (Kronvang formula in `read_nsdata`), while
Malagó Eq. (4) computes `RES = 1 / (1 + n · HL^C)` from the
hydraulic load `HL` (m; lake outflow divided by lake area), with
calibration parameter `C` and a fixed coefficient `n`.

## Parameters: nothing silently chosen

`alpha`, `A`, `B`, `C` are the paper's calibration parameters and
`n` is "a fixed coefficient retrieved from Behrendt et al.
(2003)". **Malagó et al. do not print the value of `n` anywhere in
the paper** (checked in the PDF text and the PMC HTML full text),
so `launch_green_rgrid` takes all five as required arguments and
raises if any is missing. The paper's calibrated Mediterranean
optima (its Table 3) are published values a caller may adopt —
they are recorded here for reference, not used as defaults:

| Determinand | α | A | B | C | η |
|---|---|---|---|---|---|
| N-NO3 | 5.59 | 0.88 | 4.62 | 0.59 | 0.34 |
| TN | 8.39 | 6.70 | 2.95 | 0.96 | 0.39 |
| P-PO4 | 14.98 | 8.59 | 4.24 | 0.49 | NA |
| TP | 1.89 | 0.84 | 0.76 | 0.11 | NA |

`η` is a further calibration quantity of the paper's nitrogen
source construction (the natural-area fraction used when estimating
the atmospheric-deposition load in natural areas, §2.4); it enters
the surplus/deposition layers upstream of Eqs. (1)–(4), not the
routing core, and is therefore not a parameter of this module. It
will be needed when the corridor surplus grid is built.

## Interpretation rules the implementation had to state

- **Eq. (4) applies only to lake cells.** At `HL = 0` the formula
  returns `RES = 1` (total retention), so it cannot describe a
  cell without a lake. Following the paper's construction (lakes
  positioned from HydroLAKES, at most one per cell), a cell whose
  `HydLoad` is NaN is a non-lake cell and gets `RES = 0`; a lake
  cell's `HL` follows the formula as written.
- **`x` normalisation scope.** The paper normalises by the maximum
  inverse "over the entire Mediterranean area" — the study area,
  spatially. The precipitation and irrigation fields are annual,
  so the default (`x_scope="per_year"`) takes each year's own
  spatial maximum; `x_scope="global"` takes one maximum over all
  years supplied. The choice is explicit in the call and should be
  recorded with any run's outputs.
- **Zero water input.** A cell with precipitation (+ irrigation)
  of zero has an unbounded inverse; it is assigned `x = 1` (the
  maximum-reduction end) and excluded from the maximum. The paper
  does not discuss the case (Mediterranean precipitation is
  positive everywhere).
- **Missing specific runoff is an error**, not a default: the
  paper derives `qs` for every cell from modelled discharge, so a
  NaN `SpecRunoff` raises rather than silently implying `R = 1`.

## Grid topology builder (`build_grid_catchments`)

Inputs are a regular grid's D8 drainage, from either a DEM
(`d8_from_dem`, steepest descent) or an ESRI-coded flow-direction
array (`d8_from_codes`). Output is a `green.py` catchment table —
one row per active cell, `HydroID` = flat row-major cell index,
`To_catch` = the receiving cell (or −1) — plus the
`HydroID → (row, col)` map, so a grid feeds both GREEN modes
through the same data structures the catchment port consumes.
`NrmLengthKm` is the cell's flow length (cell size, or √2 × cell
size for a diagonal step; cell size for outlets) normalised by the
grid maximum, GREENeR's definition. `assemble_grid_annual` builds
one year of Rgrid `annual_data` from 2D field arrays and is the
assembly point for real rasters.

Builder rules and stated limits:

- **Single flow direction (D8)** — no flow dispersion, matching
  the paper's emitting–receiving cell structure.
- **Steepest descent** on drop per unit distance (diagonals
  weighted by √2); a strictly positive drop is required.
- **No depression filling and no flat resolution**: pits and flats
  with no lower neighbour become outlets (−1). Production routing
  should not use raw DEM descent; the paper's own route upscales
  30 s HydroSHEDS flow accumulation/direction to the 5 arc-min
  grid following Olivera et al. (2002), which is the recommended
  source of the direction grid (see application needs below).
- **Ties** in steepest descent are broken by a fixed neighbour
  scan order (north-west to south-east), stated in the code, so
  results are deterministic.
- Cells draining **off the grid** or **into an inactive
  (nodata) cell** become outlets (−1) and act as additional
  outlets of the study grid — the behaviour the PORT_NOTES grid
  recipe already prescribed for the catchment machinery.
- Cycles raise an error (via the shared Shreve-level computation)
  rather than routing silently.

## Verification gates (`verify_rgrid.py`) — all passing

1. **Topology builder.** DEM case (2×2, one inactive cell): built
   `To_catch` = {0→2, 1→2, 2→−1} exactly as intended; normalised
   flow lengths {1/√2, 1, 1/√2} to 1e-15; Shreve levels {1, 1, 2}.
   ESRI-codes case (3×3, all cells draining to one outlet): built
   `To_catch` = {0→4, 1→4, 2→4, 3→4, 4→7, 5→4, 6→7, 7→−1, 8→7}
   exactly.
2. **Degeneracy.** A grid whose cells coincide one-to-one with a
   catchment network, run through the Rgrid routing with
   `forms="catchment"` (GREENeR factors in Eq. (1)'s arrangement),
   reproduces `green.launch_green` exactly: branching TP case with
   lake retention, max abs difference **0.0**; four-catchment TN
   case over two years, max abs difference **1.8e-15**. The
   builder supplied the topology; per-cell `NrmLengthKm` /
   `LakeFrRet` were set from the catchment table under the cell
   mapping, since the gate is about routing and accumulation.
3. **Hand-computed Malagó forms** (expected values from the
   paper's equations via `math`, independent of the module):
   soil factor with irrigation, TN (P=[500, 250], I=[300, 0],
   α=2 → x=[0.3125, 1], S=[e^−0.625, e^−2]) — abs diff **0.0**;
   soil factor TP, precipitation only (x=[0.5, 1]) — **0.0**;
   river factor Eq. (3) at qs=10, A=0.84, B=0.76 → R=0.171418 —
   **0.0**; lake factor Eq. (4) at HL=10, C=0.11, n=0.05 →
   RES=0.939485, and NaN → 0 — **0.0**; full Eq. (1) on one cell
   (SUR=120, PS=30, α=1.89, qs=5, HL=8) — **0.0**; upstream-load
   placement on a two-cell chain (expected downstream load
   150·e^−1 + 10 = 65.181916; the soil-reduced-upstream variant
   differs by > 1, so the test discriminates) — **0.0**.
4. **Source-apportionment identity under the Rgrid forms.**
   Six-cell branching grid, two surplus components plus point
   sources, a lake cell, two years (12 cell-years): per-source
   routed loads sum to the total with max abs difference **0.0**;
   the shares total agrees with `launch_green_rgrid` to 7.1e-15.
5. **Catchment path unchanged.** `green.py`'s own hand case from
   `verify_port.py` reproduces its expected value (abs diff
   3.6e-15, the same as the port's recorded run), `green.py` is
   byte-unmodified, and the port's full `verify_port.py` re-run
   completed with `verification_results.json` byte-identical.

## What a corridor application still needs

- **Irrigation field — the paper's own source chain.** Malagó
  et al. spatialise **AQUASTAT (FAO 2016) country-level net
  irrigation volumes** onto the 5 arc-min grid using the
  **SPAM (You et al. 2014) irrigated-area** distribution weighted
  by the PET − precipitation deficit, with PET from annual
  temperature by the Langbein (1949) formula (their Eq. 5). A
  corridor run needs the Canadian AQUASTAT net irrigation volume,
  a SPAM (or successor) irrigated-area grid, and a temperature
  field — none is currently held in the project. Precipitation in
  the paper is **MSWEP** (Beck et al. 2017); the project's held
  climate grids can substitute, but the substitution must be
  recorded, since `x` is a ratio of inverses and is sensitive to
  the precipitation field's dry-end values.
- **Specific runoff and hydraulic load.** The paper derives both
  from **LISFLOOD** annual discharge (2000–2014): `qs` = annual
  discharge / upstream area per cell; `HL` = lake outflow / lake
  area, with lake positions and areas from **HydroLAKES**
  (Messager et al. 2016), at most one lake per cell. The corridor
  has no held discharge grid; HYDAT-gauge-based runoff surfaces or
  a hydrological model's gridded discharge are the candidate
  substitutes, and whichever is chosen defines `qs`'s upstream
  areas through the same grid topology the builder produces.
- **Flow direction at the run resolution.** The project holds no
  corridor DEM raster. The paper's route — HydroSHEDS 30 s flow
  accumulation/direction upscaled to the run grid following
  Olivera et al. (2002) — is the documented path; the builder
  consumes the resulting direction grid (or a DEM, with the
  pit/flat limits stated above). The run resolution itself
  (paper: 5 arc-min) is a project decision still owed.
- **Surplus grid construction from held inputs.** Per the paper's
  §2.3.3.2: cropland/fodder surplus = mineral + manure fertiliser
  + N fixation (N only) + atmospheric deposition (N only) − crop
  uptake; grassland surplus = inputs − 60% of inputs (Bouwman et
  al. 2005); plus natural-area background deposition, erosion,
  urban washoff, P weathering (GLIM lithology, Hartmann et al.
  2014) and scattered dwellings in the diffuse pool; point sources
  from connected population × per-capita emissions × treatment
  efficiency. Held material that maps onto this: the SPARROW
  deposit's per-reach source columns (fertiliser P, manure P,
  point-source P) for aggregation to cells; the census/IPNI
  livestock and fertiliser trajectories built in this project for
  the temporal dimension; the held WorldCover cache
  (`maps/worldcover_cache.npz`) for land-cover classes and forest
  fractions; SoilGrids (via the maps builder's VRT route) for
  soil properties. Not yet held: crop uptake grids (SPAM-based in
  the paper), the GLIM lithology layer, and HydroLAKES for the
  corridor.
- **Calibration.** A corridor calibration needs observed annual
  loads at grid outlets (the project's station load estimates are
  the candidate set) and must choose values for α, A, B, C — and
  a sourced value of `n` from Behrendt et al. (2003), which
  Malagó et al. cite but do not print.
