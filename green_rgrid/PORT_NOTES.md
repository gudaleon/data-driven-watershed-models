# GREEN port notes — GREENeR (R, GPL-3) → Python

**Source of record:** GREENeR 1.0.2, cloned at
`~/workspace/fei_lake_ontario/codebase/GREENeR` (core file `R/rgreen.R`),
example data at `~/workspace/fei_lake_ontario/codebase/GREENeR_data`.
**Port:** `green.py` in this folder (numpy/pandas, runs under
`~/workspace/mapvenv/bin/python`). **Verification:** `verify_port.py`,
numbers in `verification_results.json`. Licence: the port is a
derivative of GPL-3 code and is GPL-3 accordingly; GREENeR remains
© its authors (Udias, Grizzetti, Vigiak, Gomez, Alfaro, Aloe, JRC).

## What was ported

| GREENeR function | Port | Notes |
|---|---|---|
| `launch_green` | `launch_green` | The core load computation (below). |
| `aggregate_loop` | `_route` | Downstream accumulation. The port routes in a true topological (Kahn) order, which reproduces the R Shreve-grouped loop exactly for valid inputs and does not depend on which Shreve convention the input carries. |
| `green_shares` | `green_shares` | Source apportionment: one run per source with all other sources zeroed; valid because GREEN is linear in the sources at fixed retention factors. |
| `region_nut_balance` | `region_nut_balance` | Diffuse / land-retention decomposition plus routed river and lake retention (one faithfully replicated quirk — see below). |
| `read_NSdata` | `read_nsdata` | Builds catch/annual tables from the GREENeR_data CSV layout, including the Kronvang lake-retention formula (7.3 for TN, 26 for TP; capped at 0.1), length normalisation, and `InvNrmRain = 50 / max(50, Rain)`. |
| `shreve` | `compute_shreve_levels` | See the Shreve caveat below. |
| `check_colnames_*`, `check_years`, `append_empty_cols`, `data_preparation` | folded into `_prepare` / `_check_years` | Same column schemas, same TP handling. |

## The core computation (per catchment *i*, per year)

```
AA_i           = exp(-alpha_P * InvNrmRain_i)              # 1 - basin (land) retention
SumDS_i        = Min + Man + Fix + Soil + Atm * (1 - ForestFraction)
SumPSCmp1_i    = Ps + atm_coeff * Atm * ForestFraction + sd_coef * Sd
CatchToRiver_i = SumDS_i * AA_i + SumPSCmp1_i
BB_i           = exp(-alpha_L * NrmLengthKm_i)             # 1 - river retention
CatchLoad_i    = (1 - LakeFrRet_i) * (CatchToRiver_i + U_i) * BB_i
```

`U_i` is the sum of `CatchLoad` over the catchments draining directly
into *i*. Parameters: `alpha_P` (basin retention coefficient, driven by
the inverse normalised rainfall), `alpha_L` (river retention
coefficient, driven by the normalised reach length), `sd_coef`
(fraction of scattered-dwelling input reaching the stream), and
`atm_coeff = 0.38` (fixed in GREENeR). Diffuse agricultural sources
(`Min, Man, Fix, Soil`, and the agricultural fraction of atmospheric
deposition) are reduced by the land retention `AA`; point sources `Ps`
enter the river directly; lake retention `LakeFrRet` and river
retention `BB` then apply to everything passing through the catchment.

## Faithfulness notes and quirks found in the R code

1. **TP background split.** GREENeR handles TP by renaming `Bg` to
   `Atm` and setting `Fix = Soil = 0` (`append_empty_cols`), so in the
   code the TP background term is split as
   `Bg*(1-FF)` through land retention plus `0.38 * Bg * FF` direct.
   The vignette's printed TP equation shows the forest-fraction term
   without the 0.38 factor. **The port follows the code**, and any
   future Stan re-expression should decide explicitly which form it
   implements.
2. **Balance retention terms.** In `aggregate_loop` with `type = "N"`,
   the data.table statements execute sequentially, so `CatchRivRet`
   and `CatchLakeRet` are computed from the *already updated*
   `CatchLoad` rather than from the pre-retention flux. The port
   replicates this ordering and flags it here; it affects only the
   diagnostic balance columns, not `CatchLoad`.
3. **"Shreve" is two different things.** `read_NSdata`'s `shreve()`
   generates topological levels (headwaters = 1, each next generation
   = max upstream level + 1). The package's built-in `catch_data_TN`
   instead carries classic Shreve magnitudes (sums). Both are valid
   processing orders for `aggregate_loop`; the port's Kahn routing is
   equivalent to both, and `compute_shreve_levels` reproduces the
   `shreve()` level definition for inputs that lack the column.
4. **GREENeR_data is not the provenance of the built-in example.**
   The annual source data built from the GREENeR_data CSVs matches
   `annual_data_TN` exactly (below), but the length file does not
   correspond to the built-in `NrmLengthKm` at all (correlation
   −0.04), and the built-in `LakeFrRet` differs from the Kronvang
   formula on the CSV lake properties by up to 0.006. Anyone using
   `read_NSdata` (R or Python) on those CSVs gets a *different
   catchment table* from the package's built-in one. This is an
   inconsistency inside the upstream repositories, recorded here so it
   is not mistaken for a port error.

## Verification (verify_port.py, verification_results.json)

The direct gate — run GREENeR in R on an example basin and compare
loads number for number — **could not be executed**: this VM currently
has no R interpreter (`which Rscript` is empty; `~/R/library` survives
but R itself was wiped by a VM restart, and standing project rules say
the agent does not reinstall R). The gate should be re-run when R is
available; the script structure anticipates adding that comparison to
`verification_results.json`. Substitute verification, all passing:

- **Hand-computed tests** (expected values computed in the test script
  directly from the equations with `math.exp`): single-catchment TN
  case, absolute difference 3.6e-15; three-catchment branching TP case
  with lake retention, max absolute difference 1.8e-15.
- **Input equivalence:** `read_nsdata` on the GREENeR_data CSVs vs the
  package's built-in Lay TN data — topology (`To_catch`) mismatches: 0
  of 189; all seven source columns differ by exactly 0.0 over 5,481
  matched catchment-years; `ForestFraction` max difference 5.0e-10;
  `InvNrmRain` max difference 2.4e-11. (Catchment scalars differ for
  the provenance reasons in note 4.)
- **Source-share identity** (an exact algebraic property of GREEN):
  per-source routed loads sum to the total load with max absolute
  difference **0.0** on all 1,323 catchment-years of the Lay TN run
  and all 1,323 of the Lay TP run; shares vs `launch_green` totals
  agree to 9.1e-13 (floating-point summation order).
- **Plausibility against observations:** Lay TN, 2003–2009, with the
  vignette's *uncalibrated* documentation parameters
  (alpha_P = 35.09, alpha_L = 0.02, sd_coef = 0.2): modelled vs the 22
  observed loads in the data give Pearson r = 0.990, RMSE 656 t/yr,
  mean bias −428 t/yr; outlet (HydroID 368447) loads 1,161–3,207 t/yr
  across 2003–2009, the same magnitude as the observed downstream
  loads. This is a plausibility check on an arbitrary parameter set,
  not a calibration claim.

## Deliberately left out

- `calib_green` / `calib_green_help` — the Latin-Hypercube calibration.
  Per the ensemble design, calibration of this member will be Bayesian
  (Stan) rather than LHS-plus-goodness-of-fit; the port's
  `launch_green` is written so its per-catchment terms map directly
  onto a Stan data block (arrays of sources, `AA`/`BB` factors, and a
  downstream index for the routing loop).
- `gof` goodness-of-fit metrics and `optimization_metrics.R`.
- All plotting and mapping functions (`rgreen_plots.R`, ~75 KB of the
  package): input plots, maps, Sankey diagrams, calibration graphics.
- The `sf` geometry reader (`read_geometry`) — the port is table-only;
  geometry joins belong to the project's Python GIS layer.

## Feeding a regular grid to this port (the GREEN-Rgrid adaptation)

Malagó et al. (2019) ran GREEN on a regular 5 arc-min grid; no code
for that variant was ever released, and GREENeR's catchment machinery
accepts it without modification once the grid is expressed as a
catchment table:

- **One row per grid cell** in `catch_data`: `HydroID` = cell index;
  `To_catch` = the cell's D8 downstream neighbour (the outlet cell of
  the whole grid gets −1; cells draining outside the study grid also
  get −1 and act as additional outlets).
- **`Shreve`** may be omitted — `compute_shreve_levels` derives the
  topological levels from `To_catch`, and routing uses Kahn order
  regardless.
- **`NrmLengthKm`** = the cell's flow length (cell size, or √2 × cell
  size for diagonal D8 steps) normalised by the basin maximum, exactly
  GREENeR's definition; alternatively the raw length can be carried
  and `alpha_L` rescaled, but the normalised form keeps GREENeR's
  parameter ranges meaningful.
- **`LakeFrRet`** per cell from a lake/reservoir raster (HydroLAKES
  depth and residence time through the same Kronvang formula
  `read_nsdata` implements; zero for non-lake cells).
- **Annual sources per cell** in `annual_data`: fertiliser and manure
  rates (kg/ha/yr × cell area), background from land cover, point
  sources assigned to the cell containing the discharge, scattered
  dwellings from population grids; `ForestFraction` from the
  land-cover grid; `InvNrmRain` from the precipitation grid with
  GREENeR's normalisation (inverse rainfall scaled so the wettest
  cell-year maps consistently — for multi-year runs, fix the
  normaliser across years so `alpha_P` stays comparable).

Two genuine differences from Malagó et al.'s Rgrid formulation must be
decided before the adaptation is called GREEN-Rgrid rather than
"GREEN on a grid": their soil reduction factor is driven by
precipitation **plus irrigation**, where GREENeR's `AA` uses rainfall
alone; and their river retention is a function of **specific runoff**,
where GREENeR's `BB` is a function of reach length. Both are one-line
changes in `_local_fluxes` / the `BB` computation once the project
rules on them.
