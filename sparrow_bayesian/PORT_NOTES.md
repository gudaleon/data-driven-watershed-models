# RSPARROW → Python port notes (classical NLLS core)

Source: USGS RSPARROW 1.1.0 (public domain; Alexander & Gorman Sanisaca,
USGS software release, DOI 10.5066/P9UAZ6FO), cloned at
`~/workspace/fei_lake_ontario/codebase/RSPARROW`.
Port: `classical_sparrow.py` (numpy / scipy / pandas; mapvenv).
The `classical_` prefix leaves the folder free for the future Stan
Bayesian version, which will share this module's data structures.

## What was ported

The model-defining core only, isolated from these R/Fortran files:

| RSPARROW source | Ported as |
|---|---|
| `R/estimateFeval.R` | `ClassicalSparrow.delivery_factors`, `.decay_factors`, `.incremental_loads`, `.residuals` |
| `src/tnoder.for` | `ClassicalSparrow._accumulate` with `condition=True` |
| `src/ptnoder.for` (via `R/predict.R`) | `._accumulate` with `condition=False`; `predict_loads` |
| `R/estimateOptimize.R` (NLMRT `nlfb`) | `.fit` using `scipy.optimize.least_squares`, bounds per parameters.csv with estimateOptimize.R's bound-relaxation rule |
| `R/estimateNLLSmetrics.R` | SEs from `MSE · (JᵀJ)⁻¹` at the solution; leverage-adjusted retransformation factor |
| `R/createDataMatrix.R` | beta layout `[SOURCE, DELIVF, STRM, RESV]`; delivery variables mean-centred over all reaches when `mean_adjust_delivery` is on |
| `R/selectParmValues.R`, `R/readParameters.R` | `Parameter.is_selected` / `.is_constant` (parmConstant ≡ init = min = max ≠ 0) |
| `R/setNLLSWeights.R` | default weights = 1 (`NLLS_weights="default"`); per-site weights accepted by the constructor |

Model equations (identical to estimateFeval.R / predict.R):

- delivery factor per source: `F[i,s] = exp(Σ_d design[s,d] · β_dlv[d] · D_adj[i,d])`
- delivered source load: `G[i] = Σ_s β_src[s] · S[i,s] · F[i,s]`
- stream decay `A[i] = Π_k exp(−T[i,k]·β_str[k])`; reservoir decay
  `B[i] = Π_m 1/(1 + R[i,m]·β_res[m])`
- incremental load `inc[i] = √A[i]·B[i]·G[i]` (the √ applies in-reach
  decay to half the reach for internally generated load);
  transport factor `carry[i] = frac[i]·A[i]·B[i]`
- accumulation in hydrologic order:
  `load[i] = inc[i] + carry[i]·node[fnode[i]]`, then
  `node[tnode[i]] += iftran[i]·load[i]`
- calibration (tnoder, ifadjust = 1): at a reach with observed load
  `depvar > 0`, residual `= log(depvar) − log(load)` (a non-positive
  predicted load is replaced by 1.0 before the log, as in the
  Fortran), and the observed load is substituted downstream
  ("monitoring adjustment" / conditioning). Residuals are weighted by
  √weight.

## What was left out (deliberately)

- The Shiny Decision Support System, interactive mapper and all
  plotting/report code (`predictMaps`, `diagnosticPlots*`, `goShiny*`,
  the Rmd diagnostics).
- Scenario execution (`predictScenarios*`) — the AST/BMP layer is a
  separate ensemble member; note the scenario engine calls the same
  accumulation ported here.
- Bootstrap estimation (`estimateBootstraps`, `predictBootstraps`) —
  the port computes the one bootstrap-derived quantity that enters
  published predictions, the retransformation factor (below).
- Data-prep and GIS utilities (hydseq generation, drainage-area
  checks, shapefile handling). The port consumes a prepared reach
  table; if no `hydseq` column is supplied it computes a valid
  hydrologic order itself (`topological_order`, Kahn's algorithm).

## Verification

Two independent gates; scripts `test_classical_sparrow.py` and
`verify_classical_sparrow.py`. (A live R re-run was not used as the
gate: RSPARROW ships Windows `.dll` binaries for its Fortran routines
and has a large CRAN dependency tree, so a full Linux re-install was
impractical in this environment. The gates below compare against
RSPARROW's own published tutorial outputs — the R original's record —
and against hand computation.)

### Gate 1 — synthetic network, hand-computed (8 reaches)

Two sources (one delivery-adjusted), one delivery variable, one
stream-decay and one reservoir-decay variable, a 0.75/0.25 flow split,
and one monitoring site. Every expected load was computed by hand in
the test file with plain Python arithmetic:

- unconditioned loads: max relative error **0.0**
- conditioned loads: max relative error **0.0**
- log residual at the site: max relative error **0.0**
- NLLS recovery with all 8 reaches observed: SSE = 7.9e-31, all five
  coefficients recovered exactly. (With only 3 of 8 reaches observed
  the 5-parameter system is underdetermined — a manifold of exact
  fits — so the recovery gate observes all reaches.)

### Gate 2 — published tutorial Model8 (MRB3 total nitrogen)

11,526 reaches, 708 calibration sites, 13 estimated parameters
(5 SOURCE, 5 DELIVF, 2 STRM, 1 RESV), delivery mean adjustment on.
Python inputs rebuilt from `UserTutorial/data/data1.csv`
(depvar = `LOAD_A_00600` at `Tagsite == 1` reaches — exactly 708
sites, matching the summary) and `Model8_design_matrix.csv`.

- **SSE at the published estimates: Python 115.6632 vs RSPARROW
  115.6632** (Model8_summary.txt) — identical at printed precision.
- Per-station log residuals vs `Model8_residuals.csv`: max abs
  difference **5.9e-04** over 708 sites (consistent with the summary
  printing coefficients to 9 significant digits).
- Per-reach loads vs `Model8_predicts_load.csv`: the published file
  is the port's accumulation times one constant,
  ratio min 1.08253281 / max 1.08253348 across all 11,526 reaches.
  That constant is RSPARROW's documented retransformation factor
  `bootcorrection = mean(exp(Resids / √(1 − leverage)))`
  (predict.R:413, estimateNLLSmetrics.R:313), applied to every
  published prediction. Recomputing the factor from the port's own
  residuals and Jacobian gives 1.083199, and after applying it the
  loads match the published file to max relative difference
  **6.2e-04** for both `pload_total` and `mpload_total`.
- **NLLS refit** from the parameters.csv starting values
  (scipy `least_squares`): SSE 115.6632, and all 13 coefficients
  recovered to within 0.04 % of the published estimates (worst case
  rchdecay1: 0.42424 vs 0.42409). Standard errors agree with the
  published values to ~1–8 % (e.g. point 0.1125 vs 0.1110;
  ndep 0.0398 vs 0.0388), the small differences coming from the
  Jacobian evaluation method (NLMRT Richardson vs scipy 2-point).

## Mapping onto the Stan Bayesian version

The Fortran accumulation is already a sequential loop, so the Stan
model can mirror it one-to-one. `SparrowNetwork` precomputes, in
Python, everything Stan needs as data:

- `order` — reach indices in topological (upstream→downstream)
  order. In Stan, iterate `for (r in 1:n_reach)` over reaches in
  this order; no recursion or graph logic is needed inside Stan.
- `fnode_idx`, `tnode_idx` — node indices compressed to
  `1..n_node`. Stan keeps a `vector[n_node] node` accumulator,
  exactly the Fortran `node()` array:
  `load[r] = inc[r] + carry[r] * node[fnode[r]];`
  `node[tnode[r]] += iftran[r] * load[r];`
- `downstream_index` — for each reach (in topological order), the
  topological-order position of the first reach leaving its tnode,
  or −1 if terminal; `n_downstream` flags flow splits (> 1), where
  `frac` apportions the flux. This supports an alternative
  reach-vector formulation in Stan (load pushed downstream reach by
  reach) and generated-quantities bookkeeping (e.g. load delivered
  to a target reach).
- All parameter-dependent pieces (`inc`, `carry`) are elementwise
  vector operations on reach arrays — delivery factors are a matrix
  product `exp((D_adj .* beta_dlv) * design')`, decays are row
  products — so in Stan they are computed once per leapfrog step
  from the parameter vector, in the same beta order
  `[SOURCE, DELIVF, STRM, RESV]` used here.
- The likelihood replaces conditioning: the Bayesian version should
  use the *unconditioned* accumulation and score
  `log(depvar[site]) ~ normal(log(load[site]), sigma)` at the 708
  site reaches (site reach indices in topological order are a data
  array). Conditioning (substituting observed loads downstream) is
  an NLLS estimation device in tnoder; it has no place inside a
  generative Bayesian model, and dropping it is the one deliberate
  semantic difference between the classical and Bayesian members.
