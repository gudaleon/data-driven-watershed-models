# Bayesian SPARROW (Stan) — build notes

The Bayesian member of the Ontario SPARROW ensemble, built to the
model design settled on 2026-11-20. Files in this folder:

| File | What it is |
|---|---|
| `sparrow.stan` | The Stan model (CmdStan 2.40.0; compiles clean with the CmdStan `make`) |
| `prepare_stan_data.py` | Builds the Stan JSON from the project's built inputs (extended network + station loads v4) |
| `run_synthetic_test.py` | The synthetic recovery test — the build's verification gate |
| `synthetic/` | Synthetic data, Stan draw CSVs and `synthetic_results*.json` for two noise realisations |
| `sparrow` | Compiled model binary (build artefact; not source, do not commit) |
| `classical_sparrow.py`, `PORT_NOTES.md` | The verified classical RSPARROW port this build sits beside |

## Design trace

* **Sources per reach** (kg P/yr from the mid-continental data1 table,
  base year 2002): point-source P with its coefficient fixed at 1;
  fertilizer P and manure P, each multiplied by an estimated delivery
  coefficient; urban land area (km²) times an estimated export rate;
  background land (forest + grass/shrub + wetland + pasture area)
  times an estimated export rate; atmospheric deposition on water
  area times an estimated rate with the informative prior below.
  Agricultural land area is not a separate source term: fertilizer and
  manure already represent agricultural inputs.
* **Land-to-stream delivery**: fertilizer and manure share the factor
  exp(beta_runoff·z_runoff + beta_perm·z_perm + beta_tile·z_tile),
  with runoff, soil permeability and tile-drainage share standardised
  over reaches in `prepare_stan_data.py`. Expected signs (runoff +,
  permeability −, tile +) are left to the data; priors are centred at
  zero.
* **Attenuation**: stream reaches exp(−ks · travel time), travel time
  from reach length / mean velocity (velocity floored at 0.02 m/s —
  data1 carries zero velocities); reservoir reaches (rchtype 1) use
  1 / (1 + beta_res · rhload). In-reach generated load takes the
  half-reach factor sqrt(A) on stream reaches and B in full on
  reservoir reaches, exactly as the classical port's
  `incremental_loads`. `rhload` is taken from data1 (populated exactly
  on reservoir reaches there) and supplied divided by its median over
  the extended network's reservoir reaches (0.0041139), so the
  coefficient's prior is on a typical-reach scale.
* **Observation model**: one observation per station — its mean
  annual load over 1990–2021 — in log space, with per-station
  `obs_se` combining the load-estimator standard error of the mean
  and the interannual variance of the station's annual loads
  (computed in `prepare_stan_data.py`; zero-load years of one 0.1 km²
  station are excluded from the interannual component), plus a fitted
  residual SD `sigma` added in quadrature inside the model.
* **Generated quantities**: predicted loads at all 74 basin outlets
  and at the station reaches, plus pointwise log-likelihood.

### The unconditioned-accumulation decision (binding)

The classical member's calibration substitutes observed loads for
predicted loads at monitoring reaches before accumulating downstream
(tnoder conditioning). That is an NLLS estimation device, not a
generative step, so the Stan model accumulates from the parameters
alone and the observations enter only through the likelihood. This is
the one deliberate semantic difference between the two members, and
it is flagged in `PORT_NOTES.md` ("Mapping onto the Stan Bayesian
version") and in the `sparrow.stan` header. The accumulation itself
mirrors the Fortran loop one-to-one: reaches in topological order,
load pushed to the one or two reaches leaving each reach's to-node
(the extended network is dendritic except for a single flow-split
node), position N + 1 as a dummy sink for terminal reaches.

### Priors

| Parameter | Prior | Median |
|---|---|---|
| alpha_fert, alpha_man (delivery coefficients) | lognormal(log 0.3, 0.5) | 0.3 |
| rate_urban | lognormal(log 300, 0.5) | 300 kg/km²/yr |
| rate_bg | lognormal(log 25, 0.5) | 25 kg/km²/yr |
| rate_dep | lognormal(log 500, 0.25) | 500 kg/km²/yr (informative; 50 mg/m²/yr, Eimers et al. 2022) |
| ks (stream decay) | lognormal(log 0.15, 0.5) | 0.15 /day |
| beta_res (reservoir, scaled rhload) | lognormal(log 0.5, 0.75) | retention 0.67 at rhload = 1 |
| beta_runoff, beta_perm, beta_tile | normal(0, 0.5) | 0 |
| sigma (residual SD, log space) | normal(0, 0.5), truncated at 0 | — |

## Data preparation — verified end to end

`prepare_stan_data.py` runs on the real inputs and produces
`~/workspace/fei_lake_ontario/sparrow/sparrow_stan_data.json`:
**23,963 reaches, 181 stations fitted, 74 basin outlets**. The JSON
was validated against the compiled model with CmdStan's gradient
diagnostic (autodiff vs finite-difference gradients agree to ~1e-6),
which exercises the full data block and one log-density evaluation
over the whole network. Recorded exclusions and repairs, all stated
in the script:

* the 98 off-network station-year rows (the 15 CLOCA stations) are
  dropped — the source network has no reaches in that shore strip;
* 639 zero-incremental-area reaches in data1 have blank land-cover
  and manure fields; blank percents are treated as zero area in the
  class (the reaches still transport upstream load);
* station 03007703002's zero-load years (above).

## Synthetic recovery test — the verification gate

Protocol (`run_synthetic_test.py`): one real sub-basin (Bronte Creek,
outlet 940051310, 191 reaches, domain 1), known parameters, loads
generated by a numpy mirror of the Stan accumulation, observations at
the basin's 7 real load stations plus synthetic sites on larger
reaches (29–31 sites), log-space observation error with obs_se 0.08
and sigma 0.15. Two noise realisations (seeds 20261121, 20261122).

* **Mirror check**: the numpy data-generating mirror reproduces the
  verified classical port's unconditioned predictions to a maximum
  relative difference of **3.7e-16** — the synthetic data come from
  the verified code path, not from a second unchecked implementation.
* **Stan fits**: 4 chains × 600 warmup + 600 sampling, delta 0.9.
  **Zero divergences** in both realisations; split-R̂ ≈ 1.000 for ks,
  alpha_fert and sigma.
* **Coverage**: 10 of 11 parameters inside their 90% posterior
  intervals in each realisation, **20 of 22 overall** — the nominal
  rate for 90% intervals. Realisation 1 missed alpha_fert (posterior
  mean 0.221 against the true 0.35, with alpha_man compensating at
  0.367 against 0.25); realisation 2 missed beta_perm (−0.018 against
  −0.20). ks recovered at 0.124 [0.055, 0.218] and 0.184 [0.076,
  0.339] against the true 0.18.

### Identifiability warnings the test exposes

1. **Fertilizer vs manure split is weakly identified.** The two
   agricultural sources share one delivery factor and their spatial
   patterns are correlated, so the posterior trades alpha_fert
   against alpha_man; their delivered sum is much better identified
   than either coefficient. The real fit should report the combined
   agricultural delivery alongside the split.
2. **Delivery-modifier signs are not guaranteed at basin scale.**
   beta_perm and beta_tile had wide intervals covering zero in
   realisation 2. On the full network (181 stations, far more
   modifier contrast) this should improve, but the modifiers are the
   parameters to watch.
3. **Land-rate parameters are wide** (rate_urban 90% interval spans
   roughly a factor of 4–5; rate_bg similar). This is honest width
   from load data alone, not a defect; it argues for reporting
   source shares with their posterior intervals, not point values.

## The ks comparison (Alex's acceptance check)

Alex's earlier MultiBUGS-to-Stan port of a SPARROW model found the
in-stream attenuation coefficient ks "way narrow" in Stan relative
to BUGS. The comparison slot is therefore built into this build's
verification, with the posterior SD of ks as the reported quantity:

* **What was run here**: the same synthetic data were fitted with a
  single-block random-walk Metropolis sampler (log-transformed
  positive parameters, identical priors and likelihood, 12,000
  iterations, post-burn acceptance 0.20 and 0.24). This is a
  Gibbs-tradition reference sampler, **not a BUGS run — no BUGS
  run was performed in this build, and none is claimed**.
  Stan/Metropolis posterior-SD ratios for ks: **1.19** (realisation
  1: SD 0.0509 vs 0.0426) and **0.82** (realisation 2: SD 0.0813 vs
  0.0989). On these small synthetic cases the ratio brackets 1 — no
  systematic narrowing of the Stan posterior relative to the
  Metropolis reference is observed.
* **What will be run with the real fit** (the documented slot): fit
  the identical model and data (a) in Stan as here and (b) in a
  Gibbs sampler — MultiBUGS, per Alex's original comparison, if it
  can be run in that environment, otherwise the Metropolis reference
  in `run_synthetic_test.py` scaled to the full network — and report
  sd_ks(Stan) / sd_ks(Gibbs reference) beside the fitted ks. A ratio
  materially below 1 would reproduce Alex's finding and would be
  reported as such, with the delivery–decay ridge (alpha_fert /
  alpha_man against ks) examined as the mechanism.

## What a real-data fit still needs

* **A timed pilot first.** Per-evaluation cost scales with reach
  count: the full network is 125× the synthetic basin's 191 reaches.
  The synthetic Stan fit (4 × 1,200 iterations) took tens of seconds;
  a full-network fit should be budgeted at hours on this VM, not
  minutes, and a 200-iteration pilot on the real JSON should be timed
  before committing to chain lengths. (Per Alex's ruling the proposal
  package carries no fit; this build stops at a validated model, a
  validated JSON and a passed recovery gate.)
* The combined-agricultural-delivery reporting (warning 1) and the
  ks comparison (above) are part of the real fit's deliverables, not
  optional extras.
* Prior sensitivity on rate_dep: its informative prior (SD 0.25 in
  log space) is the one prior that can move a source share; a
  widened-prior variant (SD 0.5) should accompany the real fit.
* Data vintage caveat carried from the design: all source and
  land-use fields are the deposit's 2002 base year, while the loads
  are 1990–2021 means. The fit is a mean-annual model on those terms,
  as designed; nothing in this build changes that.

Nothing in this folder has been committed to git.
