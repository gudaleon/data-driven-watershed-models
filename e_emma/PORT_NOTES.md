# PORT NOTES — E-EMMA Python port

## Provenance and licence

- **Method**: Kirchner, J.W. (2023), *Mixing models with multiple,
  overlapping, or incomplete end-members, quantified using time series
  of a single tracer*, Geophysical Research Letters 50,
  e2023GL104147.
- **Original code**: Kirchner, J.W. (2023), *eemma.R, an R script for
  Ensemble End-Member Mixing Analysis*, EnviDat,
  doi:10.16904/envidat.410. Version 1.0, build 2023.03.27.
  © 2023 ETH Zurich and James Kirchner. The scripts are GPL-3; the
  user's guide is CC-BY-SA.
- **This port** (`eemma.py`) is a derivative work of `eemma.R` and is
  likewise **GPL-3**. The author's original files are kept alongside
  for reference: `original_eemma.R`, `original_vignette_script.R`,
  `vignette_data.txt` (all © Kirchner / ETH Zurich, GPL-3 / as above).
- The placeholder download at
  `~/workspace/fei_lake_ontario/codebase/eemma_kirchner.zip` was an
  EnviDat HTML page, not the archive; the real resource files were
  fetched from the dataset's CKAN download URLs (package
  `48dff6cd-bfcb-4c5e-a30d-f909dc8d36d3`).

## What the method computes

Classical end-member mixing analysis needs several tracers and a
complete, non-overlapping end-member set. EEMMA instead exploits the
**time-series correlations of a single tracer**: the mixture series is
regressed on the end-member series, and the regression coefficients
*are* the mixing fractions. Around that core, the script provides:

- **Complete sets** (`complete=True`): fractions are constrained to
  sum to 1, imposed exactly by a difference transformation (subtract
  one end-member from the mixture and all others), regression in the
  transformed space, and an analytic back-transformation of the
  coefficients and their covariance matrix.
- **Incomplete sets** (`complete=False`): an intercept is fitted; the
  unmeasured sources collapse into one unknown end-member whose
  fraction is `1 − sum(f)`, and whose **tracer time series and its
  standard error** are reconstructed from the residuals by
  first-order second-moment (delta-method) error propagation. This is
  the diagnostic that detects a missing end-member.
- **Memory** (`memory=True`): the lagged mixture is added as a
  regressor; its rescaled coefficient is a first-order carry-over
  ("lagmix") fraction, with fractions and covariances propagated
  through the rescaling Jacobian.
- **Two estimators**: OLS, and an errors-in-variables (EIV) estimator
  implementing Fuller's (1987, *Measurement Error Models*, eqs.
  2.5.13/2.5.14) bias correction for measurement/sampling error in the
  end-member series. The correction factor λ is the smallest root of
  `det(mZZ − λ·Szz)`, found by a grid scan plus Brent refinement, with
  an asymptotic-root fallback; `uncert="relative"` (default) rescales
  the supplied error ratios inside the fit, `"absolute"` takes them
  literally; `alpha` (default 1) trades bias against variance.
- **Uncertainty**: analytic throughout — coefficient covariance from
  the (corrected) cross-product matrices, delta-method propagation
  through every back-transformation, one-tailed t p-values for the
  null `f ≤ 0`, as in the original.

The port replicates the R code path-for-path, including its input
validation, complete-case filtering, constant-end-member guards, the
R `matrix()` column-fill semantics in the memory Jacobian, and the
one-tailed p-value convention. R's `uniroot` is `scipy.optimize.brentq`;
`pt` is `scipy.stats.t.sf`.

## Verification

**Blocker, stated plainly:** R is not installed on this machine (no
`Rscript` anywhere on the filesystem; only the user package library
`~/R/library` survives), so the original script could not be executed
here for a direct run-to-run comparison. Per standing instructions R
was not reinstalled.

Verification therefore used the strongest available substitute: the
**R script's own printed outputs published in Kirchner's user's
guide**, which reproduces all eleven vignette analyses of
`eemma_vignette_script.r` on the demo data, plus the printed `f_cov`
matrix for the stream case. `verify_against_guide.py` compares the
port against every one of those numbers:

| Check | n values | Worst discrepancy |
|---|---|---|
| Mixing fractions `f`, 11 vignette cases | 40 | 4.7 × 10⁻⁸ |
| Standard errors `se`, 11 vignette cases | 40 | 4.9 × 10⁻⁹ |
| `f_cov`, stream complete case | 9 | 4.1 × 10⁻¹⁰ |
| Independent KKT constrained-LS cross-check (stream, complete) | 3 | 4.2 × 10⁻¹⁵ |
| Independent `lstsq` cross-check (xylem1, incomplete) | 4 | 1.0 × 10⁻¹² |

The guide prints 6–9 significant digits, so agreement at ~10⁻⁸ is
agreement at print precision — the port reproduces the R outputs on
the demo data. The two cross-checks are independent of the port's
transformation algebra (a KKT solve for the sum-to-one constraint,
and a plain least-squares fit for the incomplete case), so the
transformation/back-transformation machinery is verified by a second
route as well. Headline case, stream ~ precip+soil+gw (complete):
port `f` = (0.2068911, 0.2101869, 0.5829221),
`se` = (0.01141577, 0.04916945, 0.05032645) — identical to the guide
at every printed digit (true fractions 0.20/0.20/0.60).

**EIV caveat:** the guide publishes no EIV numbers, so the EIV path
is verified only as a smoke test. With default settings (equal
relative errors, α = 1) it runs cleanly and behaves as the method
requires — fractions close to OLS with slightly larger standard
errors, λ small and positive: stream complete gives
`f` = (0.20632, 0.19879, 0.59489), λ = 0.008921; xylem1 incomplete
(precip+soil) gives `f_unk` = 0.43095, λ = 0.008519. A numeric
EIV verification against the original remains owed if R becomes
available; run `original_eemma.R` on the vignette cases with
`method="EIV"` and compare against `eemma.py` with `method="EIV"`.

## A Bayesian EIV version in Stan — how it would differ

Kirchner's EIV is a *method-of-moments* correction: it subtracts a
λ-scaled error-covariance estimate from the observed cross-product
matrices and then proceeds as if the corrected matrices were data.
Uncertainty comes from the usual OLS-style covariance of the
corrected estimator, pushed through delta-method Jacobians. Three
consequences follow: fractions are unconstrained (they can leave
[0, 1] and, in the complete case, are forced to sum to 1 only by the
transformation, not by a prior); the error variances are inputs, not
inferred quantities (with `uncert="relative"` only their ratios
matter, and their common scale is absorbed into λ); and every
uncertainty statement is a Gaussian/delta approximation around a
point estimate.

A Bayesian errors-in-variables model in Stan would instead treat the
**true end-member series as latent variables**: each observed series
is modelled as its latent truth plus measurement error with an
explicit (possibly inferred) error variance, and the mixture is the
fraction-weighted sum of the latent truths plus its own error. That
changes the analysis in specific ways:

- **Fractions get a simplex prior** (e.g. Dirichlet) in the complete
  case, so posterior fractions lie in [0, 1] and sum to 1 by
  construction rather than by algebraic transformation; the unknown
  end-member of the incomplete case becomes a latent time series with
  its own prior/process model instead of a residual reconstruction.
- **Error variances become parameters** with priors, estimated jointly
  with the fractions (identifiability then rests on the prior and on
  the correlation structure across time, which is exactly the
  information EEMMA exploits); λ and Fuller's α disappear — there is
  no bias-correction factor because there is no biased estimator to
  correct.
- **Uncertainty is the full posterior**, not a delta approximation:
  credible intervals for fractions, for the unknown end-member's
  series (a posterior band at each time step, replacing the
  first-order `unk_std_err`), and for derived quantities come from
  the same draws, with no Jacobian bookkeeping.
- **Memory** becomes a state-space/AR term on the latent mixture
  rather than an extra regressor, so its uncertainty propagates
  through the dynamics instead of through a rescaling Jacobian.
- Costs: the latent-series formulation has O(n·k) latent parameters,
  needs priors where the moment method needs none, and replaces a
  closed-form solve with MCMC — sampling must be checked (divergences,
  effective sample size) where the script's answer is exact arithmetic.
  The two should agree closely when errors are small relative to the
  signal variance, which is the regime where Kirchner's benchmarks
  show EIV and OLS nearly coincide.

## Files

- `eemma.py` — the port (`eemma`, `eiv`, `ols`).
- `run_vignette.py` — the eleven vignette analyses on the demo data.
- `verify_against_guide.py` — the verification gate above (exit 0 =
  pass; run with `~/workspace/mapvenv/bin/python`).
- `vignette_data.txt`, `original_eemma.R`,
  `original_vignette_script.R` — the author's demo data and scripts.
