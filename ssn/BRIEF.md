# BRIEF — SSN (spatial stream network models), statistical member

**Full brief:** `~/workspace/your_files/Orchestration_Brief_Ensemble_R_Members.md`
(Toby's orchestration brief to Claude Code / Codex / Grok Build). This file is
the per-member section of that brief. The full brief's hard rules, package
format, and handoff template govern this member.

## The whole idea

George Arhonditsis proposed a regional Bayesian SPARROW model of phosphorus
export for the tributary watersheds draining to Lake Ontario between
Oakville and Cobourg. Alex adopted an ensemble design on 2026-11-20 (after
Neumann et al. 2021, Environmental Reviews): structurally different
data-driven models estimate the same phosphorus loads in parallel, and their
agreement or disagreement is part of every result. SSN/SSNM is one of the
two statistical members (with E-EMMA). The proposal package (v16,
`~/workspace/fei_lake_ontario/Ontario_SPARROW_proposal_package_v16.docx`)
carries no model fit, by Alex's ruling; this build makes the member ready to
fit at project stage.

This member is built on a machine **with R** (state R and package versions
in the return): SSN2's maintained core is an R package and Toby's VM has no
R interpreter.

## Method of record — adopt, do not reimplement

SSN2, US EPA (Dumelle, Peterson, Ver Hoef et al.),
https://github.com/USEPA/SSN2, GPL-3, CRAN v0.4.0 — the Ver Hoef & Peterson
(2010) tail-up / tail-down / Euclidean covariance family. Network
preparation by the companion package SSNbler. The legacy SSN package
(github.com/jayverhoef/SSN) is dormant and is not used.

## Inputs (data stays in the project workspace; reference by path)

- Extended network: `~/workspace/fei_lake_ontario/sparrow/extended_reaches_v1.gpkg`
  (23,963 reaches; domain 1 = Oakville–Cobourg corridor, domain 2 = adjacent
  calibration ring), attributes in `reach_attributes_v1.csv` in the same
  folder. Provenance: mid-continental SPARROW deposit, Vouk et al. 2018,
  DOI 10.4224/300.0001; source-field base year 2002.
- Stations: `~/workspace/fei_lake_ontario/sparrow/stations_snapped_v4.csv`
  (470 stations; `station_id`, `reach_sparrowid`, `on_network`, tiers).
  Note: the file is in `sparrow/`, not `maps/`.
- Loads: `~/workspace/fei_lake_ontario/sparrow/station_year_loads_v4.csv`
  (3,575 station-years, 187 stations, 1990–2021; `load_t`, `se_t`,
  `method`, `on_network`, `reach_sparrowid`). 181 of the 187 load stations
  are on-network; 15 snapped stations are off-network (Central Lake Ontario
  shore strip, where the source network holds no reaches) — exclude them
  from network modelling with the exclusion stated, never dropped silently.

## Work

1. With SSNbler, build the spatial stream network dataset: extended-network
   reaches; observation sites from the snapped stations joined to mean
   annual loads from the loads file (mean over each station's record of
   annual TP load; response = mean annual log load; the exact aggregation
   is stated in NOTES.md and scripted, never hand-computed).
2. Fit the reference SSN model (tail-up, tail-down, Euclidean components,
   component selection per SSN2's documentation).
3. Return the fitted model summary, covariance parameters with standard
   errors, and predictions at the calibration stations.

## Acceptance tests (all evidenced in `tests/`)

- SSNbler network builds with **zero topology errors**; every rejected
  reach or site is listed with identifier and reason.
- The fit **converges**; every covariance parameter is reported with a
  standard error.
- Predictions CSV keyed by `station_id` (character, exactly as input), one
  row per on-network load station in the fit, with fitted value and
  prediction standard error.
- The member runs end to end from its numbered R scripts alone, with the
  final run's log (including `sessionInfo()`) in `logs/`.

## Rules that bind this member

Code only in this subfolder — no data files committed. Evidence or it is
not done. No headline calibration claims at build stage; fitted results are
labelled with the data they were fitted on. Nothing is pushed, published,
sent or filed without Alex's explicit approval of the exact material. Toby
verifies the return independently; disagreements are reconciled against
source data, never averaged.
