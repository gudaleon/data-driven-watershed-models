# BRIEF — Whelan TOPMODEL-P, field-scale tier

**Full brief:** `~/workspace/your_files/Orchestration_Brief_Ensemble_R_Members.md`
(Toby's orchestration brief to Claude Code / Codex / Grok Build). This file is
the per-member section of that brief. The full brief's hard rules, package
format, and handoff template govern this member.

## The whole idea

George Arhonditsis proposed a regional Bayesian SPARROW model of phosphorus
export for the tributary watersheds draining to Lake Ontario between
Oakville and Cobourg. Alex adopted an ensemble design on 2026-11-20 (after
Neumann et al. 2021, Environmental Reviews): structurally different
data-driven models estimate the same phosphorus loads in parallel. The
Whelan TOPMODEL-P model is the ensemble's **field-scale tier**, to be run on
two or three agricultural catchments at project stage: it is the only member
that represents event-scale, topography-driven phosphorus export, and it is
a stochastic simulator — its outputs enter the ensemble as scenario
distributions, not as a likelihood fit. The proposal package (v16,
`~/workspace/fei_lake_ontario/Ontario_SPARROW_proposal_package_v16.docx`)
carries no model fit, by Alex's ruling; this build makes the member ready
for the project stage.

This member is built on a machine **with R** (state R and package versions
in the return): the maintained TOPMODEL core is an R package and Toby's VM
has no R interpreter.

## Method of record

- **Hydrology core (adopt):** dynatop,
  https://github.com/waternumbers/dynatop, GPL-2, CRAN — the maintained
  Dynamic TOPMODEL implementation.
- **Phosphorus layer (build from the paper):** Whelan, Hope & Fox 2002,
  Water Science and Technology 45(9): 167–176 — TOPMODEL saturation-excess
  variable source areas; dissolved P transferred vertically and laterally;
  particulate P through overland-flow sediment transport capacity; Monte
  Carlo over the soil surface P distribution; daily timestep. The paper is
  held at `~/workspace/fei_lake_ontario/literature/Whelan_Stochastic_phosphorus_model_TOPO.pdf`.
  No public code exists for the phosphorus layer (verified 2026-11-21);
  NOTES.md must say it is an equations build and trace every equation to
  the paper.

## Inputs

- The Whelan et al. 2002 paper (path above) — the phosphorus layer's only
  specification.
- dynatop's own example catchment — the hydrology verification case.
- Soil surface P distribution: **for this machinery build, use the
  distribution from Whelan et al. 2002 and flag the substitution** in
  NOTES.md and in every dependent output. The Ontario application will
  source this distribution from the OMAFRA provincial soil-test archive at
  project stage; that archive is not an input to this build.

## Work

1. Stand up dynatop on its own example catchment; record its outputs.
2. Implement the phosphorus layer on the dynatop hydrology, daily
   timestep, soil surface P distribution as a named input object.
3. Build the Monte Carlo driver: a stated number of realisations, a fixed
   seed recorded in the scripts, sampling the soil surface P distribution,
   returning the distribution of event and annual P export.

## Acceptance tests (all evidenced in `tests/`)

- With the phosphorus layer off, the hydrology **reproduces dynatop's own
  example catchment outputs**; any difference is explained in NOTES.md.
- The phosphorus layer's **mass balance closes per event** (inputs,
  exports, change in stored P) to a stated numerical tolerance, reported
  per event.
- The Monte Carlo driver runs the stated number of realisations under the
  fixed seed and returns the **distribution** of event and annual P export
  (quantiles at minimum), not a single value; the same seed reproduces the
  same distribution exactly.
- The soil surface P substitution (paper's distribution in place of the
  OMAFRA archive) is flagged wherever it applies.

## Rules that bind this member

Code only in this subfolder — no data files committed. Evidence or it is
not done. No headline calibration claims at build stage. Nothing is pushed,
published, sent or filed without Alex's explicit approval of the exact
material. Toby verifies the return independently; disagreements are
reconciled against source data, never averaged.
