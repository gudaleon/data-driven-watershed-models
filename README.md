# data-driven-watershed-models

An ensemble of structurally different phosphorus-export models, built to
run in parallel on the same catchments against the same calibration
record: the annual total phosphorus loads of the tributary watersheds
draining to Lake Ontario between Oakville and Cobourg, with an adjacent
calibration ring. The design follows Neumann et al. (2021),
Environmental Reviews: no single model with a single parameter
specification should carry a management decision, so several adequate
descriptions of different conceptual bases estimate the same loads, and
their agreement or disagreement is reported as part of every result.

## Members (one subfolder each)

| Folder | Model | Family |
|---|---|---|
| sparrow_bayesian | Bayesian SPARROW (backbone) | Hybrid semi-empirical, reach network |
| green_rgrid | GREEN-Rgrid (Malago et al. 2019) | Grid geospatial regression |
| element | ELEMeNT (Van Meter & Basu) | Legacy store + travel time |
| whelan_topmodel_p | Whelan TOPMODEL-P | Stochastic, topography-driven, field scale |
| ssn | Spatial stream network models | Statistical, on the network |
| e_emma | Ensemble end-member mixing analysis | Statistical, tracer-based |
| sparrow_ast_bmp | SPARROW AST with BMP scenarios | Scenario / decision layer |
| news2_screening | Global NEWS (NEWS2) | Screening annex, bounding estimates |

## Rules of this repository

- Code only. Data stays in the project workspace
  (~/workspace/fei_lake_ontario/); every brief names the exact input
  files and paths. Nothing in this repo is a data copy.
- Each subfolder holds a BRIEF.md: the whole idea, the inputs, the
  method of record and the acceptance tests for that member. A builder
  with no other context must be able to start from the brief alone.
- Evidence or it is not done: every build returns its test output.
- No headline calibration claims at build stage. Smoke tests run on
  synthetic data or named subsets; fitted results are labelled with
  the data they were fitted on.
- Nothing is published, sent or filed from this work without Alex
  Neumann's explicit approval of the exact material.
