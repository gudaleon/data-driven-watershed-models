# BRIEF — SPARROW AST / BMP management layer

**Full brief:** `~/workspace/your_files/Orchestration_Brief_Ensemble_R_Members.md`
(Toby's orchestration brief to Claude Code / Codex / Grok Build). This file is
the per-member section of that brief. The full brief's hard rules, package
format, and handoff template govern this member.

## The whole idea

George Arhonditsis proposed a regional Bayesian SPARROW model of phosphorus
export for the tributary watersheds draining to Lake Ontario between
Oakville and Cobourg. Alex adopted an ensemble design on 2026-11-20 (after
Neumann et al. 2021, Environmental Reviews): structurally different
data-driven models estimate the same phosphorus loads in parallel. The AST
with BMP analysis is the ensemble's **management layer**: it fits nothing of
its own. It takes the members' predicted loads by reach, applies named
best-management-practice source-reduction scenarios, and returns reduced
loads and source shares — at project stage it will run those scenarios
through the members' posterior draws so management results come out as
posterior predictive results. The proposal package (v16,
`~/workspace/fei_lake_ontario/Ontario_SPARROW_proposal_package_v16.docx`)
carries no model fit, by Alex's ruling; this build assembles the scenario
engine and proves it on the acceptance tests below.

## Method of record — assemble, do not invent

- **Starting code:** Alex's own repository `gudaleon/AST` (Python,
  "Decision tool on agriBMP", private).
- **BMP effectiveness method:** Roland et al. 2022, Journal of Soil and
  Water Conservation 77: 15–29.
- **Interface reference:** the RSPARROW Shiny Decision Support System
  scenario design. The classical RSPARROW core is already ported in
  `sparrow_bayesian/`; its PORT_NOTES.md records that scenario execution
  (`predictScenarios*`) was deliberately left out of the port because it
  belongs to this member, and that the scenario engine calls the same
  accumulation the port implements.
- "CAST" (Chesapeake Assessment Scenario Tool) is **not** a SPARROW AST
  despite the similar name — do not use it.

## Access condition

`gudaleon/AST` is private. Work from a clone Alex provides. If you cannot
access it, state that in your handoff report and build from this brief's
design alone, flagging every place where the AST repository's own
definitions could not be checked. Do not guess at its contents.

## Inputs

- A member model's predicted loads by reach, with source apportionment by
  reach, in the reach table's terms. For the machinery tests, a small
  synthetic reach table you construct is sufficient and must be labelled
  synthetic; the extended network and its attributes live at
  `~/workspace/fei_lake_ontario/sparrow/extended_reaches_v1.gpkg` and
  `reach_attributes_v1.csv` in the same folder for reference.
- BMP effectiveness ranges: Roland et al. 2022 (above), or the AST
  repository where it defines them.

## Work

A scenario engine that: reads predicted loads by reach with their source
shares; applies a named scenario (a set of BMPs, each acting on named
sources with an effectiveness range); and returns reduced loads and source
shares per scenario, by reach and aggregated to basin outlets.

## Acceptance tests (all evidenced in `tests/`)

- **Zero-reduction identity:** with zero reductions in every source, the
  engine reproduces its input loads exactly, per reach and in aggregate.
- **Single-source reduction:** a uniform reduction applied to one source
  changes only that source's share, and changes the total by exactly the
  stated amount.
- **Sourcing:** every effectiveness value in every scenario carries its
  source (Roland et al. 2022 table or page, or the AST repository location
  it came from). A scenario containing an unsourced effectiveness value
  fails.

## Rules that bind this member

Code only in this subfolder — no data files committed. Evidence or it is
not done. No headline calibration claims at build stage. Nothing is pushed,
published, sent or filed without Alex's explicit approval of the exact
material. Toby verifies the return independently; disagreements are
reconciled against source data, never averaged.
