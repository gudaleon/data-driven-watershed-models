# NEWS2 phosphorus port — notes

## What was ported

`news2.py` is a modern Python 3 port (numpy 2.5 / pandas 3.1, run with
`~/workspace/mapvenv/bin/python`) of the phosphorus forms of the Global
NEWS 2 model, from Mayorga's own GNE implementation in
`~/workspace/fei_lake_ontario/codebase/globalnews2gne` (GPL-2):

- **DIP** and **DOP** — dissolved sub-model: `globalnews.py`
  `model_dissolved()` with `Calc_PointSources()`, `Calc_WSExportFunc()`,
  `Calc_WSdifExplicit()`, `Calc_RSdifExportCoeff()`, `Calc_FEriv()`.
- **PP** — particulate sub-model: `globalnews.py` `model_particulate()`
  with `Calc_Yld_TSS_pred()` (Beusen et al. 2005 TSS regression, Ludwig
  et al. %POC relation, PP eq. 4).
- **TP** — DIP + DOP + PP, the same sum `allformmodel.py` uses for ICEP.

The NEWS2 GNE Manual (`NEWS2_GNEManual.pdf`, May 2009) was read for the
configuration/variable conventions; the equations and constants were
taken from the code and `constants.cfg`, which the manual designates as
the calibration source (Mayorga et al. 2010, Suppl. Table A2).

## What changed relative to the source

- The GNE plumbing is gone: no `vars.cfg`/`constants.cfg`/`gensetup.cfg`
  parsing, no IN/OUT global dictionaries, no `dbf.py`, no config-driven
  generic variable expansion. Input tables are read with plain pandas
  CSV readers, aligned on BasinID in the basin order of the basins
  table; the variable-to-column mapping of `vars.cfg` [IN.VARS] is kept
  as an explicit dictionary (`INPUT_COLUMNS`) in the module.
- Calibrated constants from `constants.cfg` ([CAL.DIP]: a=0.85, b=2,
  c=1, e=0.29, EC=26; [CAL.DOP]: a=0.95, c=0.010, e=0.010, EC=15) are
  module constants, cited to their source sections in comments, as are
  `R_MIN = 0.003` m/yr, the TSS-yield cap of 5000 ton km-2 yr-1, the
  lithology coefficients, the input truncation ranges and the D_DIP
  reservoir-retention cap of 0.85 — all from `globalnews.py`.
- Output is one tidy CSV: BasinID, basin area, and DIP/DOP/PP/TP yields
  (kg km-2 yr-1) and loads (Mg yr-1), plus FEws/FEriv diagnostics and
  the TSS yield. Source-contribution apportionment (SourceContrib /
  YS1 / YS2 variables) is not ported; it does not affect export values.
- No equation, constant, truncation range, cap or unit conversion was
  altered. Edge-case behaviour of the source (arid-basin zeroing at
  Rnat < 3 mm/yr, NoFlux zeroing in the PP chain, WSdif_ant floored
  at 0) is reproduced deliberately.

## Verification

Gate: the reference output bundled with the source repository,
`output/c00_NEWS2Output.csv`, produced by the original GNE code from
the bundled inputs in `inputs/preproctbls` (STN-30 basins). Running
`verify_port.py` compares the port against it basin by basin
(6,081 basins):

| Variable | Max abs. diff. | Max rel. diff. | Sum, port | Sum, reference |
|---|---|---|---|---|
| DIP yield | 5.0e-09 | 5.0e-12 | 139,414.72 | 139,414.72 |
| DIP load | 1.5e-07 | 4.8e-12 | 1,468,747.27 | 1,468,747.27 |
| DOP yield | 5.0e-11 | 5.0e-12 | 41,739.82 | 41,739.82 |
| DOP load | 4.9e-08 | 4.8e-12 | 598,972.87 | 598,972.87 |
| PP yield | 4.9e-09 | 4.8e-12 | 521,514.74 | 521,514.74 |
| PP load | 4.4e-07 | 4.9e-12 | 6,974,324.23 | 6,974,324.23 |
| TP yield | 5.2e-09 | 4.3e-12 | 702,669.28 | 702,669.28 |
| TP load | 4.3e-07 | 4.6e-12 | 9,042,044.37 | 9,042,044.37 |

Yields in kg km-2 yr-1, loads in Mg yr-1; sums are sums over basins of
the per-basin values. Every basin agrees within a relative difference
of 1e-9; the residual differences are at the printed precision of the
reference CSV (~11 significant digits), i.e. the port reproduces the
original model exactly to the precision the reference preserves.
Because a bundled reference output exists, running the original
(Python-2-era) code was not needed for the gate.

## Known limits

- Phosphorus only. DIN, DON, DOC, DSi, PN, POC and ICEP are not
  ported; the dissolved/particulate code paths are structured so the
  other forms can be added from the same source functions.
- Inputs must be NEWS2-format basin-aggregated tables (the six files
  named in `INPUT_COLUMNS`). Aggregating raw GIS data into those
  tables (the GNE GIS pre-processor, `gngis2tbls.py`) is out of scope,
  as in the parent tasking the port is a screening member: for the
  Lake Ontario corridor, basin inputs still have to be aggregated to
  the model's table format before this module can run on them.
- The constants are the year-2000 ("c00") global calibration of
  Mayorga et al. 2010. Nothing in the port re-calibrates them; at
  corridor scale the results are screening-level bounds, and the DIP
  reservoir retention and FQrem inputs must come from the c00-style
  pre-processed columns named in `vars.cfg`.
- Basins whose lithology class is missing from the source coefficient
  map (classes 1–8 only) will raise a KeyError rather than silently
  defaulting — deliberate, matching the source's closed class list.
