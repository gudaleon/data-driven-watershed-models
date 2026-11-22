# NEWS2 inputs and calibration data — extraction and comparison v1

**Purpose.** This document is organised as the comparative forcing table that the
ensemble's forcing-data comparison will be built from: one row per forcing
variable, NEWS2's value and source beside our project's value and source, with
the discrepancy stated. The livestock/manure chain is traced end to end first,
because the user has directed that every livestock number be treated as suspect
until compared, and because issues were previously found in livestock estimates
on our own side (the Vouk audit, §4).

**Companion files** (same folder):

- `news2_inputs_inventory.csv` — 33 input variables declared in `vars.cfg`
  [IN.VARS], each with its phosphorus-submodel role, input table, column, units
  and the source of the unit documentation. 17 of the 33 are consumed by the
  phosphorus sub-models (DIP, DOP, PP); the rest are N-form, DSi or context
  variables declared in the same configuration.
- `news2_calibration_constants.csv` — 41 calibration constants and code-embedded
  coefficients for the P sub-models, each with its source file.
- `news2_regional_inputs.csv` — the full basin-level input row (all six input
  tables plus the reference output row, 176 columns) for the selected regional
  basin.

**Units caveat, stated once for the whole document.** The GNE Manual does not
tabulate input units; it defers them to Mayorga et al. (2010), Supplementary
Tables A2/A3/A4, which are **not in the repository**. Units below are therefore
the code's working conventions (as implemented in `globalnews.py` and in our
verified port `news2.py`): diffuse and point-source inputs and all yields in
kg km-2 yr-1 over basin area, loads in Mg yr-1, runoff in mm yr-1, land-use
columns in percent of basin area, retention columns as fractions in [0, 1].
Where a unit could not be confirmed from the repository, the inventory CSV
says so in its `units_source` column.

---

## 1. Regional basin selection

STN-30 v6 contains **exactly one** basin that intersects the Lake Ontario
drainage and southern Ontario:

| BasinID | Name | Area (km2) | Land area (km2) | Mouth |
|---|---|---|---|---|
| 24 | SAINT LAWRENCE | 1,052,470 | 792,309 | -70.75, 47.25 (St. Lawrence estuary) |

**Selection rule:** the unique STN-30 v6 basin whose drainage contains Lake
Ontario and southern Ontario. The entire Great Lakes–St. Lawrence system is a
single lumped STN-30 basin; **no separate Lake Ontario tributary basins exist
in STN-30 v6**. Every NEWS2 number below is therefore a whole-system average
over 1.05 million km2 — the Great Lakes, the Canadian Shield headwaters and
the US tributaries included — and is never geographically aligned with our
corridor. All per-hectare-of-agricultural-land normalisations below exist
precisely to make the least-misaligned comparison that the data allow.

Our side's geographies: the Ontario SPARROW corridor domain (domain 1,
4,600 km2, Oakville–Cobourg shore) and the extended network (domain 2,
26,912 km2), both on the binational SPARROW deposit attributes (Vouk et al.
2018, `data1`; base year 2002) as extracted in
`~/workspace/fei_lake_ontario/sparrow/reach_attributes_v1.csv`; and, for the
livestock chain, the Census-of-Agriculture-based build for the Canadian Lake
of the Woods (LoW) basin (GrokBot nutrient-inputs package, locally held at
`~/workspace/lake_of_the_woods/invoice_2026_27/grokbot_np_inputs/`, with the
census machinery documented in
`~/workspace/lake_of_the_woods/proposals/census_agriculture_analysis_v1.md`).

---

## 2. The livestock/manure chain, traced end to end

### 2.1 NEWS2's chain — it ends at one basin yield

Inside the repository, the livestock chain has exactly three links:

1. **One input value per basin:** `c0TPma` (variable `WSdif_ma_P`), manure P
   input as a basin-area yield, kg P km-2 yr-1, in
   `inputs/preproctbls/c00_diffsrc_landuse_ng.csv`. Its N analogue is
   `c0TNma`. For basin 24: **TPma = 93.38, TNma = 464.74 kg km-2 yr-1**
   (manure N:P = 4.98).
2. **The explicit diffuse sum:** `WSdif_ant_P = WSdif_fe_P + WSdif_ma_P -
   WSdif_ex_P`, floored at 0 (fertilizer + manure - crop export). Basin 24:
   117.19 + 93.38 - 117.95 = **92.63 kg km-2 yr-1**; manure is 44% of the
   fertilizer + manure input.
3. **Two multipliers:** the watershed export fraction FEws (a function of
   natural runoff only; DIP FEws = 0.0668 for basin 24) and the river export
   fraction FEriv (DIP 0.355, DOP 0.988 for basin 24). The manure-attributed
   DIP yield in the reference output (`c0YS2difMa_ant_DIP`) is 0.973
   kg km-2 yr-1.

**Everything upstream of link 1 is absent from the repository.** There are no
animal populations, no per-head excretion coefficients, no manure handling
fractions (grazing / storage / application / losses) anywhere in the code,
the configuration files or the input tables. Those steps live in the
IMAGE/MNP preprocessing of Bouwman and co-workers that produced the c00
diffuse-source tables. **Consequence: NEWS2's livestock chain cannot be
audited below the basin-yield level from this deposit.** Any error in animal
numbers, excretion rates or handling fractions is invisible inside NEWS2; it
can only be caught by comparing the resulting basin yields against
independent builds — which is what §3 and §4 do.

### 2.2 Our chain — two builds, one of them already audited

- **Vouk et al. 2018 (the `data1` SPARROW deposit; base year 2002).** This is
  the source of the corridor attributes (`man_p_us_kg`, `fert_p`,
  `pcs_p02_s`). Canadian livestock from the 2001 Interpolated Census of
  Agriculture by Sub-Sub-Drainage Area (SSDA); 11 livestock classes (Table
  C-3 — **no calves and no bulls**); Ruddy et al. 2006 coefficients, with the
  **US delivery coefficients applied to all areas, Canada included** (report
  p. C-2); manure allocated to catchments by agricultural + pasture land
  fraction.
- **Census-based build (LoW nutrient-inputs package).** Census of Agriculture
  heads by animal class (2021 CD-level heads are held for all Ontario census
  divisions in `livestock_heads_by_cd_2021_AN_v1.csv`; CCS-level for LoW) ×
  per-head P coefficients (Falcone 2021 defaults, e.g. beef cows 19.35, dairy
  cows 13.46, other cattle 10.27, pigs 4.89 kg P head-1 yr-1) × 365-day herd
  correction factors from StatCan inventory series. Canadian LoW basin
  manure P: **328.5 t yr-1 (2016), 317.9 t yr-1 (2021)**; manure N:P = 3.44.

### 2.3 What the Vouk audit already found (our side's known livestock defect)

The audit at `grokbot_np_inputs/VOUK_AUDIT_AN_v4.md` rebuilt Vouk's method for
the Canadian LoW basin and compared it with the deposit's own input file:

- SPARROW deposit Canadian LoW manure P = **139.7 t yr-1 (2002)** against
  **328.5 t yr-1** from the census-based build (2016) — the deposit is low by
  **2.35x** on the direct comparison, and by **2.4–2.6x** against Vouk's own
  method made complete (audit decomposition: 139.7 → 246.5 t with the
  cross-border loss restored → 331–346 t with calves and bulls added).
- Cause 1 (43% of the rebuild): ICoA SSDA polygons 05PC/05PD cross the US
  border; manure was spread over the whole polygon's agricultural land
  including US land, and the share placed on US land was then dropped
  (376 t N / 107 t P lost). This is a **border-SSDA phenomenon**; it does not
  apply to the wholly-Canadian SSDAs of the Ontario corridor.
- Cause 2 (applies everywhere in Canada, corridor included): the calves and
  bulls classes are absent from Vouk's Table C-3, omitting a further
  85–99 t P in LoW alone.
- Ruled out by the audit's reproductions: a recoverable-manure factor, a unit
  error, double losses, a wrong year, and Canadian manure leaking into US
  catchments.

---

## 3. The comparative forcing table

Rows are forcing variables. "NEWS2" is basin 24 (SAINT LAWRENCE, year-2000
c00 inputs) throughout. "Ours" is the corridor domain 1 unless a row says
otherwise. Per-hectare rates are per hectare of agricultural land
(`c0agrpct` / `lulc_ag_perc` respectively) unless stated.

| Forcing variable | NEWS2 value | NEWS2 source | Our value | Our source | Discrepancy | Note |
|---|---|---|---|---|---|---|
| Basin area | 1,052,470 km2 | STN30v6ngNEWS.csv, BASINID 24 | 4,600 km2 (domain 1); 31,511 km2 (domains 1+2) | reach_attributes_v1.csv, demiarea sums | Not a discrepancy — different objects | NEWS2 lumps the whole Great Lakes–St. Lawrence system into one basin; no finer NEWS2 geography exists. All basin-yield comparisons below carry this caveat. |
| Agricultural land share | 19.81% | c00_diffsrc_landuse_ng.csv, c0agrpct | 21.6% (domain 1); 18.5% (domain 2) | data1 lulc_ag_perc (2002) | Shares agree within 3 points | The one well-aligned quantity: both models see roughly one fifth of the land as agricultural. This is what makes the per-hectare rows below meaningful. |
| Manure P input, basin yield | 93.38 kg km-2 yr-1 | c0TPma (IMAGE/MNP chain, not in repo) | 189.4 (domain 1); 353.9 (domain 2) kg km-2 yr-1 | data1 man_p_us_kg (Vouk 2002) | Ours 2.0x (domain 1) to 3.8x (domain 2) NEWS2 | Partly geography (corridor is the intensive shore strip; NEWS2 averages the Shield). Partly the livestock chains themselves — see per-hectare row. |
| **Manure P input, per ha agricultural land** | **4.71 kg ha-1 yr-1** | c0TPma / c0agrpct | **8.76 (domain 1); 19.14 (domain 2) kg ha-1 yr-1** | data1 man_p_us_kg / lulc_ag_perc | **Ours 1.9x (domain 1) to 4.1x (domain 2) NEWS2** | The headline livestock discrepancy. Domain 2 (Grand River) is Ontario's most intensive livestock country, so part is real; but see the next row — our data1 value is itself biased LOW for Canada, so the true gap against NEWS2 is larger, not smaller. |
| Manure P, census-based cross-check (LoW, per ha farm area) | 4.71 kg ha-1 yr-1 (per ha ag land, basin 24) | as above | 4.93 kg ha-1 yr-1 (2016, per ha total farm area, Canadian LoW basin) | Census of Agriculture heads x Falcone 2021 coefficients, 365-day herds (grokbot_np_inputs v2) | Rates look close — but denominators differ (ag land vs total farm area incl. pasture/woodland) | On a like denominator (crop + pasture land) the LoW census rate is ~7.1 kg ha-1 yr-1, i.e. ~1.5x NEWS2. And against data1 for the SAME LoW basin the census build is 2.35x higher (§2.3): the two "independent" global/regional products disagree with the census build in opposite directions. |
| Manure N:P ratio | 4.98 | c0TNma / c0TPma | 3.44 | LoW census build (manure N 1,128.5 t / P 328.5 t, 2016) | NEWS2 manure is 45% more N-rich | Excretion N:P is set by animal mix and coefficients. A 45% gap means at least one chain's class mix or coefficients are materially off; NEWS2's chain cannot be inspected below basin level (§2.1), ours can (Falcone class coefficients, held). |
| Fertilizer P input, per ha agricultural land | 5.91 kg ha-1 yr-1 | c0TPfe / c0agrpct | 16.35 (domain 1); 15.46 (domain 2) kg ha-1 yr-1 | data1 fert_p / lulc_ag_perc | Ours 2.6–2.8x NEWS2 | Same geography caveat, but a 2.7x gap on the corridor strip is larger than the land-share difference explains. NEWS2 fertilizer is year-2000 IMAGE/MNP; data1 is 2002. Census cannot arbitrate masses (it records fertilized area, not kg) — see §5. |
| Crop P export (removal), per ha agricultural land | 5.92 kg ha-1 yr-1 | c0TPex / c0agrpct | not held | — | Cannot be quantified | We hold no crop P uptake/removal field. NEWS2's export (5.92) slightly exceeds its own fertilizer + manure input per ha would suggest is sustainable alongside its manure term; the explicit diffuse remainder (fe + ma - ex) is 92.6 kg km-2 yr-1 basin-wide. Reference needed (§5). |
| Sewage P point source (excreta + detergent), basin yield | 21.79 kg km-2 yr-1, passed to the river with **zero removal** (FEpnt_DIP = 1) | c0TPsewhum + c0TPsewdet; constants.cfg [CAL.DIP] c = 1 | 39.8 (domain 1); 16.7 (domain 2) kg km-2 yr-1 | data1 pcs_p02_s | Not directly comparable | Two different objects: NEWS2's is a generated sewage load before any treatment removal — and the basin's own removal column `c0TPfrem` = 57.7% sits in the input table **wired to no P variable** (the N forms do apply their `c0TNfrem`). data1's pcs is a delivered point-source load. The unused c0TPfrem is a finding in its own right: either TPsewhum is already net of removal (undocumented in the repo) or NEWS2 overstates DIP point sources by up to ~2.4x on this term. |
| Natural runoff | 464.9 mm yr-1 | c0R (basin 24) | 289 (domain 1); 398 (domain 2) mm yr-1 | data1 runoff_mm, area-weighted | NEWS2 1.2–1.6x ours | Geography again: basin 24's runoff includes the lake surfaces and the northern tributaries. Not an error signal. |
| Annual precipitation | 854.6 mm yr-1 | c0annprecp (basin 24) | 823 (domain 1); 879 (domain 2) mm yr-1 | data1 precip_mm, area-weighted | Within 4% | The climate forcing is aligned; the divergences above are therefore in the nutrient chains, not the weather. |
| Reservoir retention, DIP (D_DIP) | 0.641 (code caps the input at 0.85) | c0Ddip0to1; Calc_FEriv cap | Reach-scale retention from hydraulic load (hload) — no single basin number | data1 hload / rhload fields | Structural, not comparable as scalars | NEWS2 applies one lumped retention to the whole St. Lawrence system, Great Lakes included; ours is per-reach. The ensemble comparison must be made on delivered loads, not on this parameter. |
| Sediment retention (D_TSS) | 0.718 | c0Dsed0to1 | — | — | — | Drives the PP sub-model; our SPARROW build has no TSS retention analogue (PP is source-driven). Noted for the GREEN/SSN members, which do carry sediment terms. |
| Discharge removal fraction (FQrem) | 0.01173 | c0FQremSclGEMS | not held | — | Cannot be quantified | Consumptive-use removal is not a forcing in our build. Small for this basin in any case. |
| Watershed export fraction, DIP (FEws) | 0.0668 at Rnat = 0.465 m yr-1 | computed (CAL.DIP a=0.85, b=2, e=0.29) | — (SPARROW estimates delivery coefficients instead) | — | — | Structural difference: NEWS2 fixes the delivery fraction globally from runoff alone; SPARROW fits it. Listed so the forcing table records that this "forcing" is a calibration choice, not data. |

---

## 4. The largest livestock-chain discrepancy

**In one paragraph:** for the one basin where all three chains can be placed
on the same ground — the Canadian Lake of the Woods basin — the SPARROW
deposit's manure P (139.7 t yr-1, Vouk 2002) is 2.35x below the
Census-of-Agriculture build (328.5 t yr-1, 2016), and the Vouk audit
attributes the gap to two named defects in the deposit's livestock chain: a
cross-border SSDA allocation that dropped 107 t of P (43% of the rebuild)
onto US land, and the complete absence of calves and bulls from its livestock
class list (a further 85–99 t P). NEWS2's chain cannot be given the same
treatment because it holds no animal-level data at all — but its surviving
basin yield implies a manure P rate (4.71 kg ha-1 of agricultural land) that
sits ~1.5x below the census build on a like denominator and carries a manure
N:P ratio 45% above it, while our corridor attributes (themselves Vouk
data, biased low for Canada by the missing classes) already run 1.9–4.1x
above NEWS2 per hectare of agricultural land. The three chains therefore
disagree pairwise by factors of roughly 1.5 to 4 on the same forcing, and
the only chain that can be taken apart animal by animal — the census
build — is the one that says both global/regional products are too low.

## 5. Reference data still needed

1. **Mayorga et al. 2010 Supplementary Tables A2/A3/A4** — the authoritative
   units and input-provenance documentation for every NEWS2 variable; the
   GNE Manual defers to them and they are not in the repository. Until they
   are in hand, all NEWS2 units rest on code convention.
2. **The IMAGE/MNP manure preprocessing behind c0TPma/c0TNma** (Bouwman and
   co-workers): animal numbers by system, excretion coefficients, and the
   manure handling fractions (grazing, storage, application, non-soil uses).
   Without it, NEWS2's livestock chain stops at the basin yield (§2.1).
3. **A corridor-scale census livestock build.** The ingredients are held —
   2021 Census of Agriculture heads by census division for all Ontario CDs
   (`livestock_heads_by_cd_2021_AN_v1.csv`) and the Falcone class
   coefficients used in the LoW build — but the corridor CD aggregation
   (the LoW v2 table is CCS-level, LoW units only) has not been built. This
   is the single most valuable missing reference: it would put a
   census-based manure P rate on the exact corridor geography.
4. **data1 manure N at reach level** (`man_n_us_kg` exists in the deposit's
   data1 file; `reach_attributes_v1.csv` carries the P field only) — needed
   for a corridor manure N:P to set against NEWS2's 4.98 and the census
   build's 3.44.
5. **Crop P removal coefficients** for the corridor's crop mix (to build our
   analogue of `c0TPex`); the census gives areas, never nutrient masses.
6. **Documentation of whether c0TPsewhum is gross or net of treatment
   removal** — i.e. whether the unused `c0TPfrem` column is redundant or a
   disconnected input (§3, sewage row). This is a question for the NEWS2
   supplement (item 1) or the Global NEWS group, not something the code
   answers.
7. **The Vouk audit extended to the corridor's SSDAs.** The border-drop
   defect does not apply (corridor SSDAs are wholly Canadian), but the
   calves/bulls omission does, and its corridor magnitude has not been
   computed; the LoW decomposition (+71% on the kept-whole rebuild from
   that one cause) is the only available guide.

## 6. What was extracted (counts)

- Input variables inventoried: **33** (17 consumed by the P sub-models),
  in `news2_inputs_inventory.csv`.
- Calibration constants extracted: **41**, in
  `news2_calibration_constants.csv` — the [CAL.DIP] and [CAL.DOP] sections
  of `constants.cfg`, plus the code-embedded TSS-regression coefficients,
  the eight Amiotte-Suchet lithology coefficients, four truncation ranges,
  the Ludwig %POC polynomial, the PP equation coefficients, the runoff
  floor and the retention/TSS caps from `globalnews.py`.
- Regional basins selected: **1** (BASINID 24, SAINT LAWRENCE), full input
  and reference-output row in `news2_regional_inputs.csv`.
- Nothing in the original code or in the port was modified, and nothing was
  committed to git.
