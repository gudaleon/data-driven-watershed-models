# NEWS2 phosphorus code audit v1 — correctness, not fidelity

**Scope.** This audit asks whether the NEWS2 phosphorus code is *right*,
not whether our port matches it. The port (`news2.py`) already reproduces
the original's bundled reference output to ~5×10⁻¹² relative difference
(see `PORT_NOTES.md`); that proves fidelity only. Here the published
model description is treated as the specification and both
implementations are graded against it.

**Sources of record.**
- Specification: Mayorga et al. (2010), *Environ. Modell. Softw.* 25,
  837–853 (equations cited by number below); Beusen et al. (2005),
  *Global Biogeochem. Cycles* 19, GB4S05 for the particulate sub-model;
  Ludwig et al. (1996) for the %POC relation. Supplementary Tables
  A2/A3 of Mayorga et al. are not in the journal PDF and could not be
  independently re-read; where a constant lives only there, this is
  stated.
- Original code: `globalnews2gne/globalnews.py` (`model_dissolved()`,
  `Calc_PointSources()`, `Calc_WSExportFunc()`, `Calc_WSdifExplicit()`,
  `Calc_RSdifExportCoeff()`, `Calc_FEriv()`, `model_particulate()`,
  `Calc_Yld_TSS_pred()`), `globalnews2gne/constants.cfg`,
  `globalnews2gne/vars.cfg`, `globalnews2gne/allformmodel.py`.
- Port: `news2.py` (`_dissolved()`, `_tss_yield_pred()`,
  `_particulate_pp()`).
- The NEWS2 GNE Manual (May 2009) documents the framework and
  designates `constants.cfg` as the calibration source of record; it
  contains no equations.
- Provenance caveat from the paper itself (Software availability):
  *"An earlier version of the model code was used to produce the
  results described here."* The released code is therefore not
  guaranteed to be the exact code behind the paper's figures.

## Equation-by-equation comparison

| Paper | Equation | Original code | Port | Verdict |
|---|---|---|---|---|
| Eq. 1 / 19 | Yld_F = FEriv_F · (RSpnt_F + RSdif_F); Ld_F = Yld_F · A (Mg/yr) | `model_dissolved()` identical; Ld = Yld·A/1000 converts kg→Mg | identical | ✓ all three agree |
| Eq. 2 / 3 | RSpnt_E from Van Drecht effluent accounting; RSpnt_F = FEpnt_F · RSpnt_E, FEpnt_F = c_F for non-DIN forms | `Calc_PointSources()`: FEpnt = `constants.cfg` c (DIP 1.0, DOP 0.010); RSpntExc/RSpntDet columns per `vars.cfg` (c0TPsewhum, c0TPsewdet) — the (1−hwfrem)·I factors of Eq. 2 are pre-computed into the input columns | identical mapping, verified against `vars.cfg` | ✓ (Eq. 2 lives in pre-processing, by design) |
| Eq. 4 / 7 | WSdif_ant,P = WSdif_fe + WSdif_ma − WSdif_ex, floored at 0; WSdif_nat,P = 0 | `Calc_WSdifExplicit()` identical, floor applied | identical | ✓ |
| Eq. 8 / 9 | RSdif_expl,F = FEws_F · WSdif_E; FEws_F = e_F · f_F(Rnat) | identical | identical | ✓ |
| Eq. 14 | f_F(Rnat) = Rnat^a_F (all forms except DIP) | `Calc_WSExportFunc()`: DOP (and DIN/DON/DOC) branch | identical | ✓ |
| Eq. 15 | f_DIP(Rnat) = 1 / (1 + (Rnat/a)^(−b)) | identical form, plus a near-zero floor (Finding F4) | identical | ✓ form; see F4 |
| Eq. 10 | RSdif_ec,F = f_F(Rnat) · EC_F (note: f, **not** FEws — no e multiplier) | `Calc_RSdifExportCoeff()` identical | identical | ✓ |
| Eq. 11 / 12 | RSdif split: explicit via FEws, EC term split Agfr / (1−Agfr) | identical | identical | ✓ (but see F6 on Agfr > 1) |
| §3.1.4 | f_F = 0 when Rnat < 0.003 m/yr (arid basins); point sources may still export | `f_F[Rnat < R_MIN] = 0`, R_MIN = 0.003; RSpnt untouched | identical | ✓ documented choice (see F8) |
| Eq. 17 / 18 | FEriv_F = (1−L_F)(1−D_F)(1−FQrem); L_F = 0 for P forms; D_DOP = 0 | `Calc_FEriv()` identical **plus an undocumented cap D_DIP ≤ 0.85** | cap reproduced | ✗ code–paper mismatch (F2); fraction guards DIN-only (F3) |
| §3.2.2 / Disc. A2 | D_F = discharge-weighted reservoir retention (pre-computed input) | input column c0Ddip0to1, then capped | identical | see F2 |
| Eq. 25 | Yld_TSS = (1 − D_TSS) · Yld_TSS,pred | `model_particulate()` identical | identical | ✓ |
| Beusen 2005 Table 2 / Disc. A4 | ln Yld_TSS,pred = 0.5522 + 0.0301·MargGrass + 0.1234·WetlnRice + 0.2933·FrnrPrecp + 0.1027·FrnrSlope + LithoChem; inputs truncated to the development-set ranges (paper Table A3); cap 5000 Mg km⁻² yr⁻¹; class 8 (ice/no overlay) → 0 | `Calc_Yld_TSS_pred()` coefficients match Beusen Table 2 at full precision (published: intercept 0.55, 0.03, 0.12, 0.29, 0.10; lithology 1.04 / 0.90 / 0.31 / −2.07 / 0.67 — code: 1.0379 / 0.9037 / 0.3125 / −2.0678 / 0.6694) | identical, incl. truncation ranges and cap | ✓ coefficients verified; class-5 coefficient is an assumption (F5) |
| Eq. 26 / Beusen Eq. 2 | TSSpc_POC = −0.160·log₁₀(TSSc)³ + 2.83·log₁₀(TSSc)² − 13.6·log₁₀(TSSc) + 20.3, TSSc = Yld_TSS,pred/Rnat (mg/L) | polynomial reproduced verbatim; **the flattening Ludwig et al. drew beyond TSSc ≈ 2250 mg/L is not implemented** | identical | ✗ known unphysical branch left active (F1) |
| Eq. 28 / Beusen Eq. 4 | ln(1000·Ld_PP,pred) = −3.096 + 1.002·ln(1000·Ld_POC,pred), loads in Mg/yr in the paper's notation (×1000 → kg/yr) | code: exp(−3.09614 + 1.001769·ln(Ld_POC in kg/yr))/1000; pre-dam POC load recovered as Ld_TSS/(1−D_TSS) with the kg conversion explicit | identical | ✓ unit chain verified; constants differ from Beusen's *published rounded* values (−3.098, 1.002) by <0.2 % in output (F7) |
| Eq. 29 / 30 | Yld_F = (1 − D_TSS) · Yld_F,pred | PP computed on loads then converted to yield (1000·Ld/A) — algebraically identical | identical | ✓ |
| — | TP = DIP + DOP + PP | `allformmodel.py`: TPyld = DIPYld + DOPYld + PPYld | identical sum | ✓ |

Calibration constants (`constants.cfg`, designated source of record by
the manual): DIP a=0.85, b=2, c=1, e=0.29, EC=26; DOP a=0.95, c=0.010,
e=0.010, EC=15. Both implementations use them unchanged. The paper's
§6 narrative is consistent with them (for DIP only e_F was recalibrated
for NEWS 2; DON/DOP parameters were fully recalibrated), but the values
themselves sit in Supplementary Table A2, which this audit could not
independently re-read (F9).

## Findings

Severity: **High** = changes published-style results materially on real
inputs. **Medium** = changes results for a identifiable subset, or
produces impossible values on plausible inputs. **Low** = latent,
cosmetic, or negligible on current inputs.

| # | Severity | Class | Location | Finding and evidence | Recommended action |
|---|---|---|---|---|---|
| F1 | **High** | Original-code error (code–paper mismatch) | `globalnews.py`, `model_particulate()` (TSSpc_POC block); port `_particulate_pp()` faithfully reproduces it | The Ludwig %POC polynomial has a minimum at TSSc ≈ 2250 mg/L (%POC ≈ 0.48 %) and *rises* thereafter — up to 6.0 % at TSSc = 10⁶ mg/L. The code's own comment states the upturn "makes no physical sense" and that Ludwig et al. drew a horizontal line beyond that point, then declines to implement the flattening on the grounds that "the Yld_TSS_pred_MAX implementation makes that unnecessary." That reasoning is wrong: TSSc = Yld_TSS,pred/**Rnat**, so the 5000 cap on yield does not bound concentration. Probe on the bundled c00 inputs: **580 of 4,286 live basins (13.5 %) have TSSc_pred > 2250 mg/L, up to 301,580 mg/L**, where %POC is inflated up to ~9× relative to the flattened value, and PP with it (PP ∝ POC^1.002). | In any ensemble use, apply the flattening (%POC held at its 2250 mg/L minimum for higher TSSc) as a documented variant and report PP both ways; do not silently inherit the rising branch. |
| F2 | Medium | Original-code error (code–paper mismatch) | `globalnews.py`, `Calc_FEriv()`; port `_dissolved()` | D_DIP is capped at 0.85 in code (`D_F[D_F > 0.85] = 0.85`; DIN is capped at 0.965). No cap appears in the paper's Eq. 17 or §3.2.2. The cap is active on real inputs: the bundled D_DIP column reaches 1.0, with **26 basins at exactly 1.0 and 62 at ≥ 0.85**; for those basins the code overrides a data-derived total retention and guarantees ≥ 15 % × (1−FQrem) transmission. | Treat the cap as a model assumption to be declared wherever NEWS2 numbers are used; in a Bayesian refit, replace the hard cap with the retention function's own asymptote or a prior. |
| F3 | Medium | Original-code error | `globalnews.py`, `Calc_FEriv()` and `Calc_WSExportFunc()`; port identical | The fraction guards FEriv ≤ 1 and FEws ≤ 1 are enforced **for DIN only**; the code comments admit the inconsistency ("logically it should be <= 1. However, that's not how the other models were originally defined"). Probes on the port: FQrem = 1.2 → FEriv = −0.20 and **negative DIP yield (−5.53 kg km⁻² yr⁻¹)**; FQrem = −0.1 or D_DIP = −0.2 → FEriv = 1.1–1.2 (export exceeds inputs); D_TSS = 1.5 → TSS yield −25.9 and PP yield −37.0 kg km⁻² yr⁻¹. The bundled c00 inputs all lie in [0, 1], so the published run is unaffected — but nothing warns a user whose pre-processing produces out-of-range fractions. | Add input-domain validation (fractions in [0, 1]) at load time in the port's corridor workflow; report violations rather than clipping silently. |
| F4 | Low | Original-code error (benign) | `globalnews.py`, `Calc_WSExportFunc()` DIP branch; port identical | The near-zero floor `Rgt0[Rnat < 0.00001] = 0.0000001`, commented as preventing the negative exponent from "blowing up", is **dead code**: every Rnat below 10⁻⁵ m/yr is also below R_MIN = 3×10⁻³ m/yr, where f is zeroed anyway, so the floor can never affect any output. (In the original, `Rgt0 = Rnat` also mutates the local runoff array in place — harmless here because the array is a fresh temporary.) The misleading comment invites future editors to remove the wrong safeguard. | Leave behaviour unchanged; correct the comment in any derived code to credit the arid mask as the operative guard. |
| F5 | Low | Original-code error (undocumented assumption) | `globalnews.py`, `Calc_Yld_TSS_pred()` lithology map; port `LITHO_CHEM_COEF` | Beusen et al. (2005, Table 2 note) state ice and **gabbros did not occur** in the 124-river development set, so lithology class 5 (gabbros) has no published coefficient. The code silently assigns class 5 the sand/sandstone reference value 0.0 (as for class 1), while class 8 is zeroed outright. No class-5 basin exists in the bundled c00 inputs (classes present: 1, 2, 3, 4, 6, 7, 8), so the assumption is untested by the reference run. | Flag class-5 basins explicitly in any new application domain and treat their TSS/PP as assumption-driven. |
| F6 | Low | Original-code error (missing input validation) | `globalnews.py`, `model_dissolved()` Agfr split; bundled input `c00_diffsrc_landuse_ng.csv` | `agric` reaches **100.01 % in 32 bundled basins**, so (1 − Agfr) < 0 and the natural-area export-coefficient source RSdif_ec_nat is slightly **negative** — an impossible source term, produced by the input table and passed through unchecked by both implementations. Magnitude is negligible (≤ 0.01 % of EC_F), but the sign error is real. | Clip Agfr to [0, 1] in pre-processing and log corrections; add a range check on all percentage inputs. |
| F7 | Low | Paper ambiguity (rounding) | Beusen et al. (2005) Eq. 4 vs `globalnews.py`, `model_particulate()` PP block | Published constants are ln PP = −3.098 + 1.002·ln POC; the code uses −3.09614 and 1.001769 (presumably full regression precision). Output difference < 0.2 % across the data range. The same pattern holds for the TSS coefficients (code carries more digits than Beusen Table 2). | None — record that the code, not the paper's rounded table, is the precision source of record. |
| F8 | Medium | Documented calibration/design choice (sharp edge) | Paper §3.1.4, §4.1; `globalnews.py` R_MIN masks in both sub-models; port identical | Arid zeroing at Rnat < 3 mm/yr is a **hard step**: probes show PP yield jumping 0 → 31.6 kg km⁻² yr⁻¹ and TSS yield 0 → 51.9 Mg km⁻² yr⁻¹ across a 0.01 mm/yr runoff increment at the threshold, and the step can be as large as the 5000 Mg km⁻² yr⁻¹ cap. The dissolved forms step too, but by ≲ 0.07 kg km⁻² yr⁻¹ because f_F(3 mm/yr) is tiny. 739 bundled basins are arid by this rule. A basin whose runoff estimate straddles 3 mm/yr flips between "no particulate export" and "full regression export". | In corridor use, report sensitivity to the runoff threshold for near-threshold basins; consider a smooth transition only as a declared variant, never silently. |
| F9 | Low | Paper ambiguity (documentation) | Mayorga et al. (2010) Supplementary Tables A2/A3 | The calibration constants and the TSS truncation ranges are specified only in the online supplement. The truncation ranges in code (MargGrass 0–84.6865; WetlnRice 0–24.7206; FrnrPrecp 0.854985–11.1682; FrnrSlope 1.1301–38.4531) are consistent with their stated provenance (min/max of the Beusen development set) but could not be verified digit-for-digit from the main paper. | Obtain the supplement before any recalibration work; until then treat `constants.cfg` + code as the operative specification, as the manual directs. |
| F10 | Info | Documented calibration/design choice | Paper §3.1.4; `Calc_PointSources()` path | Point sources export from arid basins: in probes, a basin with Rnat = 0 still exports DIP 12.0 and DOP 0.12 kg km⁻² yr⁻¹ (the sewage terms at FEriv = 1); 40 bundled arid basins have DIP yield > 0. This is deliberate (population centres near coasts) and stated in the paper. | None — but ensemble members that zero all export at zero runoff will differ from NEWS2 exactly here; expect and explain the difference. |
| F11 | Info | Documented calibration/design choice | Paper §3.1.2, §3.1.3; constants.cfg [CAL.DOP] | DOP's f = Rnat^0.95 is unbounded and FEws,DOP carries no ≤ 1 cap; the paper notes FEws "never exceeds 1" in practice — true only because e_DOP = 0.010 (FEws would exceed 1 only at Rnat > ~127 m/yr). The export-coefficient term f·EC grows without bound in runoff (at the data maximum Rnat = 6.0 m/yr it contributes ≈ 82 kg km⁻² yr⁻¹ of DOP before the Agfr split). | None for global use; at corridor scale, check that local runoff stays in the calibrated range before trusting DOP. |
| F12 | Info | Documented calibration/design choice | Paper §3.1.2 (negative balances set to zero); §4.1 (5000 cap) | WSdif_ant floored at 0 (probe: fert+manure 15, crop export 100 → diffuse explicit term 0, yield from EC + point terms only) and the TSS 5000 Mg km⁻² yr⁻¹ cap (≈ 2× the development-set maximum) are both explicit in the paper and faithfully coded. | None. |
| F13 | Low | Port difference (not an error) | Port `_particulate_pp()`, `compute_phosphorus()` | Two deliberate robustness differences from the original, neither affecting the verified reference run: (i) the port applies `nan_to_num` to PP loads, where the original relies on the NoFlux mask to rescue NaNs (e.g. D_TSS = 1 gives 0/0 in the pre-dam POC recovery in both; both end at PP = 0); (ii) a zero-area basin yields PP_Yld = NaN in the port (1000·Ld/A) — the original has the same unguarded division, and no zero-area basin exists in the bundled data (min 362.9 km²). The port also truncates float lithology classes to int, where the original would raise on a non-integer class. | Add an A > 0 guard if the port is ever pointed at user-assembled basin tables. |

**Counts.** Port errors: **0**. Original-code errors / code–paper
mismatches: **5** (F1, F2, F3, F4, F5) plus **1** missing-input-validation
finding shared with the input table (F6). Paper ambiguities: **2**
(F7, F9). Documented calibration/design choices with sharp edges or
ensemble consequences: **4** (F8, F10, F11, F12). Port differences
requiring no action: **1** (F13).

**Unit audit (all verified by hand).** Yields kg km⁻² yr⁻¹; loads
Ld = Yld·A/1000 → Mg yr⁻¹ ✓. Runoff mm/yr → m/yr at entry ✓. TSSc =
(Mg km⁻² yr⁻¹)/(m yr⁻¹) = mg/L ✓ (paper §4.2 states the same identity).
PP chain: Ld_TSS in Mg/yr, /(1−D_TSS) recovers the pre-dam load,
×1000 → kg/yr for Beusen Eq. 4, exp(·)/1000 back to Mg/yr ✓ — the
×1000 inside the paper's Eq. 28 logarithms is the same conversion, not
an extra factor. DIP mass-balance probe: modelled DIP yield never
exceeds gross explicit inputs + EC term + point sources (max ratio
0.986 over 6,081 basins) ✓. No negative yields occur on the bundled
inputs for any form ✓.

## What the literature says — and what this audit can confirm

- **Global totals / retention method.** A global retention-model
  intercomparison (*Frontiers in Water*, 2022) notes NEWS2's TP export
  (9 Tg yr⁻¹) is almost double IMAGE-GNM's (~4 Tg yr⁻¹), and attributes
  the spread mainly to hydrological inputs, spatial resolution and the
  retention method, with P retention ideally differing between lakes
  and rivers. *Audit position:* consistent with what we find
  structurally — NEWS2's P retention is a single pre-computed,
  discharge-weighted reservoir term per basin (plus the undocumented
  0.85 cap, F2), with no in-channel term for P at all (L_P = 0). We
  cannot test the global-total discrepancy from the code alone.
- **MARINA comparisons (Strokal et al.).** For Chinese basins, MARINA
  reports model error 8 % against NEWS2's 26 %, and NEWS2 DIN/DIP
  Nash–Sutcliffe of 0.42 on five large Chinese rivers — i.e. the
  basin-aggregated, globally calibrated formulation loses skill at
  regional scale, which is MARINA's stated reason for downscaling to
  sub-basins with more explicit sources. *Audit position:* we confirm
  the mechanism (one global constant set, basin-lumped inputs, no
  sub-basin structure) but cannot reproduce their error statistics
  without their observation set.
- **Particulate sub-model skill.** Beusen et al. (2005) themselves
  bound the TSS model: predictions within a factor of 4 for 80 % of
  rivers, and "should not be used to predict the fluxes … for the
  mouth of individual rivers." *Audit position:* confirmed as a
  design limit; F1 shows an additional, purely numerical reason for
  caution at the high-concentration end that the factor-of-4 statement
  does not cover.
- **Known-issue search.** No published erratum or code correction for
  NEWS2 was found; the %POC upturn (F1) and the retention caps (F2)
  appear to be documented only inside the code's own comments.

## What NEWS2 teaches the ensemble

1. **Make every retention term explicit, named and inspectable.**
   NEWS2's FEriv = (1−L)(1−D)(1−FQrem) decomposition is its best
   design feature: each loss process is a separate factor with its own
   data provenance, so a reviewer can challenge one term at a time.
   Its worst feature is the opposite — a cap (F2) that exists only in
   code. Ensemble rule: every factor and every cap lives in the
   published equation set; if it is in the code, it is in the paper.
2. **Separate explicit budgets from calibrated net terms, and say
   which is which.** The WSdif explicit balance (fertilizer + manure −
   crop export, floored at zero) carries process meaning and scenario
   leverage; the export-coefficient term is a calibrated residual
   wearing a process costume (for DIP it is labelled "weathering").
   Keeping the two apart — and the floor rule visible — is what makes
   source attribution (Eqs. 20–23) possible at all.
3. **A shared runoff-modulation function is a strong unifying device
   — and a shared point of failure.** One f_F(Rnat) drives both the
   explicit-export fraction and the export-coefficient term, which
   keeps forms consistent; but the arid step at 3 mm/yr (F8) then
   propagates into every form at once. Ensemble members should state
   their threshold behaviour as a tested property, not an emergent
   surprise.
4. **Chained empirical sub-models multiply edge cases.** PP is three
   regressions deep (TSS → %POC → PP), and the one genuine numerical
   defect found (F1) sits at the junction of two of them, in a regime
   the original authors noticed, wrote a comment about, and talked
   themselves out of fixing with a faulty argument. Junction regimes —
   where one sub-model's output leaves the next sub-model's
   calibration range — deserve explicit tests in every ensemble member.
5. **Calibration constants are part of the model, and they age.**
   NEWS2's constants live in a config file the manual designates as
   the source of record, yet the paper warns the released code is not
   the version that produced its results. An ensemble that refits
   members (e.g. in Stan) should version constants with the same rigour
   as code, and should expect globally calibrated constants (year-2000,
   in NEWS2's case) to act as informative priors at corridor scale, not
   as truths.
6. **Design for the fractions to be fractions.** NEWS2's DIN-only
   guards (F3) show what happens when range discipline is applied to
   one form and not its siblings: nothing, on the calibration data —
   and silent negative exports on the first out-of-domain input.
   Uniform input-domain validation is cheaper than one corrupted
   scenario run.

*Audit prepared from the published papers, the bundled GNE manual, the
original GNE code and the port, with numeric probes run on the port
against synthetic edge-case basins and the bundled c00 input tables.
No code was changed.*
