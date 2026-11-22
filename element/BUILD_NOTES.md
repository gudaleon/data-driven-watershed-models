# ELEMeNT-P member — build notes

## Provenance: an equations build, not the authors' code

The ELEMeNT / ELEMeNT-P code of Van Meter, Basu and co-workers is **not
open source** (verified 2026-11-21: the University of South Florida
GAIM/WBMP Model Library entry for ELEMeNT-N gives its download link as
"Contact Dr. Van Meter", and no repository or deposit exists). The files
here re-implement the published equations in Python:

- `element_p.py` — the model (one commented code base).
- `run_gates.py` — verification gates and machinery drives; writes
  `gate_results.json` and `test_outputs.csv`.

Nothing in this folder is committed to git.

## What was built

A catchment-aggregated, annual-timestep ELEMeNT-P. The published model
(Van Meter et al. 2021, *Water Resources Research* 57, e2020WR029316 —
"VM21") is distributed over a distribution of land-use trajectories
(its Eq. 3) with land-use conversion step functions, an empirically
parameterised erosion model, manure-handling losses and wastewater
inputs. This build aggregates to a single effective land unit — the
ensemble member's role is the legacy-store and lag structure, not
land-use bookkeeping — and keeps VM21's process chain:

| Step | Source |
|---|---|
| Annual P surplus Np = inputs − crop/grass removal drives the model | VM21 Eqs. 4–6 |
| Surplus split by humification coefficient h into protected (Mp) and active (Ma) soil organic P pools | VM21 Eqs. 7–8 |
| First-order mineralisation ka·Ma + kp·Mp into the mineral pool | VM21 Eq. 10, terms 1–2 |
| Mineral pool Mmin = Mdis + Msorb,r + Msorb,i; dissolved/reversibly-sorbed split from the Freundlich isotherm (bisection; the paper uses Newton–Raphson, same root) | VM21 Eqs. 9, 12, 13 |
| Irreversible precipitation of reversibly sorbed P at rate ksorb | VM21 §2.1.4 |
| Leaching flux Js = Mdis · min(Q/Vw, 1) out of the source zone | VM21 Eq. 11 |
| Gross erosion of organic + sorbed P at rate keros | VM21 Eq. 16 |
| Surface delivery after network retention (1 − Kret); the retained share is kept as a network store | VM21 Eq. 15 |
| Subsurface transport as a convolution of Js with a parametric travel-time distribution f(τ), with first-order sorption attenuation exp(−γτ); the attenuated mass is kept as a groundwater-sorbed store so the ledger closes | VM21 Eq. 1; Van Meter & Basu 2017 supplement Eq. S2.1 |
| Outlet load Mout = subsurface outflow + surface delivery + optional wastewater series | VM21 Eq. 1 |

Stores are explicit in every output year: Mp, Ma, Mdis, Msorb,r,
Msorb,i, the age-structured transit cohorts (summed and inspectable),
gw_sorbed, network, and a deficit ledger for unmet negative-surplus
drawdown. Basu et al. 2025 (*Nature Water*, wetland N removal) was used
as mechanism support only; **no nitrogen parameter or number from it
enters this phosphorus build**, and Basu et al. 2010 supplies the
stationarity behaviour checked in Gate 3, not parameters.

## Switches (choices the papers leave open at this aggregation)

- `ttd`: `gamma` (default; shape `ttd_shape`, mean `ttd_mean`),
  `exponential` (gamma, shape 1), or `triangular` — the distribution
  the Van Meter & Basu 2017 supplement itself uses for its lag
  simulations (mean L, sd 1 yr). Both readings are implemented because
  the papers support both.
- `mineral_partition`: `freundlich` (default, VM21 Eqs. 12–13) or
  `linear` (fixed dissolved fraction of the mobile mineral pool).
- `negative_surplus`: `drawdown` (default — a negative annual surplus
  mines mineral, then active, then protected stores, with any unmet
  remainder in the deficit ledger so the balance closes exactly) or
  `clamp` (truncate at zero; provided for comparison and flagged as
  not mass-conserving against the raw trajectory).
- `surface_pathway`: on (default) or off (erosion zeroed, isolating
  the subsurface legacy pathway).

## Parameters — an open limitation, stated plainly

VM21's calibrated parameter values are in its supplement Table S10,
which this project does not hold. The defaults in `element_p.py` are
therefore **named, plausible mid-range machinery defaults, not the
paper's calibrated set**: h = 0.30, ka = 0.10 /yr, kp = 0.01 /yr,
keros = 0.001 /yr, ksorb = 0.05 /yr, γ = 0.05 /yr, Kret = 0.30,
Vw = 250 mm, Q = 120 mm/yr (subsurface drainage through the source
zone, of the order of southern Ontario recharge — not total
discharge), gamma TTD with mean 30 yr and shape 2, and a scaled
Freundlich kf = 2×10⁵ with nfr = 0.5 chosen so the dissolved share of
the mobile mineral pool is a few percent, the regime of agricultural
soils in VM21. Every result below is a machinery/behaviour result
under these defaults. Corridor application requires catchment
calibration against the load record.

## Verification gates (run_gates.py, gate_results.json)

**Gate 1 — mass balance: PASS on every run.** Per-year closure
residuals are at most 4×10⁻¹³ kg P/ha (tolerance 1×10⁻⁹) and
cumulative residuals at most ~10⁻¹², on the Canada drive, the Ontario
drive, the Gate-2 history and the 2,000-year stationarity run: in
every year, cumulative surplus equals cumulative export plus the
change in all stores, the deficit ledger included.

**Gate 2 — qualitative reproduction of VM21: PASS as a behaviour
check.** Driven by a Grand-River-style surplus history (3 kg/ha/yr in
1900 rising to the paper's ~15 kg/ha/yr peak in 1978, then declining
to 6 by 2020; cumulative surplus 976.5 kg P/ha), the model exports
109.5 kg P/ha cumulatively, a **retained share of 0.888** under the
gamma TTD and **0.919** under the triangular (supplement) branch,
against the paper's calibrated finding that ~96% of net inputs since
1900 were retained in the Grand River Watershed. Export lags surplus:
the surplus peaks in 1978, the export peaks in 2018, and the
cross-correlation lag is 30 years. The signature (lag, late export
peak, large retained share) reproduces; the retained share sits below
the paper's 0.96 because the defaults are uncalibrated — this is a
behaviour check on the machinery, not a calibration claim, and the
gap is recorded rather than tuned away.

**Gate 3 — stationarity (Basu et al. 2010): PASS.** Under a constant
surplus of 8 kg/ha/yr for 2,000 years, the dynamic stores (source zone
plus transit) approach a steady state: their relative change over the
final 50 years is 0.76%. Export alone reaches 0.51 of the input rate
in the final decade and export plus the permanent retention sinks
(pathway sorption and network retention, which are cumulative stores
in VM21's accounting) reaches 0.90, with the remainder still entering
the slowly filling irreversibly-sorbed pool — the transport-limited,
legacy-dominated regime Basu et al. describe, in which loads are
sustained by the accumulated store long after inputs stabilise.

## Machinery drives on held trajectories (not corridor applications)

- **Canada, census trajectory v2, 2010–2021** (surplus = total P
  inputs − Plant Nutrition Canada workbook crop removal;
  surplus range 3.2–10.4 kg/ha/yr): final-year load 0.73 kg/ha/yr,
  source-zone store 676 kg P/ha after a 200-year spin-up at the
  first-year surplus.
- **Ontario, long trajectory v3, 1950–2007** (gross-basis balance per
  hectare of land in crops; surplus range 5.8–21.2 kg/ha/yr):
  final-year load 1.93 kg/ha/yr, source-zone store 1,295 kg P/ha.

Both drives exist to exercise the machinery on real surplus
trajectories at national/regional scale. Neither is a catchment
application, and no load from them is a result for the proposal.

## Data still needed for corridor application

1. **Catchment-scale NAPI trajectories**: annual P surplus (fertilizer
   + manure + deposition + land-applied domestic P − crop and grass
   removal, VM21 Eqs. 4–6) for each corridor catchment over the
   legacy-accumulation decades, built from the Census of Agriculture
   record — the separate corridor NAPI build, not yet done.
2. **Travel-time parameters by catchment**: a mean (and shape) for
   f(τ) per catchment, from catchment attributes or calibration
   against the annual load record; the ensemble design's lag analysis
   depends on these, and a single regional default is not defensible
   per catchment.
3. **Source-zone initialisation per catchment**: measured or modelled
   initial soil P stores (VM21 initialises from soil P data); the
   spin-up convention used in the machinery drives is a placeholder.
4. **Calibration targets**: the annual load record assembled for the
   SPARROW backbone, plus soil P accumulation data where held (VM21's
   third objective function), and the paper's own parameter ranges
   (supplement Table S10) if the supplement can be obtained.
5. Optional refinements the aggregated build omits: land-use
   trajectory distributions and conversion step functions (VM21
   Eqs. 3, 7–8), discharge-varying keros and Kret, manure-handling
   losses (VM21 §2.2), and wastewater inputs (VM21 Eq. 1, third term)
   where a catchment's point sources warrant them.
