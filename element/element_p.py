"""
ELEMeNT-P (catchment-aggregated) -- an equations build from the published papers.

This is NOT the authors' code. The ELEMeNT / ELEMeNT-P code of Van Meter,
Basu and co-workers is not open source (verified 2026-11-21: the University
of South Florida GAIM/WBMP Model Library entry for ELEMeNT-N lists the
download link as "Contact Dr. Van Meter"). This module re-implements the
published equations in Python so the ensemble has a runnable legacy-store
member. Every equation is traced in comments to its source:

  VM21  Van Meter, K.J. et al. (2021), "Beyond the Mass Balance: Watershed
        Phosphorus Legacies...", Water Resources Research 57,
        e2020WR029316. DOI 10.1029/2020WR029316.  (ELEMeNT-P)
  VM17  Van Meter, K.J. & Basu, N.B. (2017), Environ. Res. Lett. 12,
        084017, and its supplement (ELEMeNT-N structure; supplement Eq.
        S2.1 is the travel-time convolution used for the lag simulations).
  B10   Basu, N.B. et al. (2010), Geophys. Res. Lett. 37, L23404
        (biogeochemical stationarity: legacy stores make annual loads
        transport-limited, with temporally near-invariant flow-weighted
        concentrations).

SCOPE OF THIS BUILD (catchment-aggregated, annual timestep)
-----------------------------------------------------------
The published ELEMeNT-P (VM21) is distributed over a distribution of
land-use trajectories LU(s,t) (VM21 Eq. 3) and includes land-use
conversion step functions, a HYPE-style mineral-P module with a
Freundlich equilibrium solved by Newton-Raphson, an empirically
parameterised erosion model, manure-handling losses and wastewater
inputs. This build aggregates the watershed to a single effective
land unit (the ensemble member's role is the legacy/lag structure, not
land-use bookkeeping) and keeps the process chain of VM21:

  annual P surplus Np (VM21 Eqs. 4-6, supplied by the caller as a
  NAPI-style trajectory: inputs minus crop/grass removal)
    -> soil organic P: protected pool Mp and active pool Ma, split by
       the humification coefficient h, first-order mineralisation
       ka (active) and kp (protected)          [VM21 Eqs. 7-8]
    -> soil mineral P Mmin = Mdis + Msorb,r + Msorb,i, with the
       dissolved / reversibly-sorbed split from the Freundlich
       isotherm                              [VM21 Eqs. 9, 12, 13]
    -> leaching flux Js out of the source zone, proportional to the
       dissolved pool and the annual flushing ratio Q/Vw
                                               [VM21 Eq. 11]
    -> subsurface transport as a convolution of Js with a parametric
       travel-time distribution f(tau), with first-order sorption
       attenuation exp(-gamma * tau) along the pathway
                                               [VM21 Eq. 1; VM17 Eq. S2.1]
    -> surface pathway: erosion of organic + sorbed P at rate keros
       [VM21 Eq. 16], delivered after network retention Kret
                                               [VM21 Eq. 15]
    -> outlet load Mout = subsurface outflow + surface delivery
       (+ wastewater term DMSTCsw of VM21 Eq. 1, carried as an
       optional caller-supplied series, default zero).

STORES (all explicit in the output, kg P ha-1):
  Mp, Ma            soil organic P (protected / active)   [VM21 sec. 2.1.3]
  Mdis, Msorb_r,
  Msorb_i           soil mineral P (dissolved / reversibly sorbed /
                    irreversibly precipitated)            [VM21 Eq. 9]
  transit           P in subsurface transit (age-structured cohorts)
  gw_sorbed         P sorbed along subsurface pathways (the
                    exp(-gamma tau) sink of VM21 Eq. 1, kept as a
                    store so the mass balance closes exactly)
  network           P retained in the river network (Kret share of
                    the surface flux, VM21 Eq. 15)

SWITCHES (choices the papers leave open at this aggregation level;
both branches implemented where the papers support two readings --
see BUILD_NOTES.md):
  ttd               "gamma" (default), "exponential" (gamma with shape
                    1) or "triangular" (the VM17 supplement's choice
                    for its lag simulations, mean L, sd 1 yr).
  mineral_partition "freundlich" (default; VM21 Eqs. 12-13) or
                    "linear" (fixed dissolved fraction fdiss of the
                    mobile mineral pool).
  negative_surplus  "drawdown" (default: a negative annual surplus
                    mines the soil stores, in the order mineral,
                    active, protected, with any unmet remainder in a
                    deficit ledger so the balance still closes) or
                    "clamp" (negative surplus truncated to zero --
                    provided for comparison, NOT mass-conserving by
                    construction and flagged as such in the output).
  surface_pathway   True (default) / False (erosion terms zeroed;
                    isolates the subsurface legacy pathway).

PARAMETERS: VM21's calibrated values live in its supplement Table S10,
which is not held in this project; the defaults below are therefore
named, plausible mid-range values, NOT the paper's calibrated set.
They are machinery defaults for ensemble plumbing and behaviour
tests, and BUILD_NOTES.md says so. Corridor application requires
catchment calibration (see BUILD_NOTES.md, "Data still needed").

Units: kg P ha-1 for stores and kg P ha-1 yr-1 for fluxes, annual step.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import stats


# ----------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------
@dataclass
class ElementPParameters:
    # VM21 Eq. 7-8: humification split of the annual surplus between the
    # protected and active organic pools (Janssen 1984; Six et al. 2002,
    # as cited in VM21 sec. 2.1.3).
    h: float = 0.30
    # VM21 Eqs. 7-8: first-order mineralisation rate constants [1/yr].
    ka: float = 0.10   # active pool
    kp: float = 0.01   # protected pool
    # VM21 Eq. 16: gross erosion rate constant [1/yr] (in the paper an
    # empirical function of land use and discharge, de Vente & Poesen;
    # here a constant -- see BUILD_NOTES).
    keros: float = 0.001
    # VM21 sec. 2.1.4: first-order irreversible precipitation of newly
    # adsorbed P [1/yr].
    ksorb: float = 0.05
    # VM21 Eqs. 12-13: Freundlich constants. Ceq = (Msorb_r / kf) ** nfr
    # with masses in kg P ha-1 and Vw in mm, so Ceq is a scaled
    # concentration; kf carries the compensating units.
    # kf is a SCALED constant (Ceq here is in kg P ha-1 per Vw, not
    # mg/L), chosen so the dissolved share of the mobile mineral pool
    # is a few percent, the regime of agricultural soils in VM21 (the
    # paper's own kf is in its unheld supplement, Table S10).
    kf: float = 200000.0
    nfr: float = 0.50
    # VM21 Eq. 11: source-zone water volume Vw = n * s * V [mm]; with a
    # 1 m source zone (VM21 sec. 2.1), porosity ~0.45 and mean
    # saturation ~0.55, Vw ~ 250 mm.
    Vw: float = 250.0
    # Annual water flux THROUGH the source zone (subsurface drainage /
    # recharge, not total stream discharge) Q [mm/yr] for VM21 Eq. 11;
    # may be overridden by a series in run(). Southern Ontario recharge
    # is of order 100-150 mm/yr.
    Q: float = 120.0
    # VM21 Eq. 15: network (lake/reservoir/riparian) retention share of
    # the surface flux.
    Kret: float = 0.30
    # VM21 Eq. 1 / VM17 Eq. S2.1: first-order sorption attenuation along
    # subsurface pathways [1/yr].
    gamma: float = 0.05
    # Travel-time distribution: mean [yr] and (gamma) shape [-].
    ttd_mean: float = 30.0
    ttd_shape: float = 2.0
    # Switches (see module docstring).
    ttd: str = "gamma"
    mineral_partition: str = "freundlich"
    negative_surplus: str = "drawdown"
    surface_pathway: bool = True
    # Linear-partition switch branch: dissolved fraction of the mobile
    # mineral pool.
    fdiss: float = 0.10
    # Maximum tracked transit age [yr]; residual cohort mass is released
    # at this age so no mass is stranded.
    max_age: int = 500


# ----------------------------------------------------------------------
# Travel-time distribution helpers
# ----------------------------------------------------------------------
def ttd_weights(par: ElementPParameters) -> np.ndarray:
    """Annual release weights w(a), a = 0..max_age-1: the probability
    that a unit of leached P entering the subsurface in year t leaves
    in year t+a, BEFORE sorption attenuation (the attenuation is
    applied to the cohorts themselves, year by year).

    VM21 Eq. 1 convolves the source function with f(tau); VM17's
    supplement (Eq. S2.1) does the same for nitrogen and uses a
    triangular f(tau) for its lag simulations. The gamma family is the
    standard parametric choice for catchment travel times and nests
    the exponential at shape 1.
    """
    a = np.arange(par.max_age + 1, dtype=float)
    if par.ttd == "gamma" or par.ttd == "exponential":
        shape = 1.0 if par.ttd == "exponential" else par.ttd_shape
        scale = par.ttd_mean / shape
        cdf = stats.gamma.cdf(a + 1.0, shape, scale=scale)
        cdf_prev = stats.gamma.cdf(a, shape, scale=scale)
    elif par.ttd == "triangular":
        # VM17 supplement S2.1.3: triangular, mean L, sd 1 yr.
        L = par.ttd_mean
        lo, hi = L - np.sqrt(6.0), L + np.sqrt(6.0)  # sd = (hi-lo)/sqrt(24) = 1
        cdf = stats.triang.cdf(a + 1.0, 0.5, loc=lo, scale=hi - lo)
        cdf_prev = stats.triang.cdf(a, 0.5, loc=lo, scale=hi - lo)
    else:
        raise ValueError(f"unknown ttd switch: {par.ttd}")
    w = cdf - cdf_prev
    w = w / w.sum()  # lump the tail beyond max_age into the last year
    return w


# ----------------------------------------------------------------------
# Freundlich equilibrium (VM21 Eqs. 12-13)
# ----------------------------------------------------------------------
def freundlich_split(mobile: float, par: ElementPParameters) -> tuple[float, float]:
    """Split the mobile mineral pool (Mdis + Msorb_r) into dissolved and
    reversibly sorbed parts.

    VM21 Eq. 12:  Ceq = (Msorb_r / kf) ** nfr
    VM21 Eq. 13:  Mdis = Ceq * Vw
    => Msorb_r = kf * (Mdis / Vw) ** (1 / nfr).  Solved by bisection on
    Mdis in [0, mobile] (the paper uses Newton-Raphson; the root is the
    same).
    """
    if mobile <= 0.0:
        return 0.0, 0.0
    if par.mineral_partition == "linear":
        return par.fdiss * mobile, (1.0 - par.fdiss) * mobile
    if par.mineral_partition != "freundlich":
        raise ValueError(f"unknown mineral_partition switch: {par.mineral_partition}")

    def resid(mdis: float) -> float:
        return mdis + par.kf * (mdis / par.Vw) ** (1.0 / par.nfr) - mobile

    lo, hi = 0.0, mobile
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if resid(mid) > 0.0:
            hi = mid
        else:
            lo = mid
        if hi - lo < 1e-14 * max(1.0, mobile):
            break
    mdis = 0.5 * (lo + hi)
    return mdis, mobile - mdis


# ----------------------------------------------------------------------
# The model
# ----------------------------------------------------------------------
@dataclass
class ElementPResult:
    years: np.ndarray
    surplus: np.ndarray            # Np actually applied (after clamp switch)
    surplus_raw: np.ndarray        # Np as supplied
    load_subsurface: np.ndarray    # outlet load, subsurface pathway
    load_surface: np.ndarray       # outlet load, surface pathway (post-retention)
    load_total: np.ndarray         # Mout (VM21 Eq. 1, wastewater = supplied series)
    leaching: np.ndarray           # Js out of the source zone (VM21 Eq. 11)
    erosion_gross: np.ndarray      # SURFeros (VM21 Eq. 16)
    mineralisation: np.ndarray     # ka*Ma + kp*Mp (VM21 Eq. 10, terms 1-2)
    stores: dict = field(default_factory=dict)   # year-end store trajectories
    balance_residual: np.ndarray = None          # per-year closure error
    deficit: np.ndarray = None                   # unmet drawdown ledger


def run_element_p(years, surplus, par: ElementPParameters | None = None,
                  Q_series=None, wastewater=None,
                  initial: dict | None = None) -> ElementPResult:
    """Run ELEMeNT-P (aggregated) on an annual surplus trajectory.

    years      : array of calendar years (annual, consecutive)
    surplus    : Np(t) [kg P ha-1 yr-1], VM21 Eqs. 4-6 (inputs minus
                 crop/grass removal; may be negative)
    Q_series   : optional annual flushing Q(t) [mm/yr] for VM21 Eq. 11
    wastewater : optional DMSTCsw(t) series, VM21 Eq. 1 (direct outlet
                 input, not stored); default zero
    initial    : optional initial stores dict with keys Mp, Ma, Mdis,
                 Msorb_r, Msorb_i (default all zero: a pristine
                 catchment; spin-up is the caller's job, as in VM21,
                 where the source zone is initialised from measured
                 soil P)
    """
    par = par or ElementPParameters()
    years = np.asarray(years, dtype=int)
    surplus_raw = np.asarray(surplus, dtype=float)
    n = len(years)
    if len(surplus_raw) != n:
        raise ValueError("surplus must match years")
    Q = (np.full(n, par.Q) if Q_series is None
         else np.asarray(Q_series, dtype=float))
    ww = (np.zeros(n) if wastewater is None
          else np.asarray(wastewater, dtype=float))

    init = dict(Mp=0.0, Ma=0.0, Mdis=0.0, Msorb_r=0.0, Msorb_i=0.0)
    if initial:
        init.update(initial)
    Mp, Ma = init["Mp"], init["Ma"]
    Mdis, Msorb_r, Msorb_i = init["Mdis"], init["Msorb_r"], init["Msorb_i"]

    w = ttd_weights(par)
    # Age-structured transit cohorts: cohorts[a] = mass that entered the
    # subsurface a years ago and has not yet left (post-attenuation).
    cohorts = np.zeros(par.max_age + 1)
    gw_sorbed = 0.0     # cumulative sorption sink along pathways (a store)
    network = 0.0       # cumulative network-retained surface P (a store)
    deficit = 0.0       # unmet negative-surplus drawdown (a signed store)

    out = {k: np.zeros(n) for k in
           ("Mp", "Ma", "Mdis", "Msorb_r", "Msorb_i", "Mmin", "Msource",
            "transit", "gw_sorbed", "network", "deficit")}
    load_sub = np.zeros(n); load_surf = np.zeros(n)
    leach = np.zeros(n); eros = np.zeros(n); miner = np.zeros(n)
    surplus_applied = np.zeros(n); resid = np.zeros(n)

    def total_stores():
        return (Mp + Ma + Mdis + Msorb_r + Msorb_i
                + cohorts.sum() + gw_sorbed + network - deficit)

    prev_total = total_stores()
    cum_ww = 0.0
    for t in range(n):
        Np = surplus_raw[t]
        if par.negative_surplus == "clamp":
            Np_applied = max(Np, 0.0)
        elif par.negative_surplus == "drawdown":
            Np_applied = Np
        else:
            raise ValueError(f"unknown negative_surplus switch: {par.negative_surplus}")
        surplus_applied[t] = Np_applied

        # --- VM21 Eqs. 7-8: surplus enters the organic pools, split by h.
        if Np_applied >= 0.0:
            Mp += par.h * Np_applied
            Ma += (1.0 - par.h) * Np_applied
        else:
            # Drawdown branch: mine mineral first, then active, then
            # protected; any remainder goes to the deficit ledger.
            need = -Np_applied
            for name in ("mineral", "Ma", "Mp"):
                if need <= 0.0:
                    break
                if name == "mineral":
                    avail = Mdis + Msorb_r + Msorb_i
                    take = min(avail, need)
                    if avail > 0.0:
                        frac = take / avail
                        Mdis -= frac * Mdis; Msorb_r -= frac * Msorb_r
                        Msorb_i -= frac * Msorb_i
                    need -= take
                elif name == "Ma":
                    take = min(Ma, need); Ma -= take; need -= take
                else:
                    take = min(Mp, need); Mp -= take; need -= take
            deficit += need

        # --- VM21 Eq. 10 terms 1-2: first-order mineralisation.
        f_min = par.ka * Ma + par.kp * Mp
        Ma -= par.ka * Ma
        Mp -= par.kp * Mp
        miner[t] = f_min
        # Mineralised P joins the mobile mineral pool; the equilibrium
        # split is re-solved below (VM21 sec. 2.1.4).
        mobile = Mdis + Msorb_r + f_min

        # --- VM21 Eq. 16: gross erosion of organic + sorbed P.
        if par.surface_pathway:
            e_org = par.keros * (Mp + Ma)
            e_sorb = par.keros * (Msorb_r + Msorb_i)
            Mp -= par.keros * Mp; Ma -= par.keros * Ma
            # Eroded sorbed P is taken from the sorbed pools before the
            # equilibrium re-split; represent it against mobile+irrev.
            e_irrev = par.keros * Msorb_i
            e_rev = par.keros * Msorb_r
            Msorb_i -= e_irrev
            mobile -= e_rev
            eros_gross = e_org + e_sorb
        else:
            eros_gross = 0.0
        eros[t] = eros_gross
        # --- VM21 Eq. 15: network retention of the surface flux.
        load_surf[t] = (1.0 - par.Kret) * eros_gross
        network += par.Kret * eros_gross

        # --- VM21 Eqs. 12-13: equilibrium split of the mobile pool.
        Mdis, Msorb_r = freundlich_split(max(mobile, 0.0), par)

        # --- VM21 Eq. 11: leaching out of the source zone.
        flush = min(Q[t] / par.Vw, 1.0)   # Q*dt < Vw branch / capped branch
        Js = Mdis * flush
        Mdis -= Js
        leach[t] = Js

        # --- VM21 sec. 2.1.4: irreversible precipitation of a fraction
        # of the reversibly sorbed pool.
        precip = par.ksorb * Msorb_r
        Msorb_r -= precip
        Msorb_i += precip

        # --- VM21 Eq. 1 / VM17 Eq. S2.1: subsurface transit.
        # New leached mass enters as age-0 cohort; each cohort is
        # attenuated by exp(-gamma) per year in transit (the discrete
        # form of the exp(-gamma*tau) factor) and releases the share
        # implied by the TTD hazard w(a) / (1 - CDF(a)).
        cohorts[1:] = cohorts[:-1]
        cohorts[0] = Js
        cdf = np.cumsum(w)
        outflow = 0.0
        for a in range(par.max_age + 1):
            m = cohorts[a]
            if m == 0.0:
                continue
            decayed = m * (1.0 - np.exp(-par.gamma))
            m -= decayed
            gw_sorbed += decayed
            surv = 1.0 - (cdf[a - 1] if a > 0 else 0.0)
            hazard = (w[a] / surv) if surv > 0 else 1.0
            if a == par.max_age:
                hazard = 1.0
            rel = m * min(hazard, 1.0)
            cohorts[a] = m - rel
            outflow += rel
        load_sub[t] = outflow

        cum_ww += ww[t]

        # --- Book-keeping: the year's closure identity is
        #   Np_applied + ww_in_stores(0) == dStores + outlet exports
        # with wastewater passing straight to the outlet (VM21 Eq. 1).
        exports = load_sub[t] + load_surf[t]
        now_total = total_stores()
        resid[t] = (Np_applied - exports) - (now_total - prev_total)
        prev_total = now_total

        out["Mp"][t] = Mp; out["Ma"][t] = Ma
        out["Mdis"][t] = Mdis; out["Msorb_r"][t] = Msorb_r
        out["Msorb_i"][t] = Msorb_i
        out["Mmin"][t] = Mdis + Msorb_r + Msorb_i
        out["Msource"][t] = Mp + Ma + Mdis + Msorb_r + Msorb_i
        out["transit"][t] = cohorts.sum()
        out["gw_sorbed"][t] = gw_sorbed
        out["network"][t] = network
        out["deficit"][t] = deficit

    load_total = load_sub + load_surf + ww
    return ElementPResult(
        years=years, surplus=surplus_applied, surplus_raw=surplus_raw,
        load_subsurface=load_sub, load_surface=load_surf,
        load_total=load_total, leaching=leach, erosion_gross=eros,
        mineralisation=miner, stores=out, balance_residual=resid,
        deficit=out["deficit"])


# ----------------------------------------------------------------------
# Surplus-trajectory loaders for the held project data (machinery tests)
# ----------------------------------------------------------------------
def surplus_from_components(years, inputs, removal):
    """Np = inputs - removal (VM21 Eqs. 4-6 surface-balance definition)."""
    return np.asarray(inputs, float) - np.asarray(removal, float)
