// Bayesian SPARROW for the Ontario (Lake Ontario tributary) ensemble.
//
// This is the Bayesian member of the ensemble described in the Ontario
// SPARROW proposal package.  It implements the model design settled on
// 2026-11-20, on the data structures of the classical RSPARROW port in
// this folder (classical_sparrow.py; see PORT_NOTES.md, "Mapping onto
// the Stan Bayesian version").
//
// DELIBERATE SEMANTIC DIFFERENCE FROM THE CLASSICAL MEMBER
// --------------------------------------------------------
// The accumulation here is UNCONDITIONED.  RSPARROW's calibration loop
// (tnoder.for, ifadjust = 1) substitutes the observed load for the
// predicted load at monitoring reaches before accumulating downstream.
// That is a nonlinear-least-squares estimation device, not a generative
// step, so it has no place inside a Bayesian model.  Here every reach
// load is generated from the parameters alone, and the observed station
// loads enter only through the likelihood at their reaches.
//
// NETWORK REPRESENTATION
// ----------------------
// Reaches are supplied in topological (upstream-to-downstream) order.
// For each reach, `downstream1` / `downstream2` give the positions (in
// the same ordering) of the reach or reaches leaving this reach's
// to-node; the network is dendritic except for a single flow-split
// node, hence at most two.  Position N + 1 is a dummy sink that
// collects load leaving the network at terminal reaches.  The loop
// maintains `inflow[r]`, the total load delivered to reach r's
// from-node by all reaches upstream of it -- the exact analogue of the
// Fortran `node(fnode)` accumulator, pushed downstream reach by reach:
//
//   total[r] = inc[r] + carry[r] * inflow[r]
//   inflow[downstream*[r]] += total[r]
//
// where inc[r]   = (in-reach attenuation, half-reach) * G[r]
//       carry[r] = frac[r] * (full-reach attenuation)
// and G[r] is the delivered source load generated inside reach r.
// This matches the classical port's incremental_loads/_accumulate
// (ptnoder, ifadjust = 0) line for line.
//
// SOURCES (per reach, loads in kg P / yr, areas in km^2)
// ------------------------------------------------------
//   point-source P      : coefficient fixed at 1 (data: src_point)
//   fertilizer P        : x alpha_fert  (estimated delivery coefficient)
//   manure P            : x alpha_man   (estimated delivery coefficient)
//   urban land area     : x rate_urban  (estimated export rate, kg/km^2/yr)
//   background land     : x rate_bg     (forest + grass/shrub + wetland
//                                        + pasture area; estimated rate)
//   atmospheric deposition on water area : x rate_dep (estimated rate,
//                                        informative prior, see below)
// Fertilizer and manure additionally share a land-to-stream delivery
// factor exp(beta_runoff * z_runoff + beta_perm * z_perm
//            + beta_tile * z_tile),
// with the three modifiers standardised over reaches (data).  The
// design expects beta_runoff > 0, beta_perm < 0, beta_tile > 0; the
// priors are centred at zero and the signs are left to the data.
//
// ATTENUATION
// -----------
//   stream reaches    : A = exp(-ks * ttime), ttime = travel time (days)
//                       from reach length / mean velocity
//   reservoir reaches : B = 1 / (1 + beta_res * rhload), rhload = areal
//                       hydraulic load from the source data1 table,
//                       supplied scaled by its median over reservoir
//                       reaches (see prepare_stan_data.py), so a value
//                       of 1 is a typical reservoir reach
// In-reach generated load is attenuated by sqrt(A) on stream reaches
// (half-reach convention, as in RSPARROW) and by B on reservoir
// reaches; load entering from upstream is attenuated by A (resp. B)
// in full, times the reach's flow fraction `frac`.
//
// OBSERVATION MODEL
// -----------------
// Each observation is one station's mean annual load (kg/yr) over its
// record.  In log space:
//   log_obs ~ normal(log(mu[obs_reach]), sqrt(obs_se^2 + sigma^2))
// where obs_se (data, per station) combines the load-estimator
// standard error of the mean with the interannual variance of the
// station's annual loads (both computed in prepare_stan_data.py), and
// sigma is a fitted residual standard deviation.
//
// PRIORS (weakly informative lognormal unless noted; medians in brackets)
//   alpha_fert, alpha_man ~ lognormal(log 0.3, 0.5)     [0.3]
//   rate_urban            ~ lognormal(log 300, 0.5)     [300 kg/km^2/yr]
//   rate_bg               ~ lognormal(log 25, 0.5)      [25]
//   rate_dep              ~ lognormal(log 500, 0.25)    [500; informative,
//                            centred on 50 mg/m^2/yr, the Eimers et al.
//                            2022 deposition value used in the design]
//   ks                    ~ lognormal(log 0.15, 0.5)    [0.15 /day]
//   beta_res              ~ lognormal(log 0.5, 0.75)    [0.5 at the median
//                            scaled rhload, i.e. retention 1/(1+0.5) = 0.67
//                            on a typical reservoir reach]
//   beta_runoff/perm/tile ~ normal(0, 0.5)
//   sigma                 ~ normal(0, 0.5), truncated below at 0

data {
  int<lower=1> N;                       // reaches, topological order
  // sources and land areas per reach
  vector<lower=0>[N] src_point;         // kg P/yr, coefficient fixed at 1
  vector<lower=0>[N] src_fert;          // kg P/yr
  vector<lower=0>[N] src_man;           // kg P/yr
  vector<lower=0>[N] area_urban;        // km^2
  vector<lower=0>[N] area_bg;           // km^2 (forest+grass+wetland+pasture)
  vector<lower=0>[N] area_water;        // km^2
  // delivery modifiers, standardised over reaches
  vector[N] z_runoff;
  vector[N] z_perm;
  vector[N] z_tile;
  // attenuation inputs
  vector<lower=0>[N] ttime;             // days (stream reaches)
  array[N] int<lower=0, upper=1> is_res;  // 1 = reservoir reach (rchtype 1)
  vector<lower=0>[N] rhload;            // scaled areal hydraulic load (reservoir reaches, else 0)
  vector<lower=0, upper=1>[N] frac;     // flow fraction through the reach
  // topology: positions in 1..N, N + 1 = dummy sink
  array[N] int<lower=1, upper=N + 1> downstream1;
  array[N] int<lower=1, upper=N + 1> downstream2;
  // observations: station mean annual loads
  int<lower=0> n_obs;
  array[n_obs] int<lower=1, upper=N> obs_reach;
  vector[n_obs] log_obs;                // log(mean annual load, kg/yr)
  vector<lower=0>[n_obs] obs_se;        // log-space SE per station
  // generated-quantities bookkeeping: basin outlet reaches
  int<lower=0> n_outlet;
  array[n_outlet] int<lower=1, upper=N> outlet_idx;
}

parameters {
  real<lower=0> alpha_fert;
  real<lower=0> alpha_man;
  real<lower=0> rate_urban;
  real<lower=0> rate_bg;
  real<lower=0> rate_dep;
  real beta_runoff;
  real beta_perm;
  real beta_tile;
  real<lower=0> ks;
  real<lower=0> beta_res;
  real<lower=0> sigma;
}

transformed parameters {
  vector[N] mu;                         // predicted total load per reach, kg/yr
  {
    vector[N + 1] inflow = rep_vector(0, N + 1);
    real deliv;
    real G;
    real A;
    real B;
    real inc_factor;
    real tot;
    for (r in 1:N) {
      deliv = exp(beta_runoff * z_runoff[r] + beta_perm * z_perm[r]
                  + beta_tile * z_tile[r]);
      G = src_point[r]
          + deliv * (alpha_fert * src_fert[r] + alpha_man * src_man[r])
          + rate_urban * area_urban[r]
          + rate_bg * area_bg[r]
          + rate_dep * area_water[r];
      if (is_res[r] == 1) {
        B = 1.0 / (1.0 + beta_res * rhload[r]);
        inc_factor = B;
        tot = inc_factor * G + frac[r] * B * inflow[r];
      } else {
        A = exp(-ks * ttime[r]);
        inc_factor = sqrt(A);
        tot = inc_factor * G + frac[r] * A * inflow[r];
      }
      mu[r] = tot;
      inflow[downstream1[r]] += tot;
      if (downstream2[r] != downstream1[r]) {
        inflow[downstream2[r]] += tot;
      }
    }
  }
}

model {
  // priors
  alpha_fert ~ lognormal(log(0.3), 0.5);
  alpha_man ~ lognormal(log(0.3), 0.5);
  rate_urban ~ lognormal(log(300), 0.5);
  rate_bg ~ lognormal(log(25), 0.5);
  rate_dep ~ lognormal(log(500), 0.25);
  ks ~ lognormal(log(0.15), 0.5);
  beta_res ~ lognormal(log(0.5), 0.75);
  beta_runoff ~ normal(0, 0.5);
  beta_perm ~ normal(0, 0.5);
  beta_tile ~ normal(0, 0.5);
  sigma ~ normal(0, 0.5);
  // likelihood at station reaches (unconditioned accumulation)
  for (o in 1:n_obs) {
    log_obs[o] ~ normal(log(mu[obs_reach[o]]),
                        sqrt(square(obs_se[o]) + square(sigma)));
  }
}

generated quantities {
  vector[n_outlet] outlet_load;         // predicted load at basin outlets, kg/yr
  vector[n_obs] mu_site;                // predicted load at station reaches
  vector[n_obs] log_lik;                // pointwise log likelihood
  for (j in 1:n_outlet) {
    outlet_load[j] = mu[outlet_idx[j]];
  }
  for (o in 1:n_obs) {
    mu_site[o] = mu[obs_reach[o]];
    log_lik[o] = normal_lpdf(log_obs[o] | log(mu[obs_reach[o]]),
                             sqrt(square(obs_se[o]) + square(sigma)));
  }
}
