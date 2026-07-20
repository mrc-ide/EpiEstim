# Package index

## Pre-processing

- [`aggregate_inc()`](https://mrc-ide.github.io/EpiEstim/reference/aggregate_inc.md)
  : Aggregating daily incidence to longer time windows
- [`backimpute_I()`](https://mrc-ide.github.io/EpiEstim/reference/backimpute_I.md)
  : Impute unobserved generations of infection
- [`compute_si_cutoff()`](https://mrc-ide.github.io/EpiEstim/reference/compute_si_cutoff.md)
  : Index before which at most a given probability mass is captured
- [`compute_t_min()`](https://mrc-ide.github.io/EpiEstim/reference/compute_t_min.md)
  : Compute the smallest index at which joint estimation should start
- [`discr_si()`](https://mrc-ide.github.io/EpiEstim/reference/discr_si.md)
  : Compute discretized generation time distribution
- [`first_nonzero_incid()`](https://mrc-ide.github.io/EpiEstim/reference/first_nonzero_incid.md)
  : First day of non-zero incidence
- [`get_shape_R_flat()`](https://mrc-ide.github.io/EpiEstim/reference/get_shape_R_flat.md)
  : Precompute shape of posterior distribution for R
- [`get_shape_epsilon()`](https://mrc-ide.github.io/EpiEstim/reference/get_shape_epsilon.md)
  : Precompute shape of posterior distribution for epsilon
- [`si_from_data_valid_distrs()`](https://mrc-ide.github.io/EpiEstim/reference/si_from_data_valid_distrs.md)
  : Distribution names valid when using MCMC to estimate SI from data

## Estimation

- [`estimate_R()`](https://mrc-ide.github.io/EpiEstim/reference/estimate_R.md)
  : Estimate the instantaneous reproduction number
- [`estimate_R_agg()`](https://mrc-ide.github.io/EpiEstim/reference/estimate_R_agg.md)
  : Estimate instantaneous reproduction number from coarsely aggregated
  data
- [`estimate_advantage()`](https://mrc-ide.github.io/EpiEstim/reference/estimate_advantage.md)
  : Estimate instantaneous reproduction number
- [`compute_lambda()`](https://mrc-ide.github.io/EpiEstim/reference/compute_lambda.md)
  : Compute the overall infectivity
- [`overall_infectivity()`](https://mrc-ide.github.io/EpiEstim/reference/overall_infectivity.md)
  : Overall Infectivity Due To Previously Infected Individuals
- [`process_I_multivariant()`](https://mrc-ide.github.io/EpiEstim/reference/process_I_multivariant.md)
  : Process incidence input for multivariant analyses
- [`wallinga_teunis()`](https://mrc-ide.github.io/EpiEstim/reference/wallinga_teunis.md)
  : Estimate case reproduction number using the Wallinga and Teunis
  method

## Post-processing

### MCMC diagnostics

- [`check_cdt_samples_convergence()`](https://mrc-ide.github.io/EpiEstim/reference/check_cdt_samples_convergence.md)
  : Check MCMC chain convergence using the Gelman-Rubin algorithm

### Draw from posterior

- [`draw_R()`](https://mrc-ide.github.io/EpiEstim/reference/draw_R.md) :
  Draw R from marginal posterior distribution
- [`draw_epsilon()`](https://mrc-ide.github.io/EpiEstim/reference/draw_epsilon.md)
  : Draw epsilon from marginal posterior distribution
- [`sample_posterior_R()`](https://mrc-ide.github.io/EpiEstim/reference/sample_posterior_R.md)
  : Sample from the posterior R distribution

## Helpers

### Default settings

- [`default_mcmc_controls()`](https://mrc-ide.github.io/EpiEstim/reference/default_mcmc_controls.md)
  : Set default for MCMC control

- [`default_priors()`](https://mrc-ide.github.io/EpiEstim/reference/default_priors.md)
  : Set default for Gamma priors

- [`init_mcmc_params()`](https://mrc-ide.github.io/EpiEstim/reference/init_mcmc_params.md)
  : Find clever starting points for MCMC estimation

- [`make_config()`](https://mrc-ide.github.io/EpiEstim/reference/make_config.md)
  :

  Set and check parameter settings for
  [`estimate_R()`](https://mrc-ide.github.io/EpiEstim/reference/estimate_R.md)

- [`make_mcmc_control()`](https://mrc-ide.github.io/EpiEstim/reference/make_mcmc_control.md)
  : Create list of MCMC control parameters

### Plotting functions

- [`plot(`*`<estimate_R>`*`)`](https://mrc-ide.github.io/EpiEstim/reference/plot.estimate_R.md)
  [`plot(`*`<wallinga_teunis>`*`)`](https://mrc-ide.github.io/EpiEstim/reference/plot.estimate_R.md)
  :

  Plot outputs of
  [`estimate_R()`](https://mrc-ide.github.io/EpiEstim/reference/estimate_R.md)

- [`estimate_R_plots()`](https://mrc-ide.github.io/EpiEstim/reference/estimate_R_plots.md)
  : Wrapper for plot.estimate_R

### Backwards-compatibility functions

- [`DiscrSI()`](https://mrc-ide.github.io/EpiEstim/reference/DiscrSI.md)
  : Function to ensure compatibility with EpiEstim versions \<2.0
- [`EstimateR()`](https://mrc-ide.github.io/EpiEstim/reference/EstimateR.md)
  : Function to ensure compatibility with EpiEstim versions \<2.0
- [`OverallInfectivity()`](https://mrc-ide.github.io/EpiEstim/reference/OverallInfectivity.md)
  : Function to ensure compatibility with EpiEstim versions \<2.0
- [`WT()`](https://mrc-ide.github.io/EpiEstim/reference/WT.md) :
  Function to ensure compatibility with EpiEstim versions \<2.0

### External package compatibility functions

- [`coarse2estim()`](https://mrc-ide.github.io/EpiEstim/reference/coarse2estim.md)
  : Link coarseDataTools and EpiEstim

## Datasets

- [`Flu1918`](https://mrc-ide.github.io/EpiEstim/reference/Flu1918.md) :
  Data on the 1918 H1N1 influenza pandemic in Baltimore.
- [`Flu2009`](https://mrc-ide.github.io/EpiEstim/reference/Flu2009.md) :
  Data on the 2009 H1N1 influenza pandemic in a school in Pennsylvania.
- [`Measles1861`](https://mrc-ide.github.io/EpiEstim/reference/Measles1861.md)
  : Data on the 1861 measles epidemic in Hagelloch, Germany.
- [`MockRotavirus`](https://mrc-ide.github.io/EpiEstim/reference/MockRotavirus.md)
  : Mock data on a rotavirus epidemic.
- [`SARS2003`](https://mrc-ide.github.io/EpiEstim/reference/SARS2003.md)
  : Data on the 2003 SARS epidemic in Hong Kong.
- [`Smallpox1972`](https://mrc-ide.github.io/EpiEstim/reference/Smallpox1972.md)
  : Data on the 1972 smallpox epidemic in Kosovo
- [`covid_deaths_2020_uk`](https://mrc-ide.github.io/EpiEstim/reference/covid_deaths_2020_uk.md)
  : Data on the 2020-2022 SARS-CoV-2 epidemic in the UK.
- [`flu_2009_NYC_school`](https://mrc-ide.github.io/EpiEstim/reference/flu_2009_NYC_school.md)
  : Data on the 2009 H1N1 influenza pandemic in a school in New York
  city
- [`mers_2014_15`](https://mrc-ide.github.io/EpiEstim/reference/mers_2014_15.md)
  : Data on Middle East Respiratory Syndrome (MERS) in Saudi Arabia.
