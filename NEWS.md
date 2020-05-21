# EpiEstim 2.2-3

## FIXED
* Fixed bugs in draw_one_set_of_ancestries resulting from incorrect lengths and an undefined variable (issue #92) (#93, @jstockwin)
* Fixed incorrect quantiles (issue #88) (#89, @jstockwin)

# EpiEstim 2.2-2

## MISC

* Plotting no longer displays TableGrob output (#87, @zkamvar).

# EpiEstim 2.2-1

This release contains various spelling fixes for CRAN maintenance.

# EpiEstim 2.2-0

## NEW FUNCTIONS

* `sample_posterior_R()` samples values of R from the posterior distribution of
  an `estimate_R` object (#70, @acori)

## MISC

* Added a `NEWS.md` file to track changes to the package. (#74, @zkamvar)
* Added tests for plotting with vdiffr. (#74, @zkamvar)
* Remove un-used dependencies that were added during hackout3: plyr, grid, 
  and plotly (#74, @zkamvar)
* Remove compare package from suggests and use test_that version. (#74, @zkamvar)
* Bump minimum required version of coarseDataTools to 0.6-4 (#71, @zkamvar)
* Incidence objects are now handled appropriately with accessors (#65, @zkamvar)

# EpiEstim 2.1-0

* Changed function names to snake_case (only exception is that R remains
  capital letter to avoid confusion between the reproduction number R and the
  growth rate r) and to be more explicit; so `EstimateR` becomes `estimate_R`,
  `OverallInfectivity` becomes `oberall_infectivity`, `WT` becomes
  `wallinga_teunis`, and `DiscrSI` becomes `discr_si`. Names of arguments to
  these functions have also changed to snake_case. Note that compatibility
  functions have been added so that the old functions as written in EpiEstim
  1.1-0 should still work but throw a warning pointing to the newest functions. 
* Compatibility with `incidence` package: in the function `estimate_R`, the
  first argument, i.e. the incidence from which the reproduction number is
  calculated, can now be, either a vector of case counts (as in version 1.1-0) or
  an `incidence` object (see R package `incidence`).
* Accounting for imported cases: in the function `estimate_R`, the first
  argument, i.e. the incidence from which the reproduction number can now
  provide information about known imported cases: by specifying the first
  argument as either a dataframe with columns "local" and "imported", or an
  `incidence` object with two groups (local and imported, see R package
  `incidence`). This new feature is described in Thompson et al. Epidemics 2019
  (currently in review).
* Additional methods available for function `estimate_R`: in addition to
  `non_parametric_si`, `parametric_si` and `uncertain_si`, which were already
  available in EpiEstim 1.1-0, two new methods have been added: `si_from_data` or
  `si_from_sample`. These allow feeding function `estimate_R` data on observed
  serial intervals (method `si_from_data`) or posterior samples of serial
  interval distributions obtained from such data (method `si_from_sample`). These
  new features are described in Thompson et al. Epidemics 2019 (currently in
  review).
* No more plotting option inside of `estimate_R`: estimate_R now generates on
  object of class `estimate_R`, which can be plotted separately by using the
  new `estimate_R_plots` function, which also now allows to plot several R
  estimates on a single plot. 
* New argument `config` for `estimate_R` function: this is meant to minimise
  the number of arguments to function `estimate_R`; so arguments `method`,
  `t_start`, `t_end`, `n1`, `n2`, `mean_si`, `std_si`, `std_mean_si`,
  `min_mean_si`, `max_mean_si`, `std_std_si`, `min_std_si`, `max_std_si`,
  `si_distr`, `mean_prior`, `std_prior`, and `cv_posterior` are now specified as
  a group under this new `config` argument. Such a `config` argument must be of
  class `estimate_R_config` and can be obtained as a results of the new
  `make_config` function. 
* New function `make_config`, which defines settings for function `estimate_R`,
  and sets defaults where arguments are missing. In particular, if argument
  `incid` is not `NULL`, by default `config$t_start` and `config$t_end` will be
  set so that, when the configuration is used inside `estimate_R` function, the
  reproduction number is estimated by default on sliding weekly windows (in
  EpiEstim 1.1-0 there was no default for the time window of estimation of R).
* Added a vignette to illustrate main features of the package.

## NEW DATASETS: 

 - `flu_2009_NYC_school`
 - `mers_2014_15`, 
 - `MockRotavirus`

## NEW IMPORTS

 - `stats` (to use the gamma distribution; it was already used in EpiEstim 1.1-0
   but making the dependency explicit)
 - `coarseDataTools`, `fitdistrplus`, `coda` (used for the new methods
   `si_from_data` and `si_from_sample` in `estimate_R` function to estimate the
   serial interval from data). 
 - `incidence` (so that `estimate_R` can take an `incidence` object as first
   argument)
 - `graphics`, `reshape2`, `ggplot2`, `gridExtra`, `scales`, `grDevices` (to
   make new plots of outputs of `estimate_R` and `wallinga_teunis` functions)
