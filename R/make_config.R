#' Set and check parameter settings for [estimate_R()]
#'
#' Define parameters settings for [estimate_R()]. It takes a list of named items
#' as input, sets defaults where arguments are missing, and returns a list of
#' settings.
#'
#' @param ... Acceptable arguments for `...` are:
#'
#' - `t_start`: Vector of positive integers giving the starting times of each
#'   window over which the reproduction number will be estimated. These must be
#'   in ascending order, and so that for all `i`, `t_start[i] <= t_end[i]`.
#'   `t_start[1]` should be strictly after the first day with non null
#'   incidence.
#'
#' - `t_end`: Vector of positive integers giving the ending times of each window
#'   over which the reproduction number will be estimated. These must be in
#'   ascending order, and so that for all `i`, `t_start[i] <= t_end[i]`.
#'
#' - `n1`: For method "uncertain_si" and "si_from_data"; positive integer giving
#'   the size of the sample of SI distributions to be drawn (see details).
#'
#' - `n2`: For methods "uncertain_si", "si_from_data" and "si_from_sample";
#'   positive integer giving the size of the sample drawn from the posterior
#'   distribution of R for each serial interval distribution considered (see
#'   details).
#'
#' - `mean_si`: For method "parametric_si" and "uncertain_si"; positive real
#'   giving the mean serial interval (method "parametric_si") or the average
#'   mean serial interval (method "uncertain_si", see details).
#'
#' - `std_si`: For method "parametric_si" and "uncertain_si"; non negative real
#'   giving the standard deviation of the serial interval (method
#'   "parametric_si") or the average standard deviation of the serial interval
#'   (method "uncertain_si", see details).
#'
#' - `std_mean_si`: For method "uncertain_si"; standard deviation of the
#'   distribution from which mean serial intervals are drawn (see details).
#'
#' - `min_mean_si`: For method "uncertain_si"; lower bound of the distribution
#'   from which mean serial intervals are drawn (see details).
#'
#' - `max_mean_si`: For method "uncertain_si"; upper bound of the distribution
#'   from which mean serial intervals are drawn (see details).
#'
#' - `std_std_si`: For method "uncertain_si"; standard deviation of the
#'   distribution from which standard deviations of the serial interval are
#'   drawn (see details).
#'
#' - `min_std_si`: For method "uncertain_si"; lower bound of the distribution
#'   from which standard deviations of the serial interval are drawn (see
#'   details).
#'
#' - `max_std_si`: For method "uncertain_si"; upper bound of the distribution
#'   from which standard deviations of the serial interval are drawn (see
#'   details).
#'
#' - `si_distr`: For method "non_parametric_si"; vector of probabilities giving
#'   the discrete distribution of the serial interval, starting with
#'   `si_distr[1]` (probability that the serial interval is zero), which should
#'   be zero. Note that EpiEstim assumes that the serial interval is always
#'   strictly positive. 
#'
#' - `si_parametric_distr`: For method "si_from_data"; the parametric
#'   distribution to use when estimating the serial interval from data on dates
#'   of symptoms of pairs of infector/infected individuals (see details). Should
#'   be one of "gamma", "lognormal", "weibull" or "exponential", or one of these
#'   followed by "_offset_1" (e.g. "gamma_offset_1") for a serial interval
#'   shifted by 1. Other distributions can be used by estimating the serial
#'   interval with primarycensored and passing it to [primary2estim()] and
#'   method "si_from_sample".
#'
#' - `si_discr_args`: For methods "parametric_si", "uncertain_si" and
#'   "si_from_data"; a named
#'   list of additional arguments passed to [discr_si()] when discretising the
#'   serial interval, e.g. `list(dist = stats::plnorm)`. Can contain `dist`
#'   ([stats::pgamma()] or [stats::plnorm()]), `shift`, `L`, `D`, `dprimary`
#'   and `primary_args`. The resulting distribution must give zero
#'   probability to a serial interval of zero. Defaults to an empty list, which
#'   uses the defaults of [discr_si()]. For method "si_from_data", `dprimary`
#'   and `primary_args` are also used when estimating the serial interval,
#'   and `dist` and `shift` are set by `si_parametric_distr`.
#'
#' - `mcmc_control`: Deprecated. For method "si_from_data"; an object of class
#'   \code{estimate_R_mcmc_control}, as returned by function
#'   \code{make_mcmc_control}, giving the seed used to draw the sample of
#'   serial interval distributions and starting values for their estimation.
#'   Defaults to `NULL`, in which case `seed` is used and starting values are
#'   given by [si_start_values()].
#'
#' - `mean_prior`: A positive number giving the mean of the common prior
#'   distribution for all reproduction numbers (see details).
#'
#' - `std_prior`: A positive number giving the standard deviation of the
#'   common prior distribution for all reproduction numbers (see details).
#'
#' - `cv_posterior`: A positive number giving the aimed posterior coefficient
#'   of variation (see details).
#' 
#' @inheritParams estimate_R
#'
#' @details
#' Analytical estimates of the reproduction number for an epidemic over
#' predefined time windows can be obtained using function [estimate_R()], for a
#' given discrete distribution of the serial interval. `make_config()` generates 
#' configuration specifying the way the estimation will be performed.
#'
#' The more incident cases observed over a time window, the smallest the
#' posterior coefficient of variation (CV, ratio of standard deviation over
#' mean) of the reproduction number. An aimed CV can be specified in the
#' argument `cv_posterior` (default is `0.3`), and a warning will be produced if
#' the incidence within one of the time windows considered is too low to get
#' this CV.
#'
#' ## Methods
#' 
#' The methods vary in the way the serial interval distribution is specified.
#'
#' In short there are five methods to specify the serial interval distribution
#' (see below for details on each method). This is specified in the argument
#' `method` of the [estimate_R()] function. In the first two methods, a unique
#' serial interval distribution is considered, whereas in the last three, a
#' range of serial interval distributions are integrated over:
#' - "non_parametric_si": the user specifies the discrete distribution
#'   of the serial interval
#' - "parametric_si": the user specifies the mean and sd of the serial
#'   interval
#' - "uncertain_si": the mean and sd of the serial interval are each
#'   drawn from truncated normal distributions, with parameters specified by the
#'   user
#' - "si_from_data": the serial interval distribution is directly
#'   estimated, by maximum likelihood, from interval censored exposure data, with data
#'   provided by the user together with a choice of parametric distribution for
#'   the serial interval
#' - "si_from_sample": the user directly provides the sample of serial
#'   interval distribution to use for estimation of R. This can be a useful
#'   alternative to the previous method, where the estimation of the serial
#'   interval distribution could be run once, and the same estimated SI
#'   distribution then used in [estimate_R()] in different contexts, e.g. with
#'   different time windows, hence avoiding having to rerun the estimation every time
#'   [estimate_R()] is called.
#'
#' ### `method = "non_parametric_si"`
#'
#' The discrete distribution of the serial interval is directly specified in the
#' argument `si_distr`.
#'
#' ### `method = "parametric_si"`
#'
#' The mean and standard deviation of the continuous distribution of the serial
#' interval are given in the arguments `mean_si` and `std_si`. The discrete
#' distribution of the serial interval is derived automatically using
#' [discr_si()].
#'
#' ### `method = "uncertain_si"`
#'
#' Method "uncertain_si" allows accounting for uncertainty on the serial
#' interval distribution as described in Cori et al. AJE 2013. We allow the mean
#' \eqn{\mu} and standard deviation \eqn{\sigma} of the serial interval to vary
#' according to truncated normal distributions. We sample `n1` pairs of mean and
#' standard deviations,
#' \eqn{(\mu^{(1)},\sigma^{(1)}),...,(\mu^{(n_1)},\sigma^{(n_1)})}, by 
#' independently
#' sampling the mean \eqn{\mu^{(k)}} from its truncated normal distribution
#' (with mean `mean_si`, standard deviation `std_mean_si`, minimum `min_mean_si`
#' and maximum `max_mean_si`), and the standard deviation
#' \eqn{\sigma^{(k)}} from its truncated normal distribution (with mean
#' `std_si`, standard deviation `std_std_si`, minimum `min_std_si` and maximum
#' `max_std_si`). Warnings are produced when the truncated
#' normal distributions are not symmetric around the mean. For each pair
#' \eqn{(\mu^{(k)},\sigma^{(k)})}, we then draw a sample of size `n2` in the
#' posterior distribution of the reproduction number over each time window,
#' conditionally on this parametric serial interval distribution (using the 
#' `discr_si` function to generate the corresponding serial interval probability
#' mass function). 
#' After pooling across the `n1` serial interval distributions, a sample
#' of size \eqn{`n1` \times `n2`} of the joint posterior distribution of the
#' reproduction number over each time window is obtained. The posterior mean,
#' standard deviation, and 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975 quantiles
#' of the reproduction number for each time window are obtained from this
#' sample.
#'
#'### `method = "si_from_data"`
#'
#' Method "si_from_data" allows accounting for uncertainty on the serial
#' interval distribution. Unlike method "uncertain_si", where we arbitrarily
#' vary the mean and std of the SI in truncated normal distributions, here, the
#' scope of serial interval distributions considered is directly informed by
#' data on the (potentially censored) dates of symptoms of pairs of
#' infector/infected individuals. This data, specified in argument `si_data`,
#' should be a dataframe with 4 to 6 columns:
#' - `EL`: the lower bound of the symptom onset date of the infector (given as
#'   an integer)
#' - `ER`: the upper bound of the symptom onset date of the infector (given as
#'   an integer). Should be such that `ER >= EL`. If the dates are known exactly
#'   use `ER = EL`
#' - `SL`: the lower bound of the symptom onset date of the infected individual
#'   (given as an integer)
#' - `SR`: the upper bound of the symptom onset date of the infected individual
#'   (given as an integer). Should be such that `SR >= SL`. If the dates are
#'   known exactly use `SR = SL`
#' - `type` (optional): can have entries 0, 1, or 2, corresponding to doubly
#'   interval-censored, single interval-censored or exact observations,
#'   respectively, see Reich et al. Statist. Med. 2009. If not specified, this
#'   will be automatically computed from the dates
#'
#' - `OT` (optional): the time (given as an integer) up to which symptom
#'   onsets of infected individuals are observed. Like `SR` it is a continuous
#'   bound, so with daily data a pair observed up to and including day `d` has
#'   `OT = d + 1`. Should be such that `OT > SL`. If given, the estimation
#'   accounts for right truncation. When pairs of infector/infected
#'   individuals are observed during an ongoing outbreak, `OT` should be
#'   given, as otherwise the serial interval may be underestimated (see
#'   Charniga et al. PLoS Comp Biol 2024). If not given, or for entries that
#'   are `NA`, no right truncation is assumed.
#'
#' As in coarseDataTools, `EL`, `ER`, `SL` and `SR` are continuous bounds, so a
#' symptom onset known to the day `d` is given as `EL = d` and `ER = d + 1` (as
#' in the `MockRotavirus` data), and `ER = EL` means the onset time is known
#' exactly.
#'
#' Assuming a given parametric distribution for the serial interval distribution
#' (specified in `si_parametric_distr`), the serial interval is estimated
#' directly from these data by maximum likelihood, accounting for double
#' interval censoring, using [primarycensored::fitdistdoublecens()] (Abbott et
#' al., \doi{10.5281/zenodo.13632839}). A sample of `n1` serial interval
#' distributions is then drawn from the asymptotic normal distribution of the
#' parameter estimates, using `seed`. For a Bayesian estimate of the serial
#' interval, fit it with [primarycensored::pcd_cmdstan_model()] and use
#' [primary2estim()] with method "si_from_sample".
#' For each element in the sample of serial interval distributions, we
#' then draw a sample of size `n2` in the posterior distribution of the
#' reproduction number over each time window, conditionally on this serial
#' interval distribution. After pooling, a sample of size \eqn{`n1` \times `n2`}
#' of the joint posterior distribution of the reproduction number over each time
#' window is obtained. The posterior mean, standard deviation, and 0.025, 0.05,
#' 0.25, 0.5, 0.75, 0.95, 0.975 quantiles of the reproduction number for each
#' time window are obtained from this sample.
#'
#'### `method = "si_from_sample"`
#' Method "si_from_sample" also allows accounting for uncertainty on the serial
#' interval distribution. Unlike methods "uncertain_si" and "si_from_data", the
#' user directly provides (in argument `si_sample`) a sample of serial interval
#' distribution to be explored.
#'
#' @return An object of class `estimate_R_config` with components
#' `t_start`, `t_end`, `n1`, `n2`, `mean_si`, `std_si`,
#' `std_mean_si`, `min_mean_si`, `max_mean_si`, `std_std_si`, `min_std_si`, `max_std_si`,
#' `si_distr`, `si_parametric_distr`, `si_discr_args`, `mcmc_control`, `seed`, `mean_prior`, `std_prior`,
#' `cv_posterior`, which can be used as an argument of function [estimate_R()].
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ## load data on rotavirus
#' data("MockRotavirus")
#'
#' ## estimate the reproduction number (method "si_from_data")
#' ## we are not specifying the time windows, so by defaults this will estimate
#' ## R on sliding weekly windows
#' incid <- MockRotavirus$incidence
#' method <- "si_from_data"
#' config <- make_config(incid = incid,
#'                      list(si_parametric_distr = "gamma",
#'                      n1 = 500,
#'                      n2 = 50,
#'                      seed = 2))
#'
#' R_si_from_data <- estimate_R(incid,
#'                             method = method,
#'                             si_data = MockRotavirus$si_data,
#'                             config = config)
#' plot(R_si_from_data)                     
#'
#' ## you can also create the config straight within the estimate_R call,
#' ## in that case incid and method are automatically used from the estimate_R
#' ## arguments:
#' R_si_from_data <- estimate_R(incid,
#'                             method = method,
#'                             si_data = MockRotavirus$si_data,
#'                             config = make_config(
#'                      list(si_parametric_distr = "gamma",
#'                      n1 = 500,
#'                      n2 = 50,
#'                      seed = 2)))
#' plot(R_si_from_data)
#' }

make_config <- function(..., incid = NULL) {
  config <- list(...)
  if (length(config) == 1L && is.list(config[[1]])) {
    config <- config[[1]]
  }
  
  # catch if user (wrongly) specifies method
  if(!is.null(config$method)) {
    msg <- paste("`method` should be specified as an argument to",
                 "`estimate_R`, not `make_config`.")
    stop(msg)
  }

  ## SET DEFAULTS
  defaults <- list(t_start = NULL,
                   t_end = NULL,
                   n1 = 500,
                   n2 = 50,
                   mean_si = NULL,
                   std_si = NULL,
                   std_mean_si = NULL,
                   min_mean_si = NULL,
                   max_mean_si = NULL,
                   std_std_si = NULL,
                   min_std_si = NULL,
                   max_std_si = NULL,
                   si_distr = NULL,
                   si_parametric_distr = NULL,
                   si_discr_args = list(),
                   mcmc_control = NULL,
                   seed = NULL,
                   mean_prior = 5,
                   std_prior = 5,
                   cv_posterior = 0.3)

  ## MODIFY CONFIG WITH ARGUMENTS ##
  config <- modify_defaults(defaults, config)

  ## checking and processing incid
  if (!is.null(incid)) {
    incid <- process_I(incid)
    idx_raw_incid <- as.integer(rownames(incid)) > 0
    T <- sum(idx_raw_incid) # nolint: object_overwrite_linter.

    ## filling in / checking t_start and t_end
    if (is.null(config$t_start) || is.null(config$t_end)) {
      msg <- "Default config will estimate R on weekly sliding windows.
    To change this change the t_start and t_end arguments. "
      message(msg)
      config$t_start <- seq(2, T-6)
      config$t_end <- seq(8, T)
    } else {
      check_times(config$t_start, config$t_end, T)
    }
  }

  class(config) <- "estimate_R_config"
  return(config)

}

