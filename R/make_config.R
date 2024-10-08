#' Set and check parameter settings of estimate_R
#'
#' This function defines settings for estimate_R  It takes a list of named
#' items as input, set defaults where arguments are
#' missing, and return a list of settings.
#'
#' @param ... Acceptable arguments for ... are:
#'
#' \describe{
#'
#' \item{t_start}{Vector of positive integers giving the starting times of each
#' window over which the reproduction number will be estimated. These must be in
#'  ascending order, and so that for all \code{i}, \code{t_start[i]<=t_end[i]}.
#'  t_start[1] should be strictly after the first day with non null incidence.}
#'
#' \item{t_end}{Vector of positive integers giving the ending times of each
#' window over which the reproduction number will be estimated. These must be
#' in ascending order, and so that for all \code{i},
#' \code{t_start[i]<=t_end[i]}.}
#'
#' \item{n1}{For method "uncertain_si" and "si_from_data"; positive integer
#' giving the size of the sample of SI distributions to be drawn (see details).}
#'
#' \item{n2}{For methods "uncertain_si", "si_from_data" and "si_from_sample";
#' positive integer giving the size of the sample drawn from the posterior
#' distribution of R for each serial interval distribution considered (see
#' details).}
#'
#' \item{mean_si}{For method "parametric_si" and "uncertain_si" ; positive real
#' giving the mean serial interval (method "parametric_si") or the average mean
#' serial interval (method "uncertain_si", see details).}
#'
#' \item{std_si}{For method "parametric_si" and "uncertain_si" ; non negative
#' real giving the standard deviation of the serial interval
#' (method "parametric_si") or the average standard deviation of the serial
#' interval (method "uncertain_si", see details).}
#'
#' \item{std_mean_si}{For method "uncertain_si" ; standard deviation of the
#' distribution from which mean serial intervals are drawn (see details).}
#'
#' \item{min_mean_si}{For method "uncertain_si" ; lower bound of the
#' distribution from which mean serial intervals are drawn (see details).}
#'
#' \item{max_mean_si}{For method "uncertain_si" ; upper bound of the
#' distribution from which mean serial intervals are drawn (see details).}
#'
#' \item{std_std_si}{For method "uncertain_si" ; standard deviation of the
#' distribution from which standard deviations of the serial interval are drawn
#' (see details).}
#'
#' \item{min_std_si}{For method "uncertain_si" ; lower bound of the distribution
#'  from which standard deviations of the serial interval are drawn (see
#'  details).}
#'
#' \item{max_std_si}{For method "uncertain_si" ; upper bound of the distribution
#'  from which standard deviations of the serial interval are drawn (see
#'  details).}
#'
#' \item{si_distr}{For method "non_parametric_si" ; vector of probabilities
#' giving the discrete distribution of the serial interval, starting with
#' \code{si_distr[1]} (probability that the serial interval is zero), which
#' should be zero.}
#'
#' \item{si_parametric_distr}{For method "si_from_data" ; the parametric
#' distribution to use when estimating the serial interval from data on dates of
#'  symptoms of pairs of infector/infected individuals (see details).
#' Should be one of "G" (Gamma), "W" (Weibull), "L" (Lognormal), "off1G" (Gamma
#' shifted by 1), "off1W" (Weibull shifted by 1), or "off1L" (Lognormal shifted
#' by 1).}
#'
#' \item{mcmc_control}{An object of class \code{estimate_R_mcmc_control}, as
#' returned by function \code{make_mcmc_control}. }
#'
#' \item{seed}{An optional integer used as the seed for the random number
#' generator at the start of the function (then potentially reset within the
#' MCMC for method \code{si_from_data}); useful to get reproducible results.}
#'
#' \item{mean_prior}{A positive number giving the mean of the common prior
#' distribution for all reproduction numbers (see details).}
#'
#' \item{std_prior}{A positive number giving the standard deviation of the
#' common prior distribution for all reproduction numbers (see details).}
#'
#' \item{cv_posterior}{A positive number giving the aimed posterior coefficient
#' of variation (see details).}
#'
#' }
#' @param incid As in function\code{estimate_R}.
#' @param method As in function\code{estimate_R}.
#'
#' @details
#' Analytical estimates of the reproduction number for an epidemic over
#' predefined time windows can be obtained using function \code{estimate_R},
#' for a given discrete distribution of the serial interval. \code{make_config}
#' allows to generate a configuration specifying the way the estimation will
#' be performed.
#'
#' The more incident cases are observed over a time window, the smallest the
#' posterior coefficient of variation (CV, ratio of standard deviation over
#' mean) of the reproduction number.
#' An aimed CV can be specified in the argument \code{cv_posterior}
#' (default is \code{0.3}), and a warning will be produced if the incidence
#' within one of the time windows considered is too low to get this CV.
#'
#' The methods vary in the way the serial interval distribution is specified.
#'
#' In short there are five methods to specify the serial interval distribution
#' (see below for details on each method).
#' In the first two methods, a unique serial interval distribution is
#' considered, whereas in the last three, a range of serial interval
#' distributions are integrated over:
#' \itemize{
#' \item{In method "non_parametric_si" the user specifies the discrete
#' distribution of the serial interval}
#' \item{In method "parametric_si" the user specifies the mean and sd of the
#' serial interval}
#' \item{In method "uncertain_si" the mean and sd of the serial interval are
#' each drawn from truncated normal distributions, with parameters specified by
#' the user}
#' \item{In method "si_from_data", the serial interval distribution is directly
#' estimated, using MCMC, from interval censored exposure data, with data
#' provided by the user together with a choice of parametric distribution for
#' the serial interval}
#' \item{In method "si_from_sample", the user directly provides the sample of
#' serial interval distribution to use for estimation of R. This can be a useful
#'  alternative to the previous method, where the MCMC estimation of the serial
#'  interval distribution could be run once, and the same estimated SI
#'  distribution then used in estimate_R in different contexts, e.g. with
#'  different time windows, hence avoiding to rerun the MCMC everytime
#'  estimate_R is called.}
#' }
#'
#' ----------------------- \code{method "non_parametric_si"} -------------------
#'
#' The discrete distribution of the serial interval is directly specified in the
#'  argument \code{si_distr}.
#'
#' ----------------------- \code{method "parametric_si"} -----------------------
#'
#' The mean and standard deviation of the continuous distribution of the serial
#' interval are given in the arguments \code{mean_si} and \code{std_si}.
#' The discrete distribution of the serial interval is derived automatically
#' using \code{\link{discr_si}}.
#'
#' ----------------------- \code{method "uncertain_si"} -----------------------
#'
#' \code{Method "uncertain_si"} allows accounting for uncertainty on the serial
#' interval distribution as described in Cori et al. AJE 2013.
#' We allow the mean \eqn{\mu} and standard deviation \eqn{\sigma} of the serial
#'  interval to vary according to truncated normal distributions.
#' We sample \code{n1} pairs of mean and standard deviations,
#' \eqn{(\mu^{(1)},\sigma^{(1)}),...,(\mu^{(n_2)},\sigma^{(n_2)})}, by first
#' sampling the mean \eqn{\mu^{(k)}}
#' from its truncated normal distribution (with mean \code{mean_si}, standard
#' deviation \code{std_mean_si}, minimum \code{min_mean_si} and maximum
#' \code{max_mean_si}),
#' and then sampling the standard deviation \eqn{\sigma^{(k)}} from its
#' truncated normal distribution
#' (with mean \code{std_si}, standard deviation \code{std_std_si}, minimum
#' \code{min_std_si} and maximum \code{max_std_si}), but imposing that
#' \eqn{\sigma^{(k)}<\mu^{(k)}}.
#' This constraint ensures that the Gamma probability density function of the
#' serial interval is null at \eqn{t=0}.
#' Warnings are produced when the truncated normal distributions are not
#' symmetric around the mean.
#' For each pair \eqn{(\mu^{(k)},\sigma^{(k)})}, we then draw a sample of size
#' \code{n2} in the posterior distribution of the reproduction number over each
#' time window, conditionally on this serial interval distribution.
#' After pooling, a sample of size \eqn{\code{n1}\times\code{n2}} of the joint
#' posterior distribution of the reproduction number over each time window is
#' obtained.
#' The posterior mean, standard deviation, and 0.025, 0.05, 0.25, 0.5, 0.75,
#' 0.95, 0.975 quantiles of the reproduction number for each time window are
#' obtained from this sample.
#'
#' ----------------------- \code{method "si_from_data"} -----------------------
#'
#' \code{Method "si_from_data"} allows accounting for uncertainty on the serial
#' interval distribution.
#' Unlike method "uncertain_si", where we arbitrarily vary the mean and std of
#' the SI in truncated normal distributions,
#' here, the scope of serial interval distributions considered is directly
#' informed by data
#' on the (potentially censored) dates of symptoms of pairs of infector/infected
#'  individuals.
#' This data, specified in argument \code{si_data}, should be a dataframe with 5
#'  columns:
#' \itemize{
#' \item{EL: the lower bound of the symptom onset date of the infector (given as
#'  an integer)}
#' \item{ER: the upper bound of the symptom onset date of the infector (given as
#'  an integer). Should be such that ER>=EL. If the dates are known exactly use
#'   ER = EL}
#' \item{SL: the lower bound of the symptom onset date of the infected
#' individual (given as an integer)}
#' \item{SR: the upper bound of the symptom onset date of the infected
#' individual (given as an integer). Should be such that SR>=SL.
#' If the dates are known exactly use SR = SL}
#' \item{type (optional): can have entries 0, 1, or 2, corresponding to doubly
#' interval-censored, single interval-censored or exact observations,
#' respectively, see Reich et al. Statist. Med. 2009. If not specified, this
#' will be automatically computed from the dates}
#' }
#' Assuming a given parametric distribution for the serial interval distribution
#'  (specified in si_parametric_distr),
#' the posterior distribution of the serial interval is estimated directly from
#' these data using MCMC methods implemented in the package
#' \code{coarsedatatools}.
#' The argument \code{mcmc_control} is a list of characteristics which control
#' the MCMC.
#' The MCMC is run for a total number of iterations of
#' \code{mcmc_control$burnin + n1*mcmc_control$thin};
#' but the output is only recorded after the burnin, and only 1 in every
#' \code{mcmc_control$thin} iterations,
#' so that the posterior sample size is \code{n1}.
#' For each element in the posterior sample of serial interval distribution,
#' we then draw a sample of size \code{n2} in the posterior distribution of the
#' reproduction number over each time window,
#' conditionally on this serial interval distribution.
#' After pooling, a sample of size \eqn{\code{n1}\times\code{n2}} of the joint
#' posterior distribution of
#' the reproduction number over each time window is obtained.
#' The posterior mean, standard deviation, and 0.025, 0.05, 0.25, 0.5, 0.75,
#' 0.95, 0.975 quantiles of the reproduction number for each time window are
#' obtained from this sample.
#'
#' ----------------------- \code{method "si_from_sample"} ----------------------
#'
#' \code{Method "si_from_sample"} also allows accounting for uncertainty on the
#' serial interval distribution.
#' Unlike methods "uncertain_si" and "si_from_data", the user directly provides
#' (in argument \code{si_sample}) a sample of serial interval distribution to be
#'  explored.
#'
#'
#' @return An object of class \code{estimate_R_config} with components
#' t_start, t_end, n1, n2, mean_si, std_si,
#' std_mean_si, min_mean_si, max_mean_si, std_std_si, min_std_si, max_std_si,
#' si_distr, si_parametric_distr, mcmc_control, seed, mean_prior, std_prior,
#' cv_posterior, which can be used as an argument of function \code{estimate_R}.
#' @export
#'
#' @examples
#' \dontrun{
#' ## Note the following examples use an MCMC routine
#' ## to estimate the serial interval distribution from data,
#' ## so they may take a few minutes to run
#'
#' ## load data on rotavirus
#' data("MockRotavirus")
#'
#' ## estimate the reproduction number (method "si_from_data")
#' ## we are not specifying the time windows, so by defaults this will estimate
#' ## R on sliding weekly windows
#' incid <- MockRotavirus$incidence
#' method <- "si_from_data"
#' config <- make_config(incid = incid,
#'                      method = method,
#'                      list(si_parametric_distr = "G",
#'                      mcmc_control = make_mcmc_control(burnin = 1000,
#'                      thin = 10, seed = 1),
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
#'                      list(si_parametric_distr = "G",
#'                      mcmc_control = make_mcmc_control(burnin = 1000,
#'                      thin = 10, seed = 1),
#'                      n1 = 500,
#'                      n2 = 50,
#'                      seed = 2)))
#' plot(R_si_from_data)
#' }
make_config <- function(..., incid = NULL,
                        method = c("non_parametric_si", "parametric_si",
                                   "uncertain_si", "si_from_data",
                                   "si_from_sample")) {
  config <- list(...)
  if (length(config) == 1L && is.list(config[[1]])) {
    config <- config[[1]]
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
                   mcmc_control = make_mcmc_control(),
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
    T <- sum(idx_raw_incid)

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

