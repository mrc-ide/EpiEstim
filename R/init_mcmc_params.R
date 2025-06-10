################################################################################
# init_mcmc_params finds clever starting points for the MCMC to be used to
# estimate the serial interval, when using option si_from_data in estimate_R # 
################################################################################

#' \code{init_mcmc_params} Finds clever starting points for the MCMC to be used to 
#' estimate the serial interval, e.g. when using option \code{si_from_data} in 
#' \code{estimate_R}
#' 
#' \code{init_mcmc_params} Finds values of the serial interval distribution 
#' parameters, used to initialise the MCMC estimation of the serial interval 
#' distribution. Initial values are computed based on the observed mean and 
#' standard deviation of the sample from which the parameters are to be 
#' estimated.
#' 
#' @param si_data the data on dates of symptoms of pairs of infector/infected
#'   individuals to be used to estimate the serial interval distribution. This
#'   should be a dataframe with 5 columns: 
#'   * EL: the lower bound
#'   of the symptom onset date of the infector (given as an integer)
#'   * ER:
#'   the upper bound of the symptom onset date of the infector (given as an
#'   integer). Should be such that ER>=EL
#'   * SL: the lower bound of the
#'   symptom onset date of the infected individual (given as an integer)
#'   * SR: the upper bound of the symptom onset date of the infected
#'   individual (given as an integer). Should be such that SR>=SL
#'   * type
#'   (optional): can have entries 0, 1, or 2, corresponding to doubly
#'   interval-censored, single interval-censored or exact observations, 
#'   respectively, see Reich et al. Statist. Med. 2009. If not specified, this
#'   will be automatically computed from the dates
#' @param dist the parametric distribution to use for the serial interval. 
#'   Should be one of "G" (Gamma), "W" (Weibull), "L" (Lognormal), "off1G"
#'   (Gamma shifted by 1), "off1W" (Weibull shifted by 1), or "off1L" (Lognormal
#'   shifted by 1).
#' @return A vector containing the initial values for the two parameters of the
#'   distribution of the serial interval. These are the shape and scale for all
#'   but the lognormal distribution, for which it is the meanlog and sdlog.
#' @seealso  \code{\link{estimate_R}}
#' @author Anne Cori
#' @importFrom fitdistrplus fitdist
#' @export
#' @examples
#' \dontrun{
#' ## Note the following examples use an MCMC routine
#' ## to estimate the serial interval distribution from data,
#' ## so they may take a few minutes to run
#' 
#' ## load data on rotavirus
#' data("MockRotavirus")
#' 
#' ## get clever initial values for shape and scale of a Gamma distribution
#' ## fitted to the the data MockRotavirus$si_data
#' clever_init_param <- init_mcmc_params(MockRotavirus$si_data, "G")
#' 
#' ## estimate the serial interval from data using a clever starting point for 
#' ## the MCMC chain
#' SI_fit_clever <- coarseDataTools::dic.fit.mcmc(dat = MockRotavirus$si_data,
#'                              dist = "G",
#'                              init.pars = clever_init_param,
#'                              burnin = 1000,
#'                              n.samples = 5000)
#' 
#' ## estimate the serial interval from data using a random starting point for 
#' ## the MCMC chain
#' SI_fit_naive <- coarseDataTools::dic.fit.mcmc(dat = MockRotavirus$si_data,
#'                              dist = "G",
#'                              burnin = 1000,
#'                              n.samples = 5000)
#' 
#' 
#' ## use check_cdt_samples_convergence to check convergence in both situations
#' converg_diag_clever <- check_cdt_samples_convergence(SI_fit_clever@samples)
#' converg_diag_naive <- check_cdt_samples_convergence(SI_fit_naive@samples)
#' converg_diag_clever
#' converg_diag_naive
#' 
#' }
#' 
init_mcmc_params <- function(si_data, 
                             dist = c("G", "W", "L", "off1G", 
                                      "off1W", "off1L")) {
  dist <- match.arg(dist)
  naive_SI_obs <- (si_data$SR + si_data$SL) / 2 - (si_data$ER + si_data$EL) / 2
  mu <- mean(naive_SI_obs)
  sigma <- sd(naive_SI_obs)
  if (dist == "G") {
    shape <- (mu / sigma)^2
    scale <- sigma^2 / mu
    # check this is what we want
    # tmp <- rgamma(10000, shape=shape, scale = scale)
    # mean(tmp)
    # sd(tmp)
    param <- c(shape, scale)
  } else if (dist == "W") {
    fit.w <- fitdist(naive_SI_obs + 0.1, "weibull") 
    ## using +0.1 to avoid issues with zero
    shape <- fit.w$estimate["shape"]
    scale <- fit.w$estimate["scale"]
    # check this is what we want
    # tmp <- rweibull(10000, shape=shape, scale = scale)
    # mean(tmp)
    # sd(tmp)
    param <- c(shape, scale)
  } else if (dist == "L") {
    sdlog <- sqrt(log(sigma^2 / (mu^2) + 1))
    meanlog <- log(mu) - sdlog^2 / 2
    # check this is what we want
    # tmp <- rlnorm(10000, meanlog=meanlog, sdlog = sdlog)
    # mean(tmp)
    # sd(tmp)
    param <- c(meanlog, sdlog)
  } else if (dist == "off1G") {
    shape <- ((mu - 1) / sigma)^2
    if (shape <= 0) shape <- 0.001 
    ## this is to avoid issues when the mean SI is <1
    scale <- sigma^2 / (mu - 1)
    # check this is what we want
    # tmp <- 1+rgamma(10000, shape=shape, scale = scale)
    # mean(tmp)
    # sd(tmp)
    param <- c(shape, scale)
  } else if (dist == "off1W") {
    fit.w <- fitdist(naive_SI_obs - 1 + 0.1, "weibull") 
    ## using +0.1 to avoid issues with zero
    shape <- fit.w$estimate["shape"]
    scale <- fit.w$estimate["scale"]
    # check this is what we want
    # tmp <- 1+rweibull(10000, shape=shape, scale = scale)
    # mean(tmp)
    # sd(tmp)
    param <- c(shape, scale)
  } else if (dist == "off1L") {
    sdlog <- sqrt(log(sigma^2 / ((mu - 1)^2) + 1))
    meanlog <- log(mu - 1) - sdlog^2 / 2
    # check this is what we want
    # tmp <- 1+rlnorm(10000, meanlog=meanlog, sdlog = sdlog)
    # mean(tmp)
    # sd(tmp)
    param <- c(meanlog, sdlog)
  } else {
    stop(sprintf("Distribtion (%s) not supported", dist))
  }
  if (any(is.na(param))) {
    stop("NA result. Check that si_data is in the right format. ")
  }
  return(param)
}
