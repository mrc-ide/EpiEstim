#' Title ### PLEASE USE THIS ONLY FOR COMPATIBILITY...
#'
#' @param I see I in \code{estimate_r} ### TO DO change to incidence
#' @param T.Start see \code{config$t_start} in \code{estimate_r}
#' @param T.End .
#' @param method .
#' @param n1 .
#' @param n2 .
#' @param Mean.SI .
#' @param Std.SI .
#' @param Std.Mean.SI .
#' @param Min.Mean.SI .
#' @param Max.Mean.SI .
#' @param Std.Std.SI .
#' @param Min.Std.SI .
#' @param Max.Std.SI .
#' @param SI.Distr .
#' @param SI.Data .
#' @param SI.parametricDistr .
#' @param MCMC.control .
#' @param SI.Sample .
#' @param seed .
#' @param Mean.Prior .
#' @param Std.Prior .
#' @param CV.Posterior .
#' @param plot .
#' @param legend .
#'
#' @export
#'
EstimateR <- function(I, T.Start, T.End, method = c("NonParametricSI", "ParametricSI",
                                                    "UncertainSI", "SIFromData", "SIFromSample"), 
                      n1 = NULL, n2 = NULL, Mean.SI = NULL, Std.SI = NULL,
                      Std.Mean.SI = NULL, Min.Mean.SI = NULL, Max.Mean.SI = NULL,
                      Std.Std.SI = NULL, Min.Std.SI = NULL, Max.Std.SI = NULL,
                      SI.Distr = NULL, 
                      SI.Data = NULL, SI.parametricDistr = c("G", "W", "L", "off1G", "off1W", "off1L"),  
                      MCMC.control = list(init.pars = NULL, burnin = 3000, thin=10, seed = as.integer(Sys.time())), 
                      SI.Sample = NULL, 
                      seed = NULL,
                      Mean.Prior = 5, Std.Prior = 5, CV.Posterior = 0.3,
                      plot = FALSE, legend = FALSE) {
  .Deprecated("estimate_r")
  
  method_tr <- c("NonParametricSI" = "non_parametric_si", 
                 "ParametricSI" = "parametric_si",
                 "UncertainSI" = "uncertain_si", 
                 "SIFromData" = "si_from_data", 
                 "SIFromSample" = "si_from_sample")
  method <- method_tr[[method]]

  config <- list(
    t_start = T.Start,
    t_end = T.End,
    n1 = n1,
    n2 = n2,
    mean_si = Mean.SI,
    std_si = Std.SI,
    std_mean_si = Std.Mean.SI,
    min_mean_si = Min.Mean.SI,
    max_mean_si = Max.Mean.SI,
    std_std_si = Std.Std.SI,
    min_std_si = Min.Std.SI,
    max_std_si = Max.Std.SI,
    si_distr = SI.Distr,
    si_parametric_distr = SI.parametricDistr,
    mcmc_control = MCMC.control,
    seed = seed,
    mean_prior = Mean.Prior,
    std_prior = Std.Prior, ### TO DO: change to sd_prior
    cv_posterior = CV.Posterior,
    plot = plot, ### TO DO: plot outside the function (alsways!)
    legend = legend
  )
    
  estimate_r(I, method, si_data = SI.Data, si_sample = SI.Sample, 
             config = config)
}