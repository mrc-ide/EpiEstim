#' Create list of MCMC control parameters
#'
#' This function is deprecated. The serial interval is estimated by maximum
#' likelihood in [estimate_R()] with method "si_from_data", so `burnin` and
#' `thin` are not used, and a warning is given if they differ from their
#' defaults. Use `seed` in [make_config()] to make the sample of serial
#' interval distributions reproducible.
#'
#' @param burnin Not used.
#' @param thin Not used.
#' @param seed An integer used as the seed for the random number generator
#'   when drawing the sample of serial interval distributions.
#' @param init_pars Starting values for the parameters of the serial interval
#'   distribution, as given by [si_start_values()].
#'
#' @return An object of class `estimate_R_mcmc_control` with components
#' `burnin`, `thin`, `seed`, `init_pars`. This can be
#' used as an argument of function [make_config()].
#'
#' @seealso [si_start_values()]
#' @export
make_mcmc_control <- function(burnin = 3000, thin = 10, 
                              seed = as.integer(Sys.time()), 
                              init_pars = NULL){
  mcmc_control <- list(init_pars = init_pars, 
                       burnin = burnin, 
                       thin = thin, 
                       seed = seed )
  class(mcmc_control) <- "estimate_R_mcmc_control"
  return( mcmc_control )
}