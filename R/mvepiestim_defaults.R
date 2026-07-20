#' Set default for Gamma priors
#'
#' @return a list of default parameters for the priors.
#'   Values can then be manually edited as in the examples below.
#'   Users could use functions [epitrix::gamma_shapescale2mucv()] and
#'   [epitrix::gamma_mucv2shapescale()] to set the shape and scale corresponding
#'   to the desired prior mean and coefficient of variation.
#' @export
#'
#' @examples
#' priors <- default_priors()
#' # change the prior for R to have a mean of 3
#' priors$R$shape <- 3

default_priors <- function() {
  ## Flatter epsilon and R with larger variance
  ## Mean epsilon 10 and SD 10
  ## Mean R 1 and SD 5
  list(epsilon = list(shape = 1, scale = 1),
       R = list(shape = 0.04, scale = 25))
}

#' Set default for MCMC control
#'
#' @return a list of default MCMC control parameters, containing:
#' - `n_iter`: the number of iterations for the MCMC to perform
#' - `burnin`: the burnin to use; MCMC iterations will only be recorded after
#'   the burnin
#' - `thin`: MCMC iterations will only be recorded after the burnin and every 
#'   `thin` iteration
#'
#' Values can then be manually edited as in the examples below.
#'
#' @export
#'
#' @examples
#' mcmc_control <- default_mcmc_controls()
#' # change to run for 10 times longer
#' mcmc_control$n_iter <- mcmc_control$n_iter * 10

default_mcmc_controls <- function() {
  list(n_iter = 1100L,
       burnin = 10L,
       thin = 10L)
}

