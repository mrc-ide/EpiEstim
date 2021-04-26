#' Set default for Gamma priors
#'
#' @return a list of default parameters for the priors.
#'   Values can then be manually be edited as in the examples below.
#'   Users could use functions `epitrix::gamma_shapescale2mucv` and
#'   `epitrix::gamma_mucv2shapescale` to set the shape and scale corresponding
#'   to the desired prior mean and coefficient of variation.
#' @export
#'
#' @examples
#' priors <- default_priors()
#' # change the prior for R to have a mean of 3
#' priors$R$shape <- 3
#'
default_priors <- function() {
  list(epsilon = list(shape = 1, scale = 1),
       R = list(shape = 1, scale = 1))
}


#' Set default for MCMC controls
#'
#' @return a list of default MCMC control parameters, containing:
#'
#' - n_iter: the number if iterations of the MCMC to perform
#'
#' - burnin: the burnin to use; MCMC iterations will only be recorded after
#'   the burnin
#'
#' - thin: MCMC iterations will only be recorded after
#'   the burnin and every `thin` iteration
#'
#'   Values can then be manually be edited as in the examples below.
#'
#' @export
#'
#' @examples
#' mcmc_controls <- default_mcmc_controls()
#' # change to run for 10 times longer
#' mcmc_controls$n_iter <- mcmc_controls$n_iter * 10
#'

default_mcmc_controls <- function() {
  if (n_iter < 0 | !is.integer(n_iter)){
    stop("n_iter must be a positive integer")
  }
  if (burnin < 0 | !is.integer(burnin)){
    stop("burnin must be a positive integer")
  }
  if (thin < 0 | !is.integer(thin)){
    stop("thin must be a positive integer")
  }
  if (n_iter < burnin + thin){
    stop("n_iter must be greater than burnin + thin")
  }
  list(n_iter = 1100, burnin = 10, thin = 10)
}


#' Compute the overall infectivity
#'
#' @param incid a multidimensional array containing values of the incidence
#'   for each time step (1st dimension), location (2nd dimension) and
#'   pathogen/strain/variant (3rd dimension)
#'
#' @param si_distr a matrix where each column contains the probability mass
#'   function for the discrete serial interval for each of the
#'   pathogen/strain/variants, starting with the probability mass function
#'   for day 0 in the first row, which should be 0. Each column in the matrix
#'   should sum to 1
#'
#' @return a multidimensional array containing values of the overall
#'   infectivity for each time step (1st dimension), location (2nd dimension)
#'   and pathogen/strain/variant (3rd dimension). The overall infectivity for
#'   a given location and pathogen/strain/variant represents the sum of
#'   the incidence for that location and that pathogen/strain/variant at all
#'   previous time steps, weighted by the current infectivity of those
#'   past incident cases. Pre-calculating the overall infectivity makes the
#'   algorithm much faster
#'
#' @export
#'
#' @examples
#'
#' n_v <- 2
#' n_loc <- 3 # 3 locations
#' T <- 100 # 100 time steps
#' priors <- default_priors()
#' # constant incidence 10 per day everywhere
#' incid <- array(10, dim = c(T, n_loc, n_v))
#' # arbitrary serial interval, same for both variants
#' w_v <- c(0, 0.2, 0.5, 0.3)
#' si_distr <- cbind(w_v, w_v)
#' lambda <- compute_lambda(incid, si_distr)
compute_lambda <- function(incid, si_distr) {
  ## TODO: check that si_distr[1, ] == 0
  ## TODO: check that all(colSums(si_distr) == 1)
  ## TODO: check that all(si_distr >= 0)
  lambda <- array(NA, dim = dim(incid))
  for(l in seq_len(dim(incid)[2])) {
    for(v in seq_len(dim(incid)[3])) {
      lambda[, l, v] <- EpiEstim::overall_infectivity(incid[, l, v], si_distr[, v])
    }
  }
  lambda
}


#' Draw epsilon from marginal posterior distribution
#'
#' @param R a matrix with dimensions containing values of the instantaneous
#'   reproduction number for each time step (row) and location (column), for
#'   the reference pathogen/strain/variant
#'
#' @param incid a multidimensional array containing values of the incidence
#'   for each time step (1st dimension), location (2nd dimension) and
#'   pathogen/strain/variant (3rd dimension)
#'
#' @param lambda a multidimensional array containing values of the overall
#'   infectivity for each time step (1st dimension), location (2nd dimension)
#'   and pathogen/strain/variant (3rd dimension). The overall infectivity for
#'   a given location and pathogen/strain/variant represents the sum of
#'   the incidence for that location and that pathogen/strain/variant at all
#'   previous time steps, weighted by the current infectivity of those
#'   past incident cases. It can be calculated from the incidence `incid` and
#'   the distribution of the serial interval using function `compute_lambda`)
#'
#' @param priors a list of prior parameters (shape and scale of a gamma
#'   distribution) for epsilon and R; can be obtained from the function
#'   `default_priors`. The prior for R is assumed to be the same for all
#'   time steps and all locations
#'
#' @param t_min an integer >1 giving the minimum time step to consider in the
#'   estimation. Default value is 2 (as the estimation is conditional on
#'   observations at time step 1 and can therefore only start at time step 2).
#'
#' @param t_max an integer >`t_min` and <=`nrow(incid)` giving the maximum time
#'   step to consider in the estimation. Default value is `nrow(incid)`.
#'
#' @param seed a numeric value used to fix the random seed
#'
#' @return a value or vector of values for epsilon for each non reference
#'   pathogen/strain/variant, drawn from the marginal posterior distribution
#'
#' @importFrom("stats", "median", "rgamma")
#'
#' @export
#'
#' @examples
#'
#' n_loc <- 4 # 4 locations
#' n_v <- 3 # 3 strains
#' T <- 100 # 100 time steps
#' priors <- default_priors()
#' # constant incidence 10 per day everywhere
#' incid <- array(10, dim = c(T, n_loc, n_v))
#' # arbitrary serial interval, same for both variants
#' w_v <- c(0, 0.2, 0.5, 0.3)
#' si_distr <- cbind(w_v, w_v, w_v)
#' lambda <- compute_lambda(incid, si_distr)
#' # Constant reproduction number of 1
#' R <- matrix(1, nrow = T, ncol = n_loc)
#' R[1, ] <- NA # no estimates of R on first time step
#' draw_epsilon(R, incid, lambda, priors, seed = 1)
#'
draw_epsilon <- function(R, incid, lambda, priors,
                         t_min = 2, t_max = nrow(incid),
                         seed = NULL) {
  if (!is.integer(t_min) | !is.integer(t_max)){
    stop("t_min and t_max must be integers")
  }
  if (t_min < 2 | t_max < 2){
    stop("t_min and t_max must be >=2")
  }
  if(t_min > nrow(incid) | t_max > nrow(incid)){
    stop("t_min and t_max must be <= nrow(incid)")
  }
  if (!is.numeric(seed)){
    stop("seed must be numeric")
  }
  if (R < 0){
    stop("R must be >=0")
  }
  if (!is.null(seed)) set.seed(seed)
  t <- seq(t_min, t_max, 1)
  shape <- EpiEstim:::vnapply(seq(2, dim(lambda)[3]), function(e)
    sum(incid[t, , e])) + priors$epsilon$shape
  rate <- EpiEstim:::vnapply(seq(2, dim(lambda)[3]), function(e)
    sum(R[t, ] * lambda[t, , e]) + 1 / priors$epsilon$scale)
  scale <- 1 / rate
  rgamma(dim(lambda)[3] - 1, shape = shape, scale = scale)
}

#' Draw R from marginal posterior distribution
#'
#' @param epsilon a value or vector of values for the relative transmissibility
#'   of the "new" pathogen/strain/variant(s) compared to the reference
#'   pathogen/strain/variant
#'
#' @param incid a multidimensional array containing values of the incidence
#'   for each time step (1st dimension), location (2nd dimension) and
#'   pathogen/strain/variant (3rd dimension)
#'
#' @param lambda a multidimensional array containing values of the overall
#'   infectivity for each time step (1st dimension), location (2nd dimension)
#'   and pathogen/strain/variant (3rd dimension). The overall infectivity for
#'   a given location and pathogen/strain/variant represents the sum of
#'   the incidence for that location and that pathogen/strain/variant at all
#'   previous time steps, weighted by the current infectivity of those
#'   past incident cases. It can be calculated from the incidence `incid` and
#'   the distribution of the serial interval using function `compute_lambda`)
#'
#' @param priors a list of prior parameters (shape and scale of a gamma
#'   distribution) for epsilon and R; can be obtained from the function
#'   `default_priors`. The prior for R is assumed to be the same for all
#'   time steps and all locations
#'
#' @param t_min an integer >1 giving the minimum time step to consider in the
#'   estimation. Default value is 2 (as the estimation is conditional on
#'   observations at time step 1 and can therefore only start at time step 2).
#'
#' @param t_max an integer >`t_min` and <=`nrow(incid)` giving the maximum time
#'   step to consider in the estimation. Default value is `nrow(incid)`.
#'
#' @param seed a numeric value used to fix the random seed
#'
#' @return a matrix of the instantaneous reproduction number R for the reference
#'   pathogen/strain/variant for each time step (row) and each location (column)
#'   drawn from the marginal posterior distribution
#'
#' @export
#'
#' @importFrom("stats", "median", "rgamma")
#'
#' @examples
#'
#' n_v <- 2
#' n_loc <- 3 # 3 locations
#' T <- 100 # 100 time steps
#' priors <- default_priors()
#' # constant incidence 10 per day everywhere
#' incid <- array(10, dim = c(T, n_loc, n_v))
#' # arbitrary serial interval, same for both variants
#' w_v <- c(0, 0.2, 0.5, 0.3)
#' si_distr <- cbind(w_v, w_v)
#' lambda <- compute_lambda(incid, si_distr)
#' # Epsilon = 1 i.e. no transmission advantage
#' epsilon <- 1
#' draw_R(epsilon, incid, lambda, priors, seed = 1)
#'
draw_R <- function(epsilon, incid, lambda, priors,
                   t_min = 2, t_max = nrow(incid),
                   seed = NULL) {
  ## TODO: check t_min and t_max are integers, >=2 and <= nrow(incid)
  ## TODO: check seed is a numeric value
  ## TODO: check epsilon >0
  if (!is.null(seed)) set.seed(seed)
  t <- seq(t_min, t_max, 1)
  shape <- apply(incid[t, , ], c(1, 2), sum) + priors$R$shape ## TODO: precalculate this
  shape_flat <- as.numeric(shape) ## TODO: precalculate this
  rate <- lambda[t, , 1] + apply(epsilon * lambda[t, , -1], c(1, 2), sum) +
    1 / priors$R$scale
  scale <- 1 / rate
  scale_flat <- as.numeric(scale)
  R_flat <- rgamma(length(shape_flat), shape = shape_flat, scale = scale_flat)
  R_fill <- matrix(R_flat, nrow = nrow(shape), ncol = ncol(shape))
  R <- matrix(NA, nrow(incid), ncol(incid))
  R[t, ] <- R_fill
  R
}


#' Jointly estimate the instantaneous reproduction number for a reference
#'   pathogen/strain/variant and the relative transmissibility of a
#'   "new" pathogen/strain/variant
#'
#' @param incid a multidimensional array containing values of the incidence
#'   for each time step (1st dimension), location (2nd dimension) and
#'   pathogen/strain/variant (3rd dimension)
#'
#' @param si_distr a matrix with two columns, each containing the probability mass
#'   function for the discrete serial interval for each of the two
#'   pathogen/strain/variants, starting with the probability mass function
#'   for day 0 in the first row, which should be 0. each column in the matrix
#'   should sum to 1
#'
#' @param priors a list of prior parameters (shape and scale of a gamma
#'   distribution) for epsilon and R; can be obtained from the function
#'   `default_priors`. The prior for R is assumed to be the same for all
#'   time steps and all locations
#'
#' @param mcmc_control a list of default MCMC control parameters, as obtained
#'   for example from function `default_mcmc_controls`
#'
#' @param t_min an integer >1 giving the minimum time step to consider in the
#'   estimation. Default value is 2 (as the estimation is conditional on
#'   observations at time step 1 and can therefore only start at time step 2).
#'
#' @param t_max an integer >`t_min` and <=`nrow(incid)` giving the maximum time
#'   step to consider in the estimation. Default value is `nrow(incid)`.
#'
#' @param seed a numeric value used to fix the random seed
#'
#' @return a list with two elements.
#'   1) `epsilon` is a matrix containing the MCMC chain (thinned and after
#'   burnin) for the relative transmissibility of the "new"
#'   pathogen/strain/variant(s) compared to the reference
#'   pathogen/strain/variant. Each row in the matrix is a "new"
#'   pathogen/strain/variant and each column an iteration of the MCMC.
#'   2) `R_out` is an array containing the MCMC chain (thinned and after
#'   burnin) for the reproduction number for the reference
#'   pathogen/strain/variant. The first dimension of the array is time,
#'   the second location, and the third iteration of the MCMC.
#'
#' @export
#'
#' @importFrom("stats", "median", "rgamma")
#'
#' @examples
#'
#' n_v <- 2
#' n_loc <- 3 # 3 locations
#' T <- 100 # 100 time steps
#' priors <- default_priors()
#' # constant incidence 10 per day everywhere
#' incid <- array(10, dim = c(T, n_loc, n_v))
#' # arbitrary serial interval, same for both variants
#' w_v <- c(0, 0.2, 0.5, 0.3)
#' si_distr <- cbind(w_v, w_v)
#'
#' # Dummy initial values for the MCMC
#' R_init <- matrix(5, nrow = T, ncol = n_loc)
#' R_init[1, ] <- NA # no estimates of R on first time step
#' epsilon_init <- 5
#' x <- estimate_joint(incid, si_distr, priors)
#' # Plotting to check outputs
#' par(mfrow = c(2, 2))
#' plot(x$epsilon, type = "l",
#'      xlab = "Iteration", ylab = "epsilon")
#' # Compare with what we expect with constant incidence in all locations
#' abline(h = 1, col = "red")
#' plot(x$R[10, 1, ], type = "l",
#'      xlab = "Iteration", ylab = "R time 10 location 1")
#' abline(h = 1, col = "red")
#' plot(x$R[20, 2, ], type = "l",
#'      xlab = "Iteration", ylab = "R time 20 location 2")
#' abline(h = 1, col = "red")
#' plot(x$R[30, 3, ], type = "l",
#'      xlab = "Iteration", ylab = "R time 30 location 3")
#'
estimate_joint <- function(incid, si_distr, priors,
                           mcmc_control = default_mcmc_controls(),
                           t_min = 2, t_max = nrow(incid),
                           seed = NULL
) {
  ## TODO: check t_min and t_max are integers, >=2 and <= nrow(incid)
  ## TODO: check seed is a numeric value
  ## TODO: check si_distr has the right format
  ## TODO: check mcmc_control has the correct format
  if (!is.null(seed)) set.seed(seed)
  t <- seq(t_min, t_max, 1)

  T <- nrow(incid)
  n_loc <- ncol(incid)

  lambda <- compute_lambda(incid, si_distr)

  ## find clever initial values, based on ratio of reproduction numbers
  ## in first location
  R_init <- lapply(seq_len(dim(incid)[3]), function(i) suppressWarnings(
    EpiEstim::estimate_R(incid[, 1, i], method = "non_parametric_si",
                         config = EpiEstim::make_config(si_distr = si_distr[, i],
                                              t_start = t,
                                              t_end = t)))$R$'Mean(R)')
  epsilon_init <- unlist(lapply(seq(2, length(R_init)), function(i)
    median(R_init[[i]] / R_init[[1]], na.rm = TRUE)))
  epsilon_out <- matrix(NA, nrow = length(epsilon_init),
                        ncol = mcmc_control$n_iter + 1)
  epsilon_out[, 1] <- epsilon_init
  R_init <- draw_R(mcmc_control$n_iter, incid, lambda, priors,
                   t_min = t_min, t_max = t_max)
  R_out <- array(NA, dim= c(T, n_loc, mcmc_control$n_iter + 1))
  R_out[, , 1] <- R_init

  for (i in seq_len(mcmc_control$n_iter)) {
    R_out[, , i + 1] <- draw_R(epsilon_out[, i], incid, lambda, priors,
                               t_min = t_min, t_max = t_max)
    epsilon_out[, i + 1] <- draw_epsilon(R_out[, , i + 1], incid, lambda, priors,
                                       t_min = t_min, t_max = t_max)
  }

  # remove burnin and thin
  keep <- seq(mcmc_control$burnin, mcmc_control$n_iter, mcmc_control$thin)
  epsilon_out <- epsilon_out[, keep]
  R_out <- R_out[, , keep]

  list(epsilon = epsilon_out, R = R_out)

}



## TODO: check dimensions of objects is correct everywhere
## TODO: fix number of variants to be 2

