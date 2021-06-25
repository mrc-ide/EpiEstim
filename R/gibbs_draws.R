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
  ## Flatter epsilon and R with larger variace
  ## Mean epsilon 10 and SD 10
  ## Mean R 1 and SD 5
  list(epsilon = list(shape = 1, scale = 10),
       R = list(shape = 0.04, scale = 25))
}


#' Set default for MCMC control
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
#' mcmc_control<- default_mcmc_controls()
#' # change to run for 10 times longer
#' mcmc_control$n_iter <- mcmc_control$n_iter * 10
#'

default_mcmc_controls <- function() {
  list(n_iter = 1100L,
       burnin = 10L,
       thin = 10L)
}


#' Compute the overall infectivity
#'
#' @param incid a list (as obtained from function `process_I_multivariant`)
#'   of two multidimensional arrays ("local" and "imported")
#'   containing values of the incidence
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
#' incid <- process_I_multivariant(incid)
#' # arbitrary serial interval, same for both variants
#' w_v <- c(0, 0.2, 0.5, 0.3)
#' si_distr <- cbind(w_v, w_v)
#' lambda <- compute_lambda(incid, si_distr)

compute_lambda <- function(incid, si_distr) {
  if (!inherits(incid, "incid_multivariant")) {
    msg1 <- "'incid 'should be an 'incid_multivariant' object."
    msg2 <- "Use function 'process_I_multivariant' first"
    stop(paste(msg1, msg2))
  }
  if (any(si_distr[1,] != 0)){
    stop("Values in the first row of si_distr must be 0")
  }
  if (any(abs(colSums(si_distr) - 1) > 0.01)) { # allow tolerance
    stop("The sum of each column in si_distr should be equal to 1")
  }
  if (any(si_distr < 0)){
    stop("si_distr must be >=0")
  }
  lambda <- array(NA, dim = dim(incid$local))
  for(l in seq_len(dim(incid$local)[2])) {
    for(v in seq_len(dim(incid$local)[3])) {
      lambda[, l, v] <- EpiEstim::overall_infectivity(
        data.frame(local = incid$local[, l, v],
                   imported = incid$imported[, l, v]), si_distr[, v])
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
#' @param incid a multidimensional array containing values of the (local)
#'   incidence
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
#' @importFrom stats median rgamma
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
#' incid <- process_I_multivariant(incid)
#' # arbitrary serial interval, same for both variants
#' w_v <- c(0, 0.2, 0.5, 0.3)
#' si_distr <- cbind(w_v, w_v, w_v)
#' lambda <- compute_lambda(incid, si_distr)
#' # Constant reproduction number of 1
#' R <- matrix(1, nrow = T, ncol = n_loc)
#' R[1, ] <- NA # no estimates of R on first time step
#' draw_epsilon(R, incid$local, lambda, priors, seed = 1)
#'
draw_epsilon <- function(R, incid, lambda, priors,
                         t_min = 2L, t_max = nrow(incid),
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
  if(any(R[!is.na(R)] < 0)) {
    stop("R must be >= 0")
  }
  if (!is.null(seed) & !is.numeric(seed)){
    stop("seed must be numeric")
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
#' @param incid a multidimensional array containing values of the (local)
#'   incidence
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
#' @importFrom stats median rgamma
#'
#' @examples
#'
#' n_v <- 2
#' n_loc <- 3 # 3 locations
#' T <- 100 # 100 time steps
#' priors <- default_priors()
#' # constant incidence 10 per day everywhere
#' incid <- array(10, dim = c(T, n_loc, n_v))
#' incid <- process_I_multivariant(incid)
#' # arbitrary serial interval, same for both variants
#' w_v <- c(0, 0.2, 0.5, 0.3)
#' si_distr <- cbind(w_v, w_v)
#' lambda <- compute_lambda(incid, si_distr)
#' # Epsilon = 1 i.e. no transmission advantage
#' epsilon <- 1
#' draw_R(epsilon, incid$local, lambda, priors, seed = 1)
#'
draw_R <- function(epsilon, incid, lambda, priors,
                   t_min = 2L, t_max = nrow(incid),
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
  if (any(epsilon < 0)){
    stop("epsilon must be > 0")
  }
  if (!is.null(seed) & !is.numeric(seed)){
    stop("seed must be numeric")
  }
  if (!is.null(seed)) set.seed(seed)
  t <- seq(t_min, t_max, 1)
  shape <- apply(incid[t, , , drop = FALSE], c(1, 2), sum) + priors$R$shape ## TODO: precalculate this
  shape_flat <- as.numeric(shape) ## TODO: precalculate this
  ## Fix for issue 123.
  ## overall infectivity for
  temp <- lambda[t, , 1]
  idx <- seq(2, dim(incid)[3], 1)
  for(var in idx){
    ## We want lambda_1 + e_v lambda_v for all t
    temp <- temp + epsilon[var - 1] * lambda[t, , var]
  }
  rate <- temp + 1 / priors$R$scale
  scale <- 1 / rate
  scale_flat <- as.numeric(scale)
  R_flat <- rgamma(length(shape_flat), shape = shape_flat, scale = scale_flat)
  R_fill <- matrix(R_flat, nrow = nrow(shape), ncol = ncol(shape))
  R <- matrix(NA, nrow(incid), ncol(incid))
  R[t, ] <- R_fill
  R
}
##' Index before which at most a given probability
##' mass is captured
##'
##' Across a matrix of discretised probability distributions
##' (see \code{estimate_joint}
##' this function returns the largest index
##' (across all columns) such that the
##' cumulative probability mass before index is
##' 1 - \code{miss_at_most}.
##'
##'
##' @inheritParams estimate_joint
##' @param miss_at_most
##' @return
##' @author Sangeeta Bhatia
##' @export
compute_si_cutoff <- function(si_distr, miss_at_most = 0.05) {
  cutoff <- 1 - miss_at_most
  cdf <- apply(si_distr, 2, cumsum)
  idx <- apply(
    cdf, 2,
    function(col) Position(function(x) x > cutoff, col)
  )
  max(idx)
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
#'   estimation.
#'   The NULL, t_min is calculated using the function \code{compute_si_cutoff}
#'   which gets the maximum (across all variants) of the 95th percentile of the
#'   SI distribution.
#'
#'
#' @param t_max an integer >`t_min` and <=`nrow(incid)` giving the maximum time
#'   step to consider in the estimation. Default value is `nrow(incid)`.
#'
#' @param seed a numeric value used to fix the random seed
#'
#' @param incid_imported an optional multidimensional array containing values
#'   of the incidence of imported cases
#'   for each time step (1st dimension), location (2nd dimension) and
#'   pathogen/strain/variant (3rd dimension). `incid - incid_imported` is
#'   therefore the incidence of locally infected cases. If `incid_imported` is
#'   NULL this means there are no
#'   known imported cases and all cases other than on those from the first
#'   time step will be considered locally infected.
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
#' @importFrom stats median rgamma
#' @importFrom abind adrop
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
                           t_min = NULL, t_max = nrow(incid),
                           seed = NULL,
                           incid_imported = NULL
                           ) {
  if (is.null(t_min)) {
   t_min <- compute_si_cutoff(si_distr, 0.05)
  }
  if (!is.integer(t_min) | !is.integer(t_max)){
    stop("t_min and t_max must be integers")
  }
  if (t_min < 2 | t_max < 2){
    stop("t_min and t_max must be >=2")
  }
  if(t_min > nrow(incid) | t_max > nrow(incid)){
    stop("t_min and t_max must be <= nrow(incid)")
  }
  if (any(si_distr[1,] != 0)){
    stop("Values in the first row of si_distr must be 0")
  }
  if (any(abs(colSums(si_distr) - 1) > 0.01)) { # allow tolerance
    stop("The sum of each column in si_distr should be equal to 1")
  }
  if (any(si_distr < 0)){
    stop("si_distr must be >=0")
  }
  if (mcmc_control$n_iter < 0 | !is.integer(mcmc_control$n_iter)){
    stop("n_iter in mcmc_control must be a positive integer")
  }
  if (mcmc_control$burnin < 0 | !is.integer(mcmc_control$burnin)){
    stop("burnin in mcmc_control must be a positive integer")
  }
  if (mcmc_control$thin < 0 | !is.integer(mcmc_control$thin)){
    stop("thin in mcmc_control must be a positive integer")
  }
  if (mcmc_control$n_iter < mcmc_control$burnin + mcmc_control$thin){
    stop("In mcmc_control, n_iter must be greater than burnin + thin")
  }
  if (!is.null(seed) & !is.numeric(seed)){
    stop("seed must be numeric")
  }
  if (!is.null(seed)) set.seed(seed)

  t <- seq(t_min, t_max, 1)

  T <- nrow(incid)
  n_loc <- ncol(incid)

  incid <- process_I_multivariant(incid, incid_imported)

  lambda <- compute_lambda(incid, si_distr)

  ## find clever initial values, based on ratio of reproduction numbers
  ## over the whole time period, across all locations together

  R_init <- sapply(seq_len(dim(incid$local)[3]), function(i) {
    tmp_df <- data.frame(local = apply(incid$local[, , i, drop = FALSE],
                                       c(1, 3), sum)[,1],
                         imported = apply(incid$imported[, , i, drop = FALSE],
                                          c(1, 3), sum)[,1])
    suppressWarnings(
    EpiEstim::estimate_R(tmp_df,
                         method = "non_parametric_si",
                         config = EpiEstim::make_config(
                           si_distr = si_distr[, i],
                           t_start = t_min,
                           t_end = t_max)))$R$'Mean(R)'
    })

  max_transmiss <- which.max(R_init)
  # reorder variants so most transmissible is first
  incid_reordered <- array(NA, dim = dim(incid$local))
  incid_reordered[,,1] <- incid$local[,,max_transmiss]
  incid_reordered[,,-1] <- incid$local[,, -max_transmiss]

  incid$local <- incid_reordered
  ## Re-order R_init so that they are in the same
  ## order as the incidence
  R_init_reord <- R_init
  R_init_reord[1] <- R_init[max_transmiss]
  R_init_reord[-1] <- R_init[-max_transmiss]
  R_init <- R_init_reord

  ## Re-order lambda
  lambda_reordered <- lambda
  lambda_reordered[, , 1] <- lambda[, , max_transmiss]
  lambda_reordered[, , -1] <- lambda[, , -max_transmiss]
  lambda <- lambda_reordered

  epsilon_init <- unlist(lapply(seq(2, length(R_init)), function(i)
    median(R_init[[i]] / R_init[[1]], na.rm = TRUE)))
  epsilon_out <- matrix(NA, nrow = length(epsilon_init),
                        ncol = mcmc_control$n_iter + 1)
  epsilon_out[, 1] <- epsilon_init
  R_init <- draw_R(
    epsilon = epsilon_init, incid = incid$local, lambda = lambda,
    priors = priors,t_min = t_min, t_max = t_max
  )

  R_out <- array(NA, dim= c(T, n_loc, mcmc_control$n_iter + 1))
  R_out[, , 1] <- R_init
  for (i in seq_len(mcmc_control$n_iter)) {
    R_out[, , i + 1] <- draw_R(epsilon_out[, i], incid$local, lambda, priors,
                               t_min = t_min, t_max = t_max)
    epsilon_out[, i + 1] <- draw_epsilon(
      abind::adrop(R_out[, , i + 1, drop = FALSE], drop = 3),
      incid$local, lambda, priors,
      t_min = t_min, t_max = t_max)

  }

  # remove burnin and thin
  keep <- seq(mcmc_control$burnin, mcmc_control$n_iter, mcmc_control$thin)
  epsilon_out <- epsilon_out[, keep, drop = FALSE]
  R_out <- R_out[, , keep, drop = FALSE]
  ## IF we have not re-ordered, we don't need to
  ## divide. Caution: this will only work for
  ## 2 variants at the moment.
  if (max_transmiss != 1) {
    epsilon_out[1 , ] <-  1 / epsilon_out[1 , ]
    if (nrow(epsilon_out) > 1) {
      for (row in 2:nrow(epsilon_out)) {
        epsilon_out[row, ] <- epsilon_out[row, ] *  epsilon_out[1, ]
      }
    }
  }
  ## TODO: Very importamt - do we need to fix
  ## ordering of R as well i.e. we have reshuffled the incidence
  ## but R should be returned in the order of the
  ## original incidence?
  list(epsilon = epsilon_out, R = R_out)
  # Not sure if this will be the same for >2 variants
}

#' @param incid a multidimensional array containing values of the incidence
#'   for each time step (1st dimension), location (2nd dimension) and
#'   pathogen/strain/variant (3rd dimension)
#'
#' @param incid_imported an optional multidimensional array containing values
#'   of the incidence of imported cases
#'   for each time step (1st dimension), location (2nd dimension) and
#'   pathogen/strain/variant (3rd dimension). `incid - incid_imported` is
#'   therefore the incidence of locally infected cases. If `incid_imported` is
#'   NULL this means there are no
#'   known imported cases and all cases other than on those from the first
#'   time step will be considered locally infected.
#'
#' @return a list with two elements.
#'   1) `local` a multidimensional array containing values of the incidence
#'   of locally infected cases
#'   for each time step (1st dimension), location (2nd dimension) and
#'   pathogen/strain/variant (3rd dimension)
#'   2) `imported` a multidimensional array containing values of the incidence
#'   of imported cases
#'   for each time step (1st dimension), location (2nd dimension) and
#'   pathogen/strain/variant (3rd dimension)
#'
#' @export
#'
#' @examples
#' n_v <- 3 # 3 variants
#' n_loc <- 1 # 1 location
#' T <- 100 # 100 time steps
#' # constant incidence 10 per day everywhere
#' incid <- array(10, dim = c(T, n_loc, n_v))
#' process_I_multivariant(incid)
#'
process_I_multivariant <- function(incid, incid_imported = NULL) {
  if (is.null(incid_imported)) {
    incid_imported <- incid
    ## only cases at first time step are imported
    incid_imported[-1, , ] <- 0
  }
  dim1 <- dim(incid)
  dim2 <- dim(incid_imported)
  if (length(dim1) != length(dim2) || !all(dim1 == dim2)) {
    stop("'incid' and 'incid_imported' have incompatible dimensions")
  }
  incid_local <- incid - incid_imported
  res <- list(local = incid_local, imported = incid_imported)
  class(res) <- "incid_multivariant"
  res
}

## TODO: check dimensions of objects is correct everywhere
## TODO: fix number of variants to be 2

