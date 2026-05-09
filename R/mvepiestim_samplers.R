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
#'   the distribution of the serial interval using function [compute_lambda()]
#'
#' @param priors a list of prior parameters (shape and scale of a gamma
#'   distribution) for epsilon and R; can be obtained from the function
#'   [default_priors()]. The prior for R is assumed to be the same for all
#'   time steps and all locations
#'
#' @param shape_epsilon a value or vector of values of the shape of the posterior
#'   distribution of epsilon for each of the non-reference variants, as returned
#'   by function [get_shape_epsilon()]
#'
#' @param t_min an integer > 1 giving the minimum time step to consider in the
#'   estimation. Default value is 2 (as the estimation is conditional on
#'   observations at time step 1 and can therefore only start at time step 2).
#'
#' @param t_max an integer > `t_min` and <= `nrow(incid)` giving the maximum time
#'   step to consider in the estimation. Default value is `nrow(incid)`.
#'
#' @param seed a numeric value used to fix the random seed
#'
#' @return A value or vector of values for epsilon for each non reference
#'   pathogen/strain/variant, drawn from the marginal posterior distribution
#'
#' @export
#'
#' @examples
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
draw_epsilon <- function(R, incid, lambda, priors,
                         shape_epsilon = NULL,
                         t_min = 2L, t_max = nrow(incid),
                         seed = NULL) {
  if (!is.integer(t_min) || !is.integer(t_max)){
    stop("t_min and t_max must be integers", call. = FALSE)
  }
  if (t_min < 2 || t_max < 2){
    stop("t_min and t_max must be >=2", call. = FALSE)
  }
  if(t_min > nrow(incid) || t_max > nrow(incid)){
    stop("t_min and t_max must be <= nrow(incid)", call. = FALSE)
  }
  if(any(R[!is.na(R)] < 0)) {
    stop("R must be >= 0", call. = FALSE)
  }
  if (!is.null(seed) && !is.numeric(seed)){
    stop("seed must be numeric", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)
  time_idx <- seq(t_min, t_max, 1)
  if (is.null(shape_epsilon)) {
    shape_epsilon <- get_shape_epsilon(incid, lambda, priors, t_min, t_max)
  }
  eps_scale <- get_scale_epsilon(R, lambda, priors, t = time_idx)
  stats::rgamma(dim(lambda)[3] - 1, shape = shape_epsilon, scale = eps_scale)
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
#'   the distribution of the serial interval using function [compute_lambda()]
#'
#' @param priors a list of prior parameters (shape and scale of a gamma
#'   distribution) for epsilon and R; can be obtained from the function
#'   [default_priors()]. The prior for R is assumed to be the same for all
#'   time steps and all locations
#'
#' @param shape_R_flat a vector of the shape of the posterior distribution of R
#'   for each time step t and each location l
#'   (stored in element `(l-1)*(t_max - t_min + 1) + t` of the vector),
#'   as obtained from function [get_shape_R_flat()].
#'
#' @param t_min an integer > 1 giving the minimum time step to consider in the
#'   estimation. Default value is 2 (as the estimation is conditional on
#'   observations at time step 1 and can therefore only start at time step 2).
#'
#' @param t_max an integer > `t_min` and <= `nrow(incid)` giving the maximum time
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
#' @examples
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
#' draw_R(epsilon, incid$local, lambda, priors, seed = 1, t_min = 2L)
draw_R <- function(epsilon, incid, lambda, priors,
                   shape_R_flat = NULL,
                   t_min = NULL, t_max = nrow(incid),
                   seed = NULL) {
  if (!is.integer(t_min) || !is.integer(t_max)){
    stop("t_min and t_max must be integers", call. = FALSE)
  }
  if (t_min < 2 || t_max < 2){
    stop("t_min and t_max must be >=2", call. = FALSE)
  }
  if(t_min > nrow(incid) || t_max > nrow(incid)){
    stop("t_min and t_max must be <= nrow(incid)", call. = FALSE)
  }
  if (any(epsilon < 0)){
    stop("epsilon must be > 0", call. = FALSE)
  }
  if (!is.null(seed) && !is.numeric(seed)){
    stop("seed must be numeric", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)
  time_idx <- seq(t_min, t_max, 1)
  if (is.null(shape_R_flat)) {
    shape_R_flat <- get_shape_R_flat(incid, priors, t_min, t_max)
  }
  r_scale <- get_scale_R(epsilon, incid, lambda, priors, t = time_idx)
  r_scale_flat <- as.numeric(r_scale)
  R_flat <- stats::rgamma(length(shape_R_flat), shape = shape_R_flat, scale = r_scale_flat)
  R_fill <- matrix(R_flat, nrow = length(time_idx), ncol = ncol(incid))
  R_out <- matrix(NA, nrow(incid), ncol(incid))
  R_out[time_idx, ] <- R_fill
  R_out
}



#' Compute Gelman-Rubin diagnostics for `estimate_advantage` epsilon chains
#'
#' Split each epsilon MCMC chain into two sub-chains and compute the
#' Gelman-Rubin potential scale reduction factor for each chain.
#'
#' @param epsilon_mat A numeric matrix of posterior samples for epsilon, where
#'   each row is a variant-specific epsilon chain and each column is one MCMC
#'   iteration after burn-in/thinning.
#' @param scale_reduction_threshold Numeric threshold used to determine
#'   convergence from the upper confidence limit of the scale reduction factor.
#'   Rows with any `Upper C.I.` value above this threshold are flagged as not
#'   converged. Defaults to `1.1`.
#'
#' @return A named list with:
#' \itemize{
#'   \item `diagnostic`: a list of `coda::gelman.diag()` outputs (one per row of
#'   `epsilon_mat`).
#'   \item `convergence`: a logical vector with one value per row of
#'   `epsilon_mat`, where `TRUE` indicates convergence under
#'   `scale_reduction_threshold`.
#' }
#'
#' @keywords internal
compute_convergence_diagnostic <- function(
    epsilon_mat, scale_reduction_threshold = 1.1) {
  diagnostic <- lapply(seq_len(nrow(epsilon_mat)), function(row) {
    x <- epsilon_mat[row, ]
    if (length(x) %% 2 != 0) {
      eps1 <- coda::as.mcmc(x[1:((length(x) + 1) / 2)])
      eps2 <- coda::as.mcmc(x[length(eps1):length(x)])
    }
    if (length(x) %% 2 == 0) {
      eps1 <- coda::as.mcmc(x[1:(length(x) / 2)])
      eps2 <- coda::as.mcmc(x[(length(eps1) + 1):length(x)])
    }

    eps_2chain <- coda::mcmc.list(eps1, eps2)
    coda::gelman.diag(eps_2chain, confidence = 0.95)
  })

  conv_check <- unlist(lapply(diagnostic, function(x) {
    !any(x$psrf[, "Upper C.I."] > scale_reduction_threshold)
  }))

  list(diagnostic = diagnostic, convergence = conv_check)
}
