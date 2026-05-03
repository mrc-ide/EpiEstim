#' Precompute shape of posterior distribution for R
#'
#' @param incid a multidimensional array containing values of the (local)
#'   incidence
#'   for each time step (1st dimension), location (2nd dimension) and
#'   pathogen/strain/variant (3rd dimension)
#'
#' @param priors a list of prior parameters (shape and scale of a gamma
#'   distribution) for epsilon and R; can be obtained from the function
#'   `default_priors`. The prior for R is assumed to be the same for all
#'   time steps and all locations
#'
#' @param t_min an integer > 1 giving the minimum time step to consider in the
#'   estimation. Default value is 2 (as the estimation is conditional on
#'   observations at time step 1 and can therefore only start at time step 2).
#'
#' @param t_max an integer > `t_min` and <= `nrow(incid)` giving the maximum time
#'   step to consider in the estimation. Default value is `nrow(incid)`.
#'
#' @return a vector of the shape of the posterior distribution of R for
#'   each time step t and each location l
#'   (stored in element `(l-1)*(t_max - t_min + 1) + t` of the vector)
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
#' get_shape_R_flat(incid, priors)

get_shape_R_flat <- function(incid, priors, t_min = 2L, t_max = nrow(incid)) {
  t <- seq(t_min, t_max, 1)
  shape <- apply(incid[t, , , drop = FALSE], c(1, 2), sum) + priors$R$shape
  as.numeric(shape)
}

#' Precompute shape of posterior distribution for epsilon
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
#' @param t_min an integer > 1 giving the minimum time step to consider in the
#'   estimation. Default value is 2 (as the estimation is conditional on
#'   observations at time step 1 and can therefore only start at time step 2).
#'
#' @param t_max an integer > `t_min` and <= `nrow(incid)` giving the maximum time
#'   step to consider in the estimation. Default value is `nrow(incid)`.
#'
#' @return a value or vector of values of the shape of the posterior
#'   distribution of epsilon for each of the non-reference variants
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
#' get_shape_epsilon(incid$local, lambda, priors)
get_shape_epsilon <- function(incid, lambda, priors,
                              t_min = 2L, t_max = nrow(incid)) {
  t <- seq(t_min, t_max, 1)
  vnapply(seq(2, dim(lambda)[3]), function(e)
    sum(incid[t, , e])) + priors$epsilon$shape
}


get_scale_epsilon <- function(R, lambda, priors, t) {
  rate <- vnapply(seq(2, dim(lambda)[3]), function(e) {
    sum(R[t, ] * lambda[t, , e]) + 1 / priors$epsilon$scale
  })
  1 / rate
}

get_scale_R <- function(epsilon, incid, lambda, priors, t) {
  temp <- lambda[t, , 1]
  idx <- seq(2, dim(incid)[3], 1)
  for (var in idx) {
    temp <- temp + epsilon[var - 1] * lambda[t, , var]
  }
  rate <- temp + 1 / priors$R$scale
  1 / rate
}


#' Compute the overall infectivity
#'
#' @param incid a list (as obtained from function [process_I_multivariant()])
#'   of two multidimensional arrays (`local` and `imported`) containing values 
#'   of the incidence
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
#'   algorithm much faster.
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

compute_lambda <- function(incid, si_distr) {
  if (!inherits(incid, "incid_multivariant")) {
    msg1 <- "'incid 'should be an 'incid_multivariant' object. "
    msg2 <- "Use function 'process_I_multivariant' first"
    stop(msg1, msg2)
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



#' Index before which at most a given probability
#' mass is captured
#'
#' Across a matrix of discretised probability distributions
#' (see [estimate_advantage()]) this function returns the largest index
#' (across all columns) such that the cumulative probability mass before index is
#' `1 - miss_at_most`.
#'
#' @inheritParams estimate_advantage
#' @param miss_at_most numeric. Probability mass in the tail of the SI distribution
#' 
#' @return integer
#' @author Sangeeta Bhatia
#' @export

compute_si_cutoff <- function(si_distr, miss_at_most = 0.05) {
  if (any(colSums(si_distr) != 1)) {
    warning("Input SI distributions should sum to 1. Normalising now")
    si_distr <- si_distr / colSums(si_distr)
  }
  cutoff <- 1 - miss_at_most
  cdf <- apply(si_distr, 2, cumsum)
  idx <- apply(
    cdf, 2,
    function(col) Position(function(x) x > cutoff, col)
  )
  as.integer(max(idx))
}

#' First day of non-zero incidence
#' 
#' Get the first day of non-zero incidence across all variants and locations.
#' 
#' For each variant, find the first day of non-zero incidence. The maximum of 
#' these is the smallest possible point at which estimation can begin.
#' 
#' @inheritParams estimate_advantage
#' 
#' @return integer
#' @author Sangeeta Bhatia
#' 
#' @export

first_nonzero_incid <- function(incid) {
  t_min_incid <- apply(
    incid, c(2, 3),
    function(vec) Position(function(x) x > 0, vec)
  )
  if (anyNA(t_min_incid)) {
    warning(
      "For some variants/locations, incidence is
       always zero. This will cause estimate_advantage to fail."
    )
  }
  max(t_min_incid)
}

#' Compute the smallest index at which joint estimation should start
#'
#' Unless specified by the user, `t_min` in [estimate_advantage()] is computed 
#' as the sum of two indices:
#' 
#' - the first day of non-zero incidence across all locations, computed using 
#'   [first_nonzero_incid()]
#' - the 95th percentile of the probability mass function of the SI distribution 
#'   across all variants computed using [compute_si_cutoff()]
#'
#' @inheritParams compute_si_cutoff
#' @inheritParams estimate_advantage
#' 
#' @return integer
#' @author Sangeeta Bhatia
#' 
#' @export

compute_t_min <- function(incid, si_distr, miss_at_most) {
  t_min_si <- compute_si_cutoff(si_distr, 0.05)
  t_min_incid <- first_nonzero_incid(incid)
  as.integer(t_min_incid + t_min_si)
}
