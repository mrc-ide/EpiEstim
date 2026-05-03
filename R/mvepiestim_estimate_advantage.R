#' Estimate instantaneous reproduction number
#' 
#' Jointly estimate the instantaneous reproduction number for a reference
#' pathogen/strain/variant and the relative transmissibility of a "new" 
#' pathogen/strain/variant.
#'
#' @param incid a multidimensional array containing values of the incidence
#'   for each time step (1st dimension), location (2nd dimension) and
#'   pathogen/strain/variant (3rd dimension)
#'
#' @param si_distr a matrix with two columns, each containing the probability mass
#'   function for the discrete serial interval for each of the two
#'   pathogen/strain/variants, starting with the probability mass function
#'   for day 0 in the first row, which should be 0. Each column in the matrix
#'   should sum to 1.
#'
#' @param priors a list of prior parameters (shape and scale of a gamma
#'   distribution) for epsilon and R; can be obtained from the function
#'   [default_priors()]. The prior for R is assumed to be the same for all
#'   time steps and all locations
#'
#' @param mcmc_control a list of default MCMC control parameters, obtained by 
#'   default from [default_mcmc_controls()]
#'
#' @param t_min an integer > 1 giving the minimum time step to consider in the
#'   estimation.
#'   If `NULL`, `t_min` is calculated using the function [compute_si_cutoff()]
#'   which gets the maximum (across all variants) of the 95th percentile of the
#'   SI distribution.
#'
#' @param t_max an integer > `t_min` and <= `nrow(incid)` giving the maximum time
#'   step to consider in the estimation. Default value is `nrow(incid)`.
#'
#' @param seed a numeric value used to fix the random seed
#'
#' @param incid_imported an optional multidimensional array containing values
#'   of the incidence of imported cases for each time step (1st dimension), 
#'   location (2nd dimension) and pathogen/strain/variant (3rd dimension). 
#'   `incid - incid_imported` is therefore the incidence of locally infected 
#'   cases. If `incid_imported` is `NULL` this means there are no
#'   known imported cases and all cases other than on those from the first
#'   time step will be considered locally infected.
#'
#' @param precompute a boolean (defaulting to `TRUE`) deciding whether to
#'   precompute quantities or not. Using `TRUE` will make the algorithm faster
#'   
#' @param reorder_incid a boolean (defaulting to `TRUE`) deciding whether the
#'   incidence array can be internally reordered during the estimation of the
#'   transmission advantage. If `TRUE`, the most transmissible pathogen/strain/variant
#'   is temporarily assigned to `[, , 1]` of the incidence array. We recommend the
#'   default value of `TRUE` as we find this to stabilise inference.
#'
#' @return A list with the following elements:
#' - `epsilon`: a matrix containing the MCMC chain (thinned and after burnin)
#'   for the relative transmissibility of the "new" pathogen/strain/variant(s)
#'   compared to the reference. Each row corresponds to a "new"
#'   pathogen/strain/variant and each column to an MCMC iteration.
#' - `R`: an array containing the MCMC chain (thinned and after burnin) for the
#'   reproduction number for the reference pathogen/strain/variant. Dimensions
#'   of the array are time,  location, and iteration of the MCMC.
#' - `convergence`: a logical vector based on the Gelman-Rubin convergence
#'   diagnostic. Each element in `convergence` is `TRUE` when the MCMC for the
#'   corresponding epsilon has converged number of iterations specified 
#'   (otherwise `FALSE`).
#' - `diag`: a nested list of the point estimate and upper confidence limits of
#'   the Gelman-Rubin convergence diagnostics (as implemented in coda). The
#'   length of `diag` is equal to the number of rows in epsilon. Each element od
#'   `diag` is a list of length 2 where the first element is called `psrf` and
#'   is a named list of the point estimate and upper confidence limits. The
#'   second element is `NULL` and can be ignored.
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
#' # arbitrary serial interval, same for both variants
#' w_v <- c(0, 0.2, 0.5, 0.3)
#' si_distr <- cbind(w_v, w_v)
#'
#' # Dummy initial values for the MCMC
#' R_init <- matrix(5, nrow = T, ncol = n_loc)
#' R_init[1, ] <- NA # no estimates of R on first time step
#' epsilon_init <- 5
#' x <- estimate_advantage(incid, si_distr, priors)
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

estimate_advantage <- function(incid, si_distr, priors = default_priors(),
                           mcmc_control = default_mcmc_controls(),
                           t_min = NULL, t_max = nrow(incid),
                           seed = NULL,
                           incid_imported = NULL,
                           precompute = TRUE,
                           reorder_incid = TRUE) {

  if (is.null(t_min)) {
    t_min <- compute_t_min(incid, si_distr)
  }
  if (!is.integer(t_min) || !is.integer(t_max)) {
    stop("t_min and t_max must be integers")
  }
  if (t_min < 2 || t_max < 2){
    stop("t_min and t_max must be >=2")
  }
  if(t_min > nrow(incid) || t_max > nrow(incid)){
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
  if (mcmc_control$n_iter < 0 || !is.integer(mcmc_control$n_iter)){
    stop("n_iter in mcmc_control must be a positive integer")
  }
  if (mcmc_control$burnin < 0 || !is.integer(mcmc_control$burnin)){
    stop("burnin in mcmc_control must be a positive integer")
  }
  if (mcmc_control$thin < 0 || !is.integer(mcmc_control$thin)){
    stop("thin in mcmc_control must be a positive integer")
  }
  if (mcmc_control$n_iter < mcmc_control$burnin + mcmc_control$thin){
    stop("In mcmc_control, n_iter must be greater than burnin + thin")
  }
  if (!is.null(seed) && !is.numeric(seed)){
    stop("seed must be numeric")
  }
  if (!is.null(seed)) set.seed(seed)

  if (t_min > t_max) {
    stop("t_min is greater than t_max. You can specify a smaller t_min or increase t_max.")
  }
  
  if (!identical(priors, default_priors())) {
    warning("Priors where the mean of epsilon is different from 1 are not currently supported.")
  }

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

  if (reorder_incid) {
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
  }

  ## Precalculate quantities of interest
  if (precompute) {
    shape_R_flat <- get_shape_R_flat(incid$local, priors, t_min, t_max)
    shape_epsilon <- get_shape_epsilon(incid$local, lambda, priors, t_min, t_max)
  } else {
    shape_R_flat <- NULL
    shape_epsilon <- NULL
  }

  epsilon_init <- unlist(lapply(seq(2, length(R_init)), function(i)
    stats::median(R_init[[i]] / R_init[[1]], na.rm = TRUE)))
  epsilon_out <- matrix(NA, nrow = length(epsilon_init),
                        ncol = mcmc_control$n_iter + 1)
  epsilon_out[, 1] <- epsilon_init
  
  R_init <- draw_R(
    epsilon = epsilon_init, incid = incid$local, lambda = lambda,
    priors = priors, shape_R_flat = shape_R_flat, t_min = t_min, t_max = t_max
  )
  
  R_out <- array(NA, dim= c(T, n_loc, mcmc_control$n_iter + 1))
  R_out[, , 1] <- R_init
  for (i in seq_len(mcmc_control$n_iter)) {
    R_out[, , i + 1] <- draw_R(epsilon_out[, i], incid$local, lambda, priors,
                               shape_R_flat = shape_R_flat,
                               t_min = t_min, t_max = t_max)
    epsilon_out[, i + 1] <- draw_epsilon(
      abind::adrop(R_out[, , i + 1, drop = FALSE], drop = 3),
      incid$local, lambda, priors,
      shape_epsilon = shape_epsilon,
      t_min = t_min, t_max = t_max)

  }

  # remove burnin and thin
  keep <- seq(mcmc_control$burnin, mcmc_control$n_iter, mcmc_control$thin)
  epsilon_out <- epsilon_out[, keep, drop = FALSE]
  R_out <- R_out[, , keep, drop = FALSE]
  
  ## IF we have not re-ordered, we don't need to
  ## divide. Caution: this will only work for
  ## 2 variants at the moment.
  if (reorder_incid && max_transmiss != 1) {
    epsilon_out[1 , ] <-  1 / epsilon_out[1 , ]
    if (nrow(epsilon_out) > 1) {
      for (row in 2:nrow(epsilon_out)) {
        epsilon_out[row, ] <- epsilon_out[row, ] *  epsilon_out[1, ]
      }
    }
    R_out <- R_out / as.vector(epsilon_out)
  }
  
  # Add in convergence check (gelman diagnostic)
  # Split epsilon into 2 chains
  diag <- lapply(
    seq_len(nrow(epsilon_out)), function(row) {
      x <- epsilon_out[row, ]
      if(length(x)%%2 != 0){ # if length is odd
        eps1 <- coda::as.mcmc(x[1:((length(x)+1)/2)])
        eps2 <- coda::as.mcmc(x[length(eps1):length(x)])
      }
      if(length(x)%%2==0){ # if length is even
        eps1 <- coda::as.mcmc(x[1:(length(x)/2)])
        eps2 <- coda::as.mcmc(x[(length(eps1)+1):length(x)])
      }

      eps_2chain <- coda::mcmc.list(eps1,eps2)
      coda::gelman.diag(eps_2chain, confidence = 0.95)

    }
  )
  conv_check <- lapply(diag, function(x) {
  # Are any of the scale reduction factors >1.1?
   if (any(x$psrf[, "Upper C.I."] > 1.1)) {
     message("The Gelman-Rubin algorithm suggests the MCMC may not have converged
                  within the number of iterations specified.")
     convergence <- FALSE
   } else {
     convergence <- TRUE
   }
    convergence
  }
  )
  conv_check <- unlist(conv_check)

  list(epsilon = epsilon_out, R = R_out, convergence = conv_check, diag = diag #, max_transmiss = max_transmiss
       ) 
  
}

