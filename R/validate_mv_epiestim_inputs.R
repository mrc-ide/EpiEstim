##' Check incidence input for MV-EpiEstim
##'
##' Check that incid is a 3-dimensional array with non-negatuve enteries.
##' 
##' @inheritParams estimate_advantage
##' @returns Silently returns TRUE if the checks are passed, otherwise throws an
##' error.
##' @author Sangeeta Bhatia
##' @internal
check_incidence <- function(incid) {
  if (any(incid < 0)) {
    stop("incid must be >=0")
  }
  if (!is.array(incid) || length(dim(incid)) != 3) {
    stop("incid must be a 3-dimensional array with dimensions time, location, variant")
  }
  invisible(TRUE)
}


##' Check priors  for MV-EpiEstim
##'
##' Check that priors is a list of the correct format and that the mean of the
##' prior for epsilon is 1 (as currently only priors with mean of epsilon equal
##' to 1 is supported).
##' 
##' @inheritParams estimate_advantage
##' @returns Silently returns TRUE if the checks are passed, otherwise throws an
##' error.
##' @author Sangeeta Bhatia
check_priors <- function(priors) {
  
  if (!identical(priors, default_priors())) {
    warning("Priors where the mean of epsilon is different from 1 are not currently supported.")
  }
  invisible(TRUE)
}

##' Check MCMC control parameters for MV-EpiEstim
##'
##' Check that (1) mcmc_control is a list of the correct format, (2) n_iter, burnin
##' and thin are positive integers and (3) n_iter is greater than burnin + thin.
##' 
##' @inheritParams estimate_advantage
##' @returns Silently returns TRUE if the checks are passed, otherwise throws an
##' error.
##' @author Sangeeta Bhatia
##' @internal
check_mcmc_control <- function(mcmc_control) {

  if (mcmc_control$n_iter < 0 || !is.integer(mcmc_control$n_iter)) {
    stop("n_iter in mcmc_control must be a positive integer")
  }
  if (mcmc_control$burnin < 0 || !is.integer(mcmc_control$burnin)) {
    stop("burnin in mcmc_control must be a positive integer")
  }
  if (mcmc_control$thin < 0 || !is.integer(mcmc_control$thin)) {
    stop("thin in mcmc_control must be a positive integer")
  }
  if (mcmc_control$n_iter < mcmc_control$burnin + mcmc_control$thin) {
    stop("In mcmc_control, n_iter must be greater than burnin + thin")
  }
  invisible(TRUE)
}

##' Validate t_min and t_max inputs for MV-EpiEstim
##'
##' Check that (1) t_min and t_max are integers, (2) that they are >= 2, that 
##' they are <= nrow(incid) and (3) t_min is not greater than t_max.
##' 
##' @inheritParams estimate_advantage
##' @returns Silently returns TRUE if the checks are passed, otherwise throws an
##' error.
##' @author Sangeeta Bhatia
##' @internal
check_t_min_t_max <- function(t_min, t_max, incid) {
    if (!is.integer(t_min) || !is.integer(t_max)) {
      stop("t_min and t_max must be integers")
    }
    if (t_min < 2 || t_max < 2) {
      stop("t_min and t_max must be >=2")
    }
    if (t_min > nrow(incid) || t_max > nrow(incid)) {
      stop("t_min and t_max must be <= nrow(incid)")
    }
    if (t_min > t_max) {
      stop("t_min is greater than t_max. You can specify a smaller t_min or increase t_max.")
    }
    invisible(TRUE)
}


##' Validate seed input for MV-EpiEstim
##'
##' Validate that the seed is not null and is numeric
##' 
##' @inheritParams estimate_advantage
##' @return Silently returns TRUE if the checks are passed, otherwise throws an error.
##' @author Sangeeta Bhatia
check_seed <- function(seed) {
  if (!is.null(seed) && !is.numeric(seed)) {
    stop("supplied seed is not a valid integer")
  }
  if (!is.null(seed)) set.seed(seed)
  invisible(TRUE)
}


check_estimate_advantage_inputs <- function(...) {
  estimate_advantage_args <- list(...)
  arg_names <- names(estimate_advantage_args)
  if ("incid" %in% arg_names) 
    check_incidence(estimate_advantage_args$incid)

  if ("si" %in% arg_names)
    check_si_distr(estimate_advantage_args$si_distr)

  if ("priors" %in% arg_names)
    check_priors(estimate_advantage_args$priors)
  
  if ("mcmc_control" %in% arg_names)
    check_mcmc_control(estimate_advantage_args$mcmc_control)
  
  if (all(c("t_min", "t_max", "incid") %in% arg_names)) {
    check_t_min_t_max(
      estimate_advantage_args$t_min, estimate_advantage_args$t_max,
      estimate_advantage_args$incid
    )
  }
  if ("seed" %in% arg_names)
    check_seed(estimate_advantage_args$seed)
}
