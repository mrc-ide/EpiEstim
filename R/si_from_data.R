## Internal functions used by estimate_R with method "si_from_data".
##
## The serial interval distribution is fitted to the si_data by maximum
## likelihood using primarycensored::fitdistdoublecens (Abbott et al.), which
## accounts for double interval censoring and, if an OT column is given,
## right truncation (see Charniga et al. PLoS Comp Biol 2024). A sample of
## parameter sets is then drawn from the asymptotic normal distribution of the
## estimates and each is discretised using discr_si.

## Details of a serial interval distribution given its EpiEstim name, e.g.
## "gamma" or "gamma_offset_1" (see si_distributions)
si_fit_distr <- function(distr) {
  old_names <- c(G = "gamma", W = "weibull", L = "lognormal",
                 off1G = "gamma_offset_1", off1W = "weibull_offset_1",
                 off1L = "lognormal_offset_1")
  if (distr %in% names(old_names)) {
    distr <- old_names[[distr]]
  }
  shift <- as.numeric(endsWith(distr, "_offset_1"))
  entry <- si_distribution_from_alias(sub("_offset_1$", "", distr))
  if (is.null(entry)) {
    stop("Unsupported distribution name: ", distr, call. = FALSE)
  }
  c(list(epiestim_name = distr, shift = shift), entry)
}

## Arguments of discr_si used when discretising the fitted serial interval.
## dist and shift are set by si_parametric_distr.
si_from_data_discr_args <- function(config, fit_distr) {
  discr_args <- config$si_discr_args
  if (is.null(discr_args)) {
    discr_args <- list()
  }
  ## dist, shift and the lower bound of L are set by si_parametric_distr
  check_si_discr_args(
    discr_args[setdiff(names(discr_args), c("dist", "shift", "L"))]
  )
  if (!is.null(discr_args$dist) &&
        !identical(discr_args$dist, fit_distr$pdist)) {
    stop("si_discr_args$dist conflicts with si_parametric_distr. ",
         "For method si_from_data the distribution is set by ",
         "si_parametric_distr.", call. = FALSE)
  }
  if (!is.null(discr_args$shift) && discr_args$shift != fit_distr$shift) {
    stop("si_discr_args$shift conflicts with si_parametric_distr. ",
         "For method si_from_data the shift is set by si_parametric_distr ",
         "(1 for the offset distributions, 0 otherwise).", call. = FALSE)
  }
  discr_args$dist <- NULL
  discr_args$shift <- NULL
  if (is.null(discr_args$dprimary)) {
    discr_args$dprimary <- stats::dunif
  }
  if (is.null(discr_args$primary_args)) {
    discr_args$primary_args <- list()
  }
  discr_args
}

## Interval censored data in the format used by
## primarycensored::fitdistdoublecens().
## Following coarseDataTools, EL, ER, SL and SR are continuous bounds, so a
## date known to the day is EL = d, ER = d + 1, and EL = ER (or SL = SR) is an
## exactly known time. An exact primary time gives pwindow = 0 and an exact
## secondary time gives left = right, which contributes a density. OT is also
## a continuous bound: secondary onsets are observed up to time OT.
si_data_to_censdata <- function(si_data, shift) {
  upper <- rep(Inf, nrow(si_data))
  if ("OT" %in% names(si_data)) {
    upper <- ifelse(
      is.na(si_data$OT), Inf, si_data$OT - si_data$EL - shift
    )
  }
  data.frame(
    left = si_data$SL - si_data$EL - shift,
    right = si_data$SR - si_data$EL - shift,
    pwindow = si_data$ER - si_data$EL,
    D = upper
  )
}

## Draw n parameter sets from the asymptotic normal distribution of the
## maximum likelihood estimates.
## The maximum likelihood estimate theta_hat is asymptotically normal with
## covariance V, the inverse Hessian returned by fitdistrplus. For a positive
## parameter we draw log(theta) instead, which by the multivariate delta method
## is asymptotically normal with mean log(theta_hat) and covariance J V J,
## where J = diag(1 / theta_hat) for the positive parameters (and 1 for the
## others). Exponentiating keeps the draws positive, and using the full
## covariance keeps the correlation between parameters, which is strong for
## e.g. the shape and scale of a Gamma distribution. This follows
## extract_mle_draws() in R/fit-utils.R of
## https://github.com/epinowcast/primarycensored-paper.
draw_si_params <- function(estimate, vcov, n, positive) {
  mu <- ifelse(positive, log(estimate), estimate)
  jacobian <- ifelse(positive, 1 / estimate, 1)
  vcov_log <- outer(jacobian, jacobian) * vcov
  vcov_log <- (vcov_log + t(vcov_log)) / 2
  eigen_values <- tryCatch(
    eigen(vcov_log, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) NA
  )
  if (anyNA(eigen_values) ||
        min(eigen_values) <= sqrt(.Machine$double.eps)) {
    stop("The uncertainty in the serial interval parameters could not be ",
         "estimated from si_data, so no sample of serial interval ",
         "distributions can be drawn. The fitted distribution may be a poor ",
         "fit to the data. Try different starting values with ",
         "make_mcmc_control(init_pars = ...) or a different ",
         "si_parametric_distr.", call. = FALSE)
  }
  z <- matrix(stats::rnorm(n * length(mu)), nrow = n)
  draws <- sweep(z %*% chol(vcov_log), 2, mu, "+")
  draws[, positive] <- exp(draws[, positive])
  colnames(draws) <- names(estimate)
  as.data.frame(draws)
}

## Estimate the serial interval from si_data and return a sample of discrete
## serial interval distributions
si_sample_from_data <- function(si_data, config) {
  fit_distr <- si_fit_distr(config$si_parametric_distr)
  discr_args <- si_from_data_discr_args(config, fit_distr)

  mcmc_control <- config$mcmc_control
  default_control <- make_mcmc_control()
  if (!is.null(mcmc_control) &&
        (!identical(mcmc_control$burnin, default_control$burnin) ||
           !identical(mcmc_control$thin, default_control$thin))) {
    warning("burnin and thin in mcmc_control are ignored as the serial ",
            "interval is estimated by maximum likelihood.", call. = FALSE)
  }

  censdata <- si_data_to_censdata(si_data, fit_distr$shift)
  if (is.null(mcmc_control$init_pars)) {
    start_pars <- si_start_values(si_data, fit_distr$epiestim_name)
  } else {
    start_pars <- stats::setNames(
      as.list(unname(mcmc_control$init_pars)), fit_distr$params
    )
  }
  fit <- primarycensored::fitdistdoublecens(
    censdata,
    distr = fit_distr$name,
    start = start_pars,
    dprimary = discr_args$dprimary,
    primary_args = discr_args$primary_args
  )
  converged <- fit$convergence == 0
  if (!converged) {
    warning("The maximum likelihood estimation of the serial interval ",
            "did not converge.", call. = FALSE)
  }

  si_sample <- primary2estim(
    fit,
    dist = fit_distr$pdist,
    n = config$n1,
    shift = fit_distr$shift,
    si_discr_args = discr_args,
    seed = if (is.null(mcmc_control)) config$seed else mcmc_control$seed
  )$si_sample

  list(si_sample = si_sample, converged = converged)
}
