## Internal functions used by estimate_R with method "si_from_data".
##
## The serial interval distribution is fitted to the si_data by maximum
## likelihood using primarycensored::fitdistdoublecens (Abbott et al.), which
## accounts for double interval censoring and, if an OT column is given,
## right truncation (see Charniga et al. PLoS Comp Biol 2024). A sample of
## parameter sets is then drawn from the asymptotic normal distribution of the
## estimates and each is discretised using discr_si.

## Details of a serial interval distribution given its EpiEstim name
si_fit_distr <- function(distr) {
  distr <- convert_distr_name_for_mcmc(distr)
  fam <- sub("^off1", "", distr)
  list(
    name = distr,
    shift = as.numeric(startsWith(distr, "off1")),
    distr = switch(fam, G = "gamma", W = "weibull", L = "lnorm"),
    pdist = switch(fam,
      G = stats::pgamma, W = stats::pweibull, L = stats::plnorm
    ),
    params = switch(fam,
      G = c("shape", "scale"), W = c("shape", "scale"),
      L = c("meanlog", "sdlog")
    )
  )
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
## primarycensored::fitdistdoublecens.
## Dates are daily, so a date known exactly (EL = ER or SL = SR) is read as a
## one day interval.
si_data_to_censdata <- function(si_data, shift) {
  one_day <- si_data$ER == si_data$EL | si_data$SR == si_data$SL
  if (any(one_day)) {
    message(
      sum(one_day), " entries of si_data have EL = ER or SL = SR. ",
      "Since dates are daily these are read as one day intervals ",
      "(ER = EL + 1 or SR = SL + 1)."
    )
  }
  ER <- pmax(si_data$ER, si_data$EL + 1)
  SR <- pmax(si_data$SR, si_data$SL + 1)
  upper <- rep(Inf, nrow(si_data))
  if ("OT" %in% names(si_data)) {
    upper <- ifelse(
      is.na(si_data$OT), Inf, si_data$OT + 1 - si_data$EL - shift
    )
  }
  data.frame(
    left = si_data$SL - si_data$EL - shift,
    right = SR - si_data$EL - shift,
    pwindow = ER - si_data$EL,
    D = upper
  )
}

## Draw n parameter sets from the asymptotic normal distribution of the
## maximum likelihood estimates. Positive parameters are sampled on the log
## scale (using the delta method) so that draws stay within their support.
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
  if (!identical(mcmc_control$burnin, default_control$burnin) ||
      !identical(mcmc_control$thin, default_control$thin)) {
    warning("burnin and thin in mcmc_control are ignored as the serial ",
            "interval is estimated by maximum likelihood.", call. = FALSE)
  }

  censdata <- si_data_to_censdata(si_data, fit_distr$shift)
  start_pars <- stats::setNames(
    as.list(unname(mcmc_control$init_pars)), fit_distr$params
  )
  fit <- primarycensored::fitdistdoublecens(
    censdata,
    distr = fit_distr$distr,
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
    seed = mcmc_control$seed
  )$si_sample

  list(si_sample = si_sample, converged = converged)
}
