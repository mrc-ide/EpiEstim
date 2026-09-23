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
    ),
    positive = switch(fam,
      G = c(TRUE, TRUE), W = c(TRUE, TRUE), L = c(FALSE, TRUE)
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

## Discretise a sample of serial interval distributions.
## Following the approach of the distspec package, the support is set from
## the 0.999 quantile of the primary censored serial interval (largest across
## the sample) and each distribution is right truncated at the end of this
## support. For distributions without an offset P(SI = 0) is set to 0 by
## truncating below 1.
## Returns a matrix with one column per distribution, as used by
## estimate_R with method "si_from_sample".
si_sample_from_params <- function(fit_distr, samples, discr_args = list()) {
  if (is.null(discr_args$dprimary)) {
    discr_args$dprimary <- stats::dunif
  }
  if (is.null(discr_args$primary_args)) {
    discr_args$primary_args <- list()
  }
  params <- lapply(seq_len(nrow(samples)), function(i) {
    stats::setNames(as.list(unlist(samples[i, ])), fit_distr$params)
  })
  q_max <- max(vnapply(params, function(p) {
    do.call(
      primarycensored::qprimarycensored,
      c(
        list(
          p = 0.999, pdist = fit_distr$pdist, pwindow = 1,
          dprimary = discr_args$dprimary,
          primary_args = discr_args$primary_args
        ),
        p,
        list(check = FALSE)
      )
    )
  }))
  max_value <- ceiling(q_max + fit_distr$shift)
  k <- seq(0, max_value)
  lower <- if (fit_distr$shift == 0) 1 else -Inf
  if (!is.null(discr_args$L)) {
    lower <- max(lower, discr_args$L)
  }
  upper <- max_value + 1
  if (!is.null(discr_args$D)) {
    upper <- min(upper, discr_args$D)
  }
  si_args <- si_discr_defaults(list(
    dist = fit_distr$pdist, shift = fit_distr$shift, L = lower, D = upper,
    dprimary = discr_args$dprimary, primary_args = discr_args$primary_args
  ))
  t(discr_si_param_draws(k, params, si_args))
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
    warning("burnin and thin in mcmc_control are ignored. The serial ",
            "interval is now estimated by maximum likelihood rather than ",
            "MCMC.", call. = FALSE)
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

  if (!is.null(mcmc_control$seed)) {
    set.seed(mcmc_control$seed)
  }
  samples <- draw_si_params(
    fit$estimate[fit_distr$params],
    fit$vcov[fit_distr$params, fit_distr$params],
    config$n1,
    fit_distr$positive
  )
  si_sample <- si_sample_from_params(fit_distr, samples, discr_args)

  list(si_sample = si_sample, converged = converged)
}
