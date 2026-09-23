#' Compute discretized generation time distribution
#'
#' Compute the discrete distribution of the serial interval,
#' assuming by default that the serial interval is a shifted Gamma
#' distributed, with shift 1.
#'
#' Assuming that the serial interval is a shifted Gamma distributed with mean
#' \eqn{\mu}, standard deviation \eqn{\sigma} and shift \eqn{1},
#' the discrete probability \eqn{w_k} that the serial interval is equal to
#' \eqn{k} is (see supplement of Cori et al. AJE 2013):
#'
#' \deqn{w_k = kF_{\{\mu-1,\sigma\}}(k)+(k-2)F_{\{\mu-1,\sigma\}}
#' (k-2)-2(k-1)F_{\{\mu-1,\sigma\}}(k-1)\\
#' +(\mu-1)(2F_{\{\mu-1+\frac{\sigma^2}{\mu-1},
#' \sigma\sqrt{1+\frac{\sigma^2}{\mu-1}}\}}(k-1)-
#' F_{\{\mu-1+\frac{\sigma^2}{\mu-1},
#' \sigma\sqrt{1+\frac{\sigma^2}{\mu-1}}\}}(k-2)-
#' F_{\{\mu-1+\frac{\sigma^2}{\mu-1},
#' \sigma\sqrt{1+\frac{\sigma^2}{\mu-1}}\}}(k))}
#'
#' where \eqn{F_{\{\mu,\sigma\}}} is the cumulative density function of a Gamma
#' distribution with mean \eqn{\mu} and standard deviation \eqn{\sigma}.
#'
#' This is the probability that a continuous delay, starting at a time
#' uniformly distributed within a day and shifted by `shift` days (one by
#' default), ends on day \eqn{k}. It is computed using the
#' [primarycensored][primarycensored::primarycensored-package] package
#' (Abbott et al.), which also allows other delay distributions (`dist`), a
#' different `shift`, other distributions of the primary event time within
#' the day (`dprimary`), and truncation of the serial interval below `L` and at or
#' above `D`. When `L` or `D` are finite, \eqn{w_k} is normalised over
#' \eqn{L \le k < D}.
#'
#' @param k Positive integer, or vector of positive integers for which the
#' discrete distribution is desired.
#' @param mu A positive real giving the mean of the Gamma distribution.
#' Must be greater than `shift` (1 by default).
#' @param sigma A non-negative real giving the standard deviation of the Gamma
#' distribution.
#' @param dist The distribution of the serial interval. Either one of
#' "gamma" (the default), "lognormal" or "weibull", parameterised by `mu` and
#' `sigma`, or a cumulative distribution function (e.g. [stats::pgamma()])
#' whose parameters are passed through `...`. In the latter case `mu` and
#' `sigma` should not be specified, and the distribution applies to the
#' serial interval minus `shift`.
#' @param shift A non-negative real giving the shift of the serial interval
#' distribution. The default of 1 means that the serial interval is at
#' least one day.
#' @param L The serial interval is truncated below `L`, so that
#' \eqn{w_k = 0} for \eqn{k < L}. Defaults to `-Inf`, meaning no truncation.
#' @param D The serial interval is truncated at `D`, so that
#' \eqn{w_k = 0} for \eqn{k \ge D}. Defaults to `Inf`, meaning no truncation.
#' @param dprimary The probability density function of the time of the
#' primary event within the day. Defaults to [stats::dunif()]. See
#' [primarycensored::dprimarycensored()] for other options, such as
#' [primarycensored::dexpgrowth()].
#' @param primary_args A list of additional arguments passed to `dprimary`.
#' @param ... Parameters of `dist` when `dist` is a function.
#'
#' @return Gives the discrete probability \eqn{w_k} that the serial interval is
#' equal to \eqn{k}. This is not normalised over `k`, so it sums to less than 1
#' if `k` does not cover the support of the distribution. Use `D` to truncate
#' and normalise the distribution.
#'
#' @seealso [overall_infectivity()], [estimate_R()],
#' [primarycensored::dprimarycensored()]
#'
#' @author Anne Cori
#'
#' @references
#' Cori, A. et al. A new framework and software to estimate time-varying
#' reproduction numbers during epidemics (AJE 2013).
#'
#' Abbott, S. et al. [primarycensored][primarycensored::primarycensored-package]:
#' Primary Event Censored Distributions.
#' \doi{10.5281/zenodo.13632839}
#'
#' @export
#'
#' @examples
#' ## Computing the discrete serial interval of influenza
#' mean_flu_si <- 2.6
#' sd_flu_si <- 1.5
#' dicrete_si_distr <- discr_si(seq(0, 20), mean_flu_si, sd_flu_si)
#' plot(seq(0, 20), dicrete_si_distr, type = "h",
#'      lwd = 10, lend = 1, xlab = "time (days)", ylab = "frequency")
#' title(main = "Discrete distribution of the serial interval of influenza")
#'
#' ## Using a lognormal serial interval with no shift, truncated below 1
#' discr_si(seq(0, 20), mean_flu_si, sd_flu_si, dist = "lognormal",
#'          shift = 0, L = 1)

discr_si <- function(k, mu, sigma, dist = "gamma", shift = 1, L = -Inf,
                     D = Inf, dprimary = stats::dunif, primary_args = list(),
                     ...)
{
  if (is.function(dist)) {
    if (!missing(mu) || !missing(sigma)) {
      stop("mu and sigma should not be specified when dist is a function; ",
           "pass the parameters of dist through ... instead.", call. = FALSE)
    }
    pdist <- dist
    dist_args <- list(...)
  } else {
    check_si_moments(mu, sigma, shift)
    pdist <- si_pdist(dist)
    dist_args <- si_dist_args(dist, mu - shift, sigma)
  }
  check_si_support(k, shift, L, D)

  res <- numeric(length(k))
  in_support <- k >= L & k < D
  if (any(in_support)) {
    res[in_support] <- do.call(
      primarycensored::dprimarycensored,
      c(
        list(
          x = k[in_support] - shift, pdist = pdist, pwindow = 1,
          L = L - shift, D = D - shift, dprimary = dprimary,
          primary_args = primary_args
        ),
        dist_args,
        list(check = FALSE)
      )
    )
  }
  res <- vnapply(res, function(e) max(0, e))

  return(res)
}

## Call discr_si with the extra arguments given in config$si_discr_args
discr_si_config <- function(k, mu, sigma, si_discr_args = NULL) {
  if (is.null(si_discr_args)) {
    si_discr_args <- list()
  }
  check_si_discr_args(si_discr_args)
  do.call(discr_si, c(list(k = k, mu = mu, sigma = sigma), si_discr_args))
}

## Discretise the serial interval for several draws of its mean and standard
## deviation. Gives the same result as calling discr_si for each draw but
## builds the primarycensored object once and only updates its parameters.
## Returns a matrix with one row per draw.
discr_si_draws <- function(k, mu, sigma, si_discr_args = NULL) {
  if (is.null(si_discr_args)) {
    si_discr_args <- list()
  }
  check_si_discr_args(si_discr_args)
  si_args <- utils::modifyList(
    list(
      dist = "gamma", shift = 1, L = -Inf, D = Inf, dprimary = stats::dunif,
      primary_args = list()
    ),
    si_discr_args
  )
  check_si_moments(mu, sigma, si_args$shift)
  check_si_support(k, si_args$shift, si_args$L, si_args$D)

  ## work on the delay scale, i.e. the serial interval minus the shift
  delay_L <- si_args$L - si_args$shift
  delay_D <- si_args$D - si_args$shift
  in_support <- k >= si_args$L & k < si_args$D
  lower <- k[in_support] - si_args$shift
  upper <- pmin(lower + 1, delay_D)
  cdf_points <- sort(unique(c(
    lower, upper, delay_L[is.finite(delay_L)], delay_D[is.finite(delay_D)]
  )))
  pos_cdf_points <- cdf_points[cdf_points > 0 & is.finite(cdf_points)]

  pdist <- si_pdist(si_args$dist)
  pcens <- do.call(
    primarycensored::new_pcens,
    c(
      list(pdist = pdist, dprimary = si_args$dprimary,
           primary_args = si_args$primary_args),
      si_dist_args(si_args$dist, mu[1] - si_args$shift, sigma[1])
    )
  )

  res <- matrix(0, nrow = length(mu), ncol = length(k))
  for (i in seq_along(mu)) {
    pcens$args <- si_dist_args(si_args$dist, mu[i] - si_args$shift, sigma[i])
    ## the delay is non-negative so the censored CDF is 0 at and below 0
    cdf <- numeric(length(cdf_points))
    cdf[cdf_points == Inf] <- 1
    if (length(pos_cdf_points) > 0) {
      cdf[cdf_points > 0 & is.finite(cdf_points)] <-
        primarycensored::pcens_cdf(pcens, pos_cdf_points, pwindow = 1)
    }
    cdf_L <- if (is.finite(delay_L)) cdf[match(delay_L, cdf_points)] else 0
    cdf_D <- if (is.finite(delay_D)) cdf[match(delay_D, cdf_points)] else 1
    pmf <- (cdf[match(upper, cdf_points)] - cdf[match(lower, cdf_points)]) /
      (cdf_D - cdf_L)
    res[i, in_support] <- pmax(pmf, 0)
  }

  return(res)
}

## Check the extra arguments to discr_si given in config$si_discr_args
check_si_discr_args <- function(si_discr_args) {
  allowed <- c("dist", "shift", "L", "D", "dprimary", "primary_args")
  if (!is.list(si_discr_args) ||
        (length(si_discr_args) > 0 && is.null(names(si_discr_args)))) {
    stop("si_discr_args must be a named list.", call. = FALSE)
  }
  unknown <- setdiff(names(si_discr_args), allowed)
  if (length(unknown) > 0) {
    stop("si_discr_args contains unsupported arguments: ",
         toString(unknown), ". Supported arguments are: ",
         toString(allowed), ".", call. = FALSE)
  }
  if (!is.null(si_discr_args$dist) && !is.character(si_discr_args$dist)) {
    stop("si_discr_args$dist must be one of 'gamma', 'lognormal' or ",
         "'weibull'.", call. = FALSE)
  }
  support <- utils::modifyList(list(shift = 1, L = -Inf), si_discr_args)
  if (support$shift < 1 && support$L < 1) {
    stop("si_discr_args gives a non-zero probability of a serial interval ",
         "of zero, which EpiEstim does not allow. Use shift >= 1 or L >= 1.", call. = FALSE)
  }
  invisible(NULL)
}

## Cumulative distribution function for a named serial interval distribution
si_pdist <- function(dist) {
  pdists <- list(
    gamma = stats::pgamma, lognormal = stats::plnorm, weibull = stats::pweibull
  )
  if (!is.character(dist) || length(dist) != 1 || !dist %in% names(pdists)) {
    stop("dist must be one of 'gamma', 'lognormal' or 'weibull', ",
         "or a cumulative distribution function.", call. = FALSE)
  }
  pdists[[dist]]
}

## Check the mean and standard deviation of the serial interval
check_si_moments <- function(mu, sigma, shift) {
  if (any(sigma < 0)) {
    stop("sigma must be >=0.", call. = FALSE)
  }
  if (any(mu <= shift)) {
    stop("mu must be >", shift, call. = FALSE)
  }
  invisible(NULL)
}

## Check the values, shift and truncation of the serial interval
check_si_support <- function(k, shift, L, D) {
  if (any(k < 0)) {
    stop("all values in k must be >=0.", call. = FALSE)
  }
  if (shift < 0) {
    stop("shift must be >=0.", call. = FALSE)
  }
  if (L >= D) {
    stop("L must be smaller than D.", call. = FALSE)
  }
  invisible(NULL)
}

## Parameters of a named distribution with mean mu and standard deviation sigma
si_dist_args <- function(dist, mu, sigma) {
  switch(
    dist,
    gamma = list(shape = (mu / sigma)^2, scale = sigma^2 / mu),
    lognormal = {
      sdlog <- sqrt(log(sigma^2 / mu^2 + 1))
      list(meanlog = log(mu) - sdlog^2 / 2, sdlog = sdlog)
    },
    weibull = {
      cv2 <- (sigma / mu)^2
      shape <- exp(stats::uniroot(
        function(log_shape) {
          shape <- exp(log_shape)
          exp(lgamma(1 + 2 / shape) - 2 * lgamma(1 + 1 / shape)) - 1 - cv2
        },
        interval = c(-5, 7), tol = 1e-12
      )$root)
      list(shape = shape, scale = mu / gamma(1 + 1 / shape))
    }
  )
}
