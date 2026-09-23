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
#' default), ends on day \eqn{k}. It is computed using the primarycensored
#' package (Abbott et al.), which also allows other delay distributions
#' (`dist`), a different
#' `shift`, other distributions of the primary event time within the day
#' (`dprimary`), and truncation of the serial interval below `L` and at or
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
#' Abbott, S. et al. primarycensored: Primary Event Censored Distributions.
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
           "pass the parameters of dist through ... instead.")
    }
    pdist <- dist
    dist_args <- list(...)
  } else {
    if (sigma < 0) {
      stop("sigma must be >=0.")
    }
    if (mu <= shift) {
      stop("mu must be >", shift)
    }
    pdist <- si_pdist(dist)
    dist_args <- si_dist_args(dist, mu - shift, sigma)
  }
  if (any(k < 0)) {
    stop("all values in k must be >=0.")
  }
  if (shift < 0) {
    stop("shift must be >=0.")
  }
  if (L >= D) {
    stop("L must be smaller than D.")
  }

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

## Check the extra arguments to discr_si given in config$si_discr_args
check_si_discr_args <- function(si_discr_args) {
  allowed <- c("dist", "shift", "L", "D", "dprimary", "primary_args")
  if (!is.list(si_discr_args) ||
        (length(si_discr_args) > 0 && is.null(names(si_discr_args)))) {
    stop("si_discr_args must be a named list.")
  }
  unknown <- setdiff(names(si_discr_args), allowed)
  if (length(unknown) > 0) {
    stop("si_discr_args contains unsupported arguments: ",
         toString(unknown), ". Supported arguments are: ",
         toString(allowed), ".")
  }
  if (!is.null(si_discr_args$dist) && !is.character(si_discr_args$dist)) {
    stop("si_discr_args$dist must be one of 'gamma', 'lognormal' or ",
         "'weibull'.")
  }
  support <- utils::modifyList(list(shift = 1, L = -Inf), si_discr_args)
  if (support$shift < 1 && support$L < 1) {
    stop("si_discr_args gives a non-zero probability of a serial interval ",
         "of zero, which EpiEstim does not allow. Use shift >= 1 or L >= 1.")
  }
  invisible(NULL)
}

## Cumulative distribution function for a named serial interval distribution
si_pdist <- function(dist) {
  switch(
    dist,
    gamma = stats::pgamma,
    lognormal = stats::plnorm,
    weibull = stats::pweibull,
    stop("dist must be one of 'gamma', 'lognormal' or 'weibull', ",
         "or a cumulative distribution function.")
  )
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
