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
#' @param dist The cumulative distribution function of the serial interval
#' minus `shift`. Defaults to [stats::pgamma()]. With [stats::pgamma()] or
#' [stats::plnorm()] the distribution can be given by `mu` and `sigma`.
#' Otherwise, or instead, pass the parameters of `dist` through `...`, e.g.
#' `dist = stats::pweibull, shape = 2, scale = 3`.
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
#' @param ... Parameters of `dist`, used when `mu` and `sigma` are not given.
#'
#' @return Gives the discrete probability \eqn{w_k} that the serial interval is
#' equal to \eqn{k}. This is not normalised over `k`, so it sums to less than 1
#' if `k` does not cover the support of the distribution. Use `L` and `D` to
#' truncate and normalise the distribution.
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
#' discr_si(seq(0, 20), mean_flu_si, sd_flu_si, dist = stats::plnorm,
#'          shift = 0, L = 1)
#'
#' ## Using a Weibull serial interval given by its shape and scale
#' discr_si(seq(0, 20), dist = stats::pweibull, shape = 1.5, scale = 2)

discr_si <- function(k, mu, sigma, dist = stats::pgamma, shift = 1, L = -Inf,
                     D = Inf, dprimary = stats::dunif, primary_args = list(),
                     ...)
{
  if (!is.function(dist)) {
    stop("dist must be a cumulative distribution function, e.g. ",
         "stats::pgamma.", call. = FALSE)
  }
  dist_args <- list(...)
  if (!missing(mu) || !missing(sigma)) {
    if (missing(mu) || missing(sigma)) {
      stop("Both mu and sigma must be given together.", call. = FALSE)
    }
    if (length(dist_args) > 0) {
      stop("Specify either mu and sigma or the parameters of dist, ",
           "not both.", call. = FALSE)
    }
    check_si_moments(mu, sigma, shift)
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
          x = k[in_support] - shift, pdist = dist, pwindow = 1,
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
## deviation. Gives the same result as calling discr_si for each draw.
## Returns a matrix with one row per draw.
discr_si_draws <- function(k, mu, sigma, si_discr_args = NULL) {
  if (is.null(si_discr_args)) {
    si_discr_args <- list()
  }
  check_si_discr_args(si_discr_args)
  si_args <- si_discr_defaults(si_discr_args)
  check_si_moments(mu, sigma, si_args$shift)
  params <- lapply(seq_along(mu), function(i) {
    si_dist_args(si_args$dist, mu[i] - si_args$shift, sigma[i])
  })
  discr_si_param_draws(k, params, si_args)
}

## Defaults of the extra arguments of discr_si
si_discr_defaults <- function(si_discr_args = list()) {
  defaults <- lapply(
    formals(discr_si)[c("dist", "shift", "L", "D", "dprimary", "primary_args")],
    eval
  )
  utils::modifyList(defaults, si_discr_args)
}

## Discretise the serial interval for several draws of the parameters of
## si_args$dist. params is a list with one named list of parameters per draw
## and si_args a complete list of the extra arguments of discr_si (see
## si_discr_defaults). The primarycensored object is built once and its
## parameters updated for each draw. Returns a matrix with one row per draw.
discr_si_param_draws <- function(k, params, si_args) {
  check_si_support(k, si_args$shift, si_args$L, si_args$D)
  res <- matrix(0, nrow = length(params), ncol = length(k))
  in_support <- k >= si_args$L & k < si_args$D
  if (!any(in_support) || length(params) == 0) {
    return(res)
  }

  pcens <- do.call(
    primarycensored::new_pcens,
    c(
      list(pdist = si_args$dist, dprimary = si_args$dprimary,
           primary_args = si_args$primary_args),
      params[[1]]
    )
  )
  x <- k[in_support] - si_args$shift
  for (i in seq_along(params)) {
    pcens <- do.call(stats::update, c(list(pcens), params[[i]]))
    pmf <- primarycensored::pcens_pmf(
      pcens, x, pwindow = 1,
      L = si_args$L - si_args$shift, D = si_args$D - si_args$shift
    )
    res[i, in_support] <- pmax(pmf, 0)
  }

  return(res)
}

## Check the extra arguments to discr_si given in config$si_discr_args
check_si_discr_args <- function(si_discr_args) {
  allowed <- c("dist", "shift", "L", "D", "dprimary", "primary_args")
  if (!is.list(si_discr_args) ||
        (length(si_discr_args) > 0 &&
           (is.null(names(si_discr_args)) || !all(nzchar(names(si_discr_args)))))) {
    stop("si_discr_args must be a named list.", call. = FALSE)
  }
  unknown <- setdiff(names(si_discr_args), allowed)
  if (length(unknown) > 0) {
    stop("si_discr_args contains unsupported arguments: ",
         toString(unknown), ". Supported arguments are: ",
         toString(allowed), ".", call. = FALSE)
  }
  if (!is.null(si_discr_args$dist) &&
        is.na(si_moment_dist(si_discr_args$dist))) {
    stop("si_discr_args$dist must be stats::pgamma or stats::plnorm, as the ",
         "serial interval is given by its mean and standard deviation.",
         call. = FALSE)
  }
  support <- utils::modifyList(list(shift = 1, L = -Inf), si_discr_args)
  if (support$shift < 1 && support$L < 1) {
    stop("si_discr_args gives a non-zero probability of a serial interval ",
         "of zero, which EpiEstim does not allow. Use shift >= 1 or L >= 1.",
         call. = FALSE)
  }
  invisible(NULL)
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

## Name of a distribution that can be given by its mean and standard
## deviation, or NA
si_moment_dist <- function(dist) {
  if (identical(dist, stats::pgamma)) {
    "gamma"
  } else if (identical(dist, stats::plnorm)) {
    "lognormal"
  } else {
    NA_character_
  }
}

## Parameters of dist with mean mu and standard deviation sigma
si_dist_args <- function(dist, mu, sigma) {
  moment_dist <- si_moment_dist(dist)
  if (is.na(moment_dist)) {
    stop("mu and sigma can only be used with dist = stats::pgamma or ",
         "stats::plnorm. For other distributions pass the parameters of ",
         "dist through ... instead.", call. = FALSE)
  }
  if (moment_dist == "gamma") {
    return(list(shape = (mu / sigma)^2, scale = sigma^2 / mu))
  }
  sdlog <- sqrt(log(sigma^2 / mu^2 + 1))
  list(meanlog = log(mu) - sdlog^2 / 2, sdlog = sdlog)
}
