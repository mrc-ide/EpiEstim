# Frozen copy of the discr_si() implementation from EpiEstim <= 3.0.0.
# Used to check that the primarycensored based implementation reproduces
# the original closed form solution for the shifted Gamma.
# nolint start: condition_call_linter. Kept as in EpiEstim <= 3.0.0.
old_discr_si <- function(k, mu, sigma) {
  if (sigma < 0) {
    stop("sigma must be >=0.")
  }
  if (mu <= 1) {
    stop("mu must be >1")
  }
  if (any(k < 0)) {
    stop("all values in k must be >=0.")
  }

  a <- ((mu - 1) / sigma)^2
  b <- sigma^2 / (mu - 1)

  cdf_gamma <- function(k, a, b) stats::pgamma(k, shape = a, scale = b)

  res <- k * cdf_gamma(k, a, b) +
    (k - 2) * cdf_gamma(k - 2, a, b) - 2 * (k - 1) * cdf_gamma(k - 1, a, b)
  res <- res + a * b * (2 * cdf_gamma(k - 1, a + 1, b) -
                          cdf_gamma(k - 2, a + 1, b) - cdf_gamma(k, a + 1, b))
  res <- vapply(res, function(e) max(0, e), numeric(1))

  return(res)
}
# nolint end

# Mean of a discrete distribution on k
pmf_mean <- function(k, w) sum(k * w)

# Original discr_si applied to several draws, with the interface of the
# internal discr_si_draws
old_discr_si_draws <- function(k, mu, sigma, si_discr_args = NULL) {
  t(vapply(
    seq_along(mu), function(i) old_discr_si(k, mu[i], sigma[i]),
    numeric(length(k))
  ))
}
