## Serial interval distributions supported by name in estimate_R() with
## method "si_from_data", and without extra arguments in primary2estim().
## Entries are named as in primarycensored::pcd_distributions, which gives
## their CDF (pdist) and descriptive name (aliases). Each entry gives:
## - params: the parameters of the CDF
## - positive: which parameters are positive, and so drawn on the log scale
## - from_stan: the parameters of the CDF from draws of the parameters of the
##   primarycensored Stan model (a matrix with columns params[1], ...)
## - start: starting values for the fit from the mean and standard deviation
##   of the serial interval
si_distributions <- list(
  gamma = list(
    params = c("shape", "scale"),
    positive = c(TRUE, TRUE),
    ## Stan's gamma is parameterised by shape and rate
    from_stan = function(draws) {
      data.frame(shape = draws[, 1], scale = 1 / draws[, 2])
    },
    start = function(mu, sigma) {
      list(shape = (mu / sigma)^2, scale = sigma^2 / mu)
    }
  ),
  lnorm = list(
    params = c("meanlog", "sdlog"),
    positive = c(FALSE, TRUE),
    from_stan = function(draws) {
      data.frame(meanlog = draws[, 1], sdlog = draws[, 2])
    },
    start = function(mu, sigma) {
      sdlog <- sqrt(log(sigma^2 / mu^2 + 1))
      list(meanlog = log(mu) - sdlog^2 / 2, sdlog = sdlog)
    }
  ),
  weibull = list(
    params = c("shape", "scale"),
    positive = c(TRUE, TRUE),
    from_stan = function(draws) {
      data.frame(shape = draws[, 1], scale = draws[, 2])
    },
    ## shape from the approximation shape = cv^-1.086 (Justus et al. 1978)
    start = function(mu, sigma) {
      shape <- (sigma / mu)^-1.086
      list(shape = shape, scale = mu / gamma(1 + 1 / shape))
    }
  ),
  exp = list(
    params = "rate",
    positive = TRUE,
    from_stan = function(draws) data.frame(rate = draws[, 1]),
    start = function(mu, sigma) list(rate = 1 / mu)
  )
)

## Registry entry of primarycensored::pcd_distributions for each supported
## distribution, with its CDF and descriptive name
si_distribution_registry <- function() {
  pdists <- primarycensored::pcd_distributions
  pdists <- pdists[pdists$name %in% names(si_distributions), ]
  lapply(stats::setNames(seq_len(nrow(pdists)), pdists$name), function(i) {
    list(
      name = pdists$name[i],
      alias = pdists$aliases[i],
      pdist = getExportedValue("stats", pdists$pdist[i])
    )
  })
}

## Details of a supported distribution given its CDF, or NULL
si_distribution_from_cdf <- function(dist) {
  for (entry in si_distribution_registry()) {
    if (identical(entry$pdist, dist)) {
      return(c(entry, si_distributions[[entry$name]]))
    }
  }
  NULL
}

## Details of a supported distribution given its descriptive name
## (e.g. "gamma" or "lognormal"), or NULL
si_distribution_from_alias <- function(alias) {
  for (entry in si_distribution_registry()) {
    if (identical(entry$alias, alias)) {
      return(c(entry, si_distributions[[entry$name]]))
    }
  }
  NULL
}

## Names of the serial interval distributions supported by name
si_distribution_aliases <- function() {
  unname(vapply(si_distribution_registry(), function(entry) entry$alias, ""))
}

#' Starting values for the estimation of the serial interval
#'
#' Computes starting values for the parameters of the serial interval
#' distribution, used to initialise the maximum likelihood estimation of the
#' serial interval in [estimate_R()] with method "si_from_data". The
#' starting values match the mean and standard deviation of naive serial
#' intervals computed from the midpoints of the intervals in `si_data`.
#'
#' @param si_data The data on dates of symptoms of pairs of infector/infected
#'   individuals, as described in [estimate_R()].
#' @param dist The parametric distribution of the serial interval, as in
#'   `si_parametric_distr` in [make_config()], e.g. "gamma" or
#'   "gamma_offset_1".
#'
#' @return A named list of starting values for the parameters of the
#'   distribution, e.g. `shape` and `scale` for the Gamma distribution.
#'
#' @seealso [estimate_R()], [make_config()]
#' @export
#' @examples
#' data("MockRotavirus")
#' si_start_values(MockRotavirus$si_data, "gamma")
si_start_values <- function(si_data, dist) {
  fit_distr <- si_fit_distr(dist)
  naive_si <- (si_data$SR + si_data$SL) / 2 - (si_data$ER + si_data$EL) / 2
  ## avoid issues when the mean serial interval is below the shift
  si_mean <- max(mean(naive_si) - fit_distr$shift, 0.1)
  start_values <- fit_distr$start(si_mean, stats::sd(naive_si))
  if (anyNA(unlist(start_values))) {
    stop("NA result. Check that si_data is in the right format.",
         call. = FALSE)
  }
  start_values
}
