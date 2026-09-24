#' Link primarycensored and EpiEstim
#'
#' [primary2estim()] transforms a serial interval distribution estimated with
#' [primarycensored][primarycensored::primarycensored-package] into a sample
#' of discrete serial interval distributions in the right format for input into
#' [estimate_R()] with method "si_from_sample". This allows the serial interval
#' to be estimated once and reused.
#'
#' Three inputs are supported:
#' - A maximum likelihood fit from [primarycensored::fitdistdoublecens()]. `n`
#'   parameter sets are drawn from the asymptotic normal distribution of the
#'   estimates, with positive parameters drawn on the log scale. This is the
#'   approach used by [estimate_R()] with method "si_from_data".
#' - A Bayesian fit from the model returned by
#'   [primarycensored::pcd_cmdstan_model()] (a `CmdStanMCMC` object). `n`
#'   evenly spaced posterior draws of the delay parameters are used, so this
#'   is the way to use a Bayesian fit of the serial interval. The draws are
#'   read with the fit's own `draws()` method, which needs the posterior
#'   package, installed alongside cmdstanr when the fit was made.
#' - A data frame of parameter draws, with one row per draw and one column per
#'   parameter of `dist`, named as the arguments of `dist` (e.g. `shape` and
#'   `scale` for [stats::pgamma()]).
#'
#' Fits from [fitdistrplus::fitdist()] and [fitdistrplus::fitdistcens()] are
#' also accepted. For most epidemiological data, a fit with
#' [primarycensored::fitdistdoublecens()] is preferred as it accounts for the
#' censoring of both the primary and secondary events, and for right
#' truncation. Only the parameter estimates and their covariance are used, and
#' their names are checked against the arguments of `dist`.
#'
#' Distributions beyond those supported by name in [estimate_R()] with method
#' "si_from_data" can be used by giving `log_params` for a maximum likelihood
#' fit, or `param_map` for a Stan fit.
#'
#' Each parameter set is discretised with [discr_si()]. The support runs up to
#' the largest 0.999 quantile of the primary censored serial interval across
#' the sample, and each distribution is right truncated at the end of this
#' support. When `shift` is smaller than 1, the probability of a serial
#' interval of zero is set to zero by truncating below 1.
#'
#' @param x A fit from [primarycensored::fitdistdoublecens()], a fit from the
#'   model returned by [primarycensored::pcd_cmdstan_model()], a fit from
#'   fitdistrplus, or a data frame of parameter draws (see details).
#' @param dist The cumulative distribution function of the serial interval
#'   minus `shift`, used for the fit, e.g. [stats::pgamma()]. Any cumulative
#'   distribution function of a non-negative delay can be used, with
#'   `log_params` or `param_map` for distributions not supported by name in
#'   [estimate_R()] with method "si_from_data".
#' @param n A positive integer giving the number of serial interval
#'   distributions to draw. For a Stan fit with fewer posterior draws, all
#'   draws are used.
#' @param shift A non-negative real giving the shift of the serial interval,
#'   i.e. the fit is to the serial interval minus `shift`. Use 1 for a fit to
#'   serial intervals shifted by one day.
#' @param si_discr_args A named list of additional arguments to [discr_si()],
#'   which can contain `L`, `D`, `dprimary` and `primary_args`. `dprimary` and
#'   `primary_args` should match those used for the fit.
#' @param seed An integer used as the seed for the random number generator
#'   when drawing from a [primarycensored::fitdistdoublecens()] fit.
#' @param log_params For a maximum likelihood fit; the names of the parameters
#'   that are positive, and so drawn on the log scale. Defaults to the positive
#'   parameters of the distributions supported by name in [estimate_R()] with
#'   method "si_from_data" (e.g. `shape` and `scale` for [stats::pgamma()]),
#'   and otherwise to the parameters with a positive estimate.
#' @param param_map For a Stan fit; a function taking a matrix of posterior
#'   draws of the parameters of the primarycensored Stan model (with columns
#'   `params[1]`, `params[2]`, ...) and returning a data frame of the
#'   corresponding parameters of `dist`, with one column per parameter. Only
#'   needed for distributions other than those supported by name in
#'   [estimate_R()] with method "si_from_data", e.g.
#'   `function(draws) data.frame(rate = draws[, 1])`.
#' @param ... Not used.
#'
#' @return A list with two elements:
#'   - `si_sample`: a matrix where each column gives one distribution of the
#'      serial interval to be explored, as used by [estimate_R()] with method
#'      "si_from_sample".
#'   - `si_parametric_distr`: the name of the parametric distribution of the
#'      serial interval, with "_offset_1" added when `shift` is 1.
#'
#' @seealso [estimate_R()], [discr_si()]
#'
#' @references
#' Abbott, S. et al. [primarycensored][primarycensored::primarycensored-package]:
#' Primary Event Censored Distributions. \doi{10.5281/zenodo.13632839}
#'
#' @importFrom fitdistrplus fitdist
#' @export
#' @examples
#' ## load data on rotavirus
#' data("MockRotavirus")
#' si_data <- MockRotavirus$si_data
#'
#' ## estimate the serial interval by maximum likelihood, accounting for
#' ## double interval censoring
#' censdata <- data.frame(
#'   left = si_data$SL - si_data$EL,
#'   right = si_data$SR - si_data$EL,
#'   pwindow = si_data$ER - si_data$EL,
#'   D = Inf
#' )
#' fit <- primarycensored::fitdistdoublecens(
#'   censdata, distr = "gamma", start = list(shape = 2, scale = 1)
#' )
#'
#' ## turn this into a sample of serial interval distributions
#' si_sample <- primary2estim(fit, dist = stats::pgamma, n = 100,
#'                            seed = 1)$si_sample
#'
#' ## use estimate_R to estimate the reproduction number
#' ## based on these estimates of the serial interval
#' R_si_from_sample <- estimate_R(
#'   MockRotavirus$incidence,
#'   method = "si_from_sample",
#'   si_sample = si_sample,
#'   config = make_config(list(n2 = 50))
#' )
#'
#' plot(R_si_from_sample)
primary2estim <- function(x, ...) {
  UseMethod("primary2estim")
}

#' @rdname primary2estim
#' @export
primary2estim.default <- function(x, ...) {
  stop("x must be a fit from primarycensored::fitdistdoublecens(), a fit ",
       "from primarycensored::pcd_cmdstan_model() or a data frame of ",
       "parameter draws.", call. = FALSE)
}

#' @rdname primary2estim
#' @export
primary2estim.fitdist <- function(x, dist, n = 1000, shift = 0,
                                  si_discr_args = list(), seed = NULL,
                                  log_params = NULL, ...) {
  if (missing(dist) || !is.function(dist)) {
    stop("dist must be given as a cumulative distribution function, e.g. ",
         "stats::pgamma.", call. = FALSE)
  }
  estimate <- x$estimate
  if (!is.numeric(estimate) || is.null(names(estimate))) {
    stop("x$estimate must be a named numeric vector of parameter estimates.",
         call. = FALSE)
  }
  unknown <- setdiff(names(estimate), names(formals(dist)))
  if (length(unknown) > 0) {
    stop("The parameters in x$estimate must be arguments of dist. Unknown ",
         "parameters: ", toString(unknown), ".", call. = FALSE)
  }
  fit_vcov <- x$vcov
  if (!is.null(fit_vcov) &&
        (!is.matrix(fit_vcov) ||
           !setequal(rownames(fit_vcov), names(estimate)))) {
    stop("x$vcov must be a matrix with rows and columns named as the ",
         "parameters in x$estimate.", call. = FALSE)
  }
  if (!is.null(fit_vcov)) {
    fit_vcov <- fit_vcov[names(estimate), names(estimate), drop = FALSE]
  }
  if (is.null(log_params)) {
    entry <- si_distribution_from_cdf(dist)
    log_params <- if (is.null(entry)) {
      names(estimate)[estimate > 0]
    } else {
      entry$params[entry$positive]
    }
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  samples <- draw_si_params(
    estimate, fit_vcov, n, names(estimate) %in% log_params
  )
  primary2estim(
    samples, dist = dist, shift = shift, si_discr_args = si_discr_args
  )
}

#' @rdname primary2estim
#' @export
primary2estim.fitdistcens <- primary2estim.fitdist

#' @rdname primary2estim
#' @export
primary2estim.CmdStanMCMC <- function(x, dist, n = 1000, shift = 0,
                                      si_discr_args = list(),
                                      param_map = NULL, ...) {
  if (missing(dist) || !is.function(dist)) {
    stop("dist must be given as a cumulative distribution function, e.g. ",
         "stats::pgamma.", call. = FALSE)
  }
  if (is.null(param_map)) {
    entry <- si_distribution_from_cdf(dist)
    if (is.null(entry)) {
      stop("There is no default mapping from the parameters of the Stan ",
           "model to those of dist. Give it with param_map.", call. = FALSE)
    }
    param_map <- entry$from_stan
  }
  if (!is.function(x$draws)) {
    stop("x must have a draws() method, as a fit from ",
         "primarycensored::pcd_cmdstan_model() does.", call. = FALSE)
  }
  draws <- x$draws(variables = "params", format = "draws_matrix")
  draws <- matrix(
    as.numeric(draws), nrow = nrow(draws),
    dimnames = list(NULL, colnames(draws))
  )
  n_draws <- nrow(draws)
  keep <- unique(round(seq(1, n_draws, length.out = min(n, n_draws))))
  samples <- as.data.frame(param_map(draws[keep, , drop = FALSE]))
  primary2estim(
    samples, dist = dist, shift = shift, si_discr_args = si_discr_args
  )
}

#' @rdname primary2estim
#' @export
primary2estim.data.frame <- function(x, dist, shift = 0,
                                     si_discr_args = list(), ...) {
  if (missing(dist) || !is.function(dist)) {
    stop("dist must be given as a cumulative distribution function, e.g. ",
         "stats::pgamma.", call. = FALSE)
  }
  unknown <- setdiff(names(x), names(formals(dist)))
  if (length(unknown) > 0) {
    stop("The columns of x must be arguments of dist. Unknown columns: ",
         toString(unknown), ".", call. = FALSE)
  }
  discr_args <- primary2estim_discr_args(si_discr_args)
  si_args <- si_sample_support(x, dist, shift, discr_args)
  params <- lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))
  si_sample <- t(discr_si_param_draws(si_args$k, params, si_args))

  list(
    si_sample = si_sample,
    si_parametric_distr = si_cdf_name(dist, shift)
  )
}

## Check the extra arguments to discr_si given to primary2estim and fill in
## the defaults of the primary event distribution
primary2estim_discr_args <- function(si_discr_args) {
  if (is.null(si_discr_args)) {
    si_discr_args <- list()
  }
  if (!is.null(si_discr_args$dist) || !is.null(si_discr_args$shift)) {
    stop("si_discr_args cannot contain dist or shift. Use the dist and ",
         "shift arguments instead.", call. = FALSE)
  }
  check_si_discr_args(si_discr_args[setdiff(names(si_discr_args), "L")])
  if (is.null(si_discr_args$dprimary)) {
    si_discr_args$dprimary <- stats::dunif
  }
  if (is.null(si_discr_args$primary_args)) {
    si_discr_args$primary_args <- list()
  }
  si_discr_args
}

## Support and truncation of a sample of serial interval distributions.
## Following the approach of the distspec package, the support runs to the
## largest 0.999 quantile of the primary censored serial interval across the
## sample and each distribution is right truncated at the end of this support.
## Returns the complete list of extra arguments of discr_si, with the support
## in k.
si_sample_support <- function(samples, dist, shift, discr_args) {
  q_max <- max(vnapply(seq_len(nrow(samples)), function(i) {
    do.call(
      primarycensored::qprimarycensored,
      c(
        list(
          p = 0.999, pdist = dist, pwindow = 1,
          dprimary = discr_args$dprimary,
          primary_args = discr_args$primary_args
        ),
        as.list(samples[i, , drop = FALSE]),
        list(check = FALSE)
      )
    )
  }))
  max_value <- ceiling(q_max + shift)
  lower <- if (shift < 1) 1 else -Inf
  if (!is.null(discr_args$L)) {
    lower <- max(lower, discr_args$L)
  }
  upper <- max_value + 1
  if (!is.null(discr_args$D)) {
    upper <- min(upper, discr_args$D)
  }
  si_args <- si_discr_defaults(list(
    dist = dist, shift = shift, L = lower, D = upper,
    dprimary = discr_args$dprimary, primary_args = discr_args$primary_args
  ))
  si_args$k <- seq(0, max_value)
  si_args
}

## Name of a serial interval distribution given its CDF and shift
si_cdf_name <- function(dist, shift) {
  pdists <- primarycensored::pcd_distributions
  pdists <- pdists[!is.na(pdists$pdist), ]
  found <- vapply(pdists$pdist, function(p) {
    fn <- tryCatch(
      getExportedValue("stats", p),
      error = function(e) {
        tryCatch(
          getExportedValue("primarycensored", p),
          error = function(e) NULL
        )
      }
    )
    identical(fn, dist)
  }, logical(1))
  name <- if (any(found)) pdists$aliases[found][1] else "custom"
  if (shift == 1) {
    name <- paste0(name, "_offset_1")
  }
  name
}
