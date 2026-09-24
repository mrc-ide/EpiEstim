MockRotavirus <- mock_rotavirus()

# A minimal stand in for a CmdStanMCMC fit from
# primarycensored::pcd_cmdstan_model(), returning the given draws of params
mock_cmdstan_fit <- function(params) {
  colnames(params) <- paste0("params[", seq_len(ncol(params)), "]")
  structure(
    list(draws = function(variables = NULL, format = NULL, ...) params),
    class = c("CmdStanMCMC", "CmdStanFit", "R6")
  )
}

# Expected sample of serial interval distributions for parameter draws, on
# the support used by primary2estim
expected_si_sample <- function(samples, dist, k, shift = 0, L = 1, ...) {
  vapply(seq_len(nrow(samples)), function(i) {
    do.call(discr_si, c(
      list(k = k, dist = dist, shift = shift, L = L, D = length(k), ...),
      as.list(samples[i, , drop = FALSE])
    ))
  }, numeric(length(k)))
}

test_that("primary2estim of a fitdistdoublecens fit matches si_from_data", {
  si_data <- process_si_data(MockRotavirus$si_data)
  config <- quiet_config(si_parametric_distr = "gamma")
  censdata <- si_data_to_censdata(si_data, 0)
  start <- init_mcmc_params(si_data, "gamma")
  fit <- primarycensored::fitdistdoublecens(
    censdata, distr = "gamma",
    start = list(shape = start[1], scale = start[2])
  )
  out <- primary2estim(fit, dist = stats::pgamma, n = 100, seed = 1)
  res <- suppressWarnings(run_si_from_data(MockRotavirus$si_data, config))
  expect_named(out, c("si_sample", "si_parametric_distr"))
  expect_identical(out$si_parametric_distr, "gamma")
  expect_identical(ncol(out$si_sample), 100L)
  support <- seq_len(nrow(out$si_sample))
  expect_equal(
    unname(t(out$si_sample)), unname(res$si_distr[, support]),
    tolerance = 1e-10
  )
  expect_true(all(res$si_distr[, -support] == 0))
})

test_that("primary2estim of parameter draws matches discr_si", {
  samples <- data.frame(shape = c(2, 3), scale = c(1, 1.5))
  out <- primary2estim(samples, dist = stats::pgamma)
  k <- seq(0, nrow(out$si_sample) - 1)
  expect_equal(
    out$si_sample, expected_si_sample(samples, stats::pgamma, k),
    tolerance = 1e-10
  )
  expect_true(all(out$si_sample[1, ] == 0))
  expect_equal(colSums(out$si_sample), c(1, 1), tolerance = 1e-10)
  expect_identical(out$si_parametric_distr, "gamma")
})

test_that("primary2estim supports any CDF and primary distribution", {
  samples <- data.frame(shape = c(1.5, 2.5), scale = c(3, 5))
  primary <- list(
    dprimary = primarycensored::dexpgrowth, primary_args = list(r = 0.3)
  )
  out <- primary2estim(
    samples, dist = stats::pweibull, si_discr_args = primary
  )
  k <- seq(0, nrow(out$si_sample) - 1)
  expected <- expected_si_sample(
    samples, stats::pweibull, k,
    dprimary = primarycensored::dexpgrowth, primary_args = list(r = 0.3)
  )
  expect_equal(out$si_sample, expected, tolerance = 1e-10)
  expect_identical(out$si_parametric_distr, "weibull")
})

test_that("primary2estim supports offset distributions", {
  samples <- data.frame(meanlog = c(1, 1.2), sdlog = c(0.5, 0.4))
  out <- primary2estim(samples, dist = stats::plnorm, shift = 1)
  k <- seq(0, nrow(out$si_sample) - 1)
  expect_equal(
    out$si_sample,
    expected_si_sample(samples, stats::plnorm, k, shift = 1, L = -Inf),
    tolerance = 1e-10
  )
  expect_true(all(out$si_sample[1, ] == 0))
  expect_true(all(out$si_sample[2, ] > 0))
  expect_identical(out$si_parametric_distr, "lognormal_offset_1")
})

test_that("primary2estim uses L and D from si_discr_args", {
  samples <- data.frame(shape = c(2, 3), scale = c(1, 1.5))
  out <- primary2estim(
    samples, dist = stats::pgamma, si_discr_args = list(L = 2, D = 5)
  )
  k <- seq(0, nrow(out$si_sample) - 1)
  expect_true(all(out$si_sample[k < 2 | k >= 5, ] == 0))
  expect_equal(colSums(out$si_sample), c(1, 1), tolerance = 1e-10)
})

test_that("primary2estim of a Stan fit uses its posterior draws", {
  shape <- c(2, 2.5, 3, 3.5, 4, 4.5)
  rate <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5)
  fit <- mock_cmdstan_fit(cbind(shape, rate))
  out <- primary2estim(fit, dist = stats::pgamma, n = 3)
  samples <- data.frame(
    shape = shape[c(1, 4, 6)], scale = 1 / rate[c(1, 4, 6)]
  )
  k <- seq(0, nrow(out$si_sample) - 1)
  expect_identical(ncol(out$si_sample), 3L)
  expect_equal(
    out$si_sample, expected_si_sample(samples, stats::pgamma, k),
    tolerance = 1e-10
  )
  all_draws <- primary2estim(fit, dist = stats::pgamma, n = 100)
  expect_identical(ncol(all_draws$si_sample), 6L)
})

test_that("primary2estim maps Stan parameters for lognormal and Weibull", {
  draws <- cbind(c(1, 1.2), c(0.5, 0.4))
  lnorm_out <- primary2estim(mock_cmdstan_fit(draws), dist = stats::plnorm)
  expect_equal(
    lnorm_out$si_sample,
    primary2estim(
      data.frame(meanlog = draws[, 1], sdlog = draws[, 2]),
      dist = stats::plnorm
    )$si_sample
  )
  weibull_out <- primary2estim(mock_cmdstan_fit(draws), dist = stats::pweibull)
  expect_equal(
    weibull_out$si_sample,
    primary2estim(
      data.frame(shape = draws[, 1], scale = draws[, 2]),
      dist = stats::pweibull
    )$si_sample
  )
})

test_that("primary2estim errors without dist", {
  samples <- data.frame(shape = c(2, 3), scale = c(1, 1.5))
  fit <- mock_cmdstan_fit(cbind(c(2, 3), c(1, 0.5)))
  expect_error(primary2estim(samples), "dist")
  expect_error(primary2estim(fit), "dist")
})

test_that("primary2estim errors for unsupported inputs", {
  fit <- mock_cmdstan_fit(cbind(c(2, 3)))
  expect_error(primary2estim(fit, dist = stats::pexp), "pgamma")
  samples <- data.frame(shape = c(2, 3), rate = c(1, 0.5))
  expect_error(primary2estim(samples, dist = stats::pweibull), "rate")
  expect_error(
    primary2estim(
      data.frame(shape = 2, scale = 1), dist = stats::pgamma,
      si_discr_args = list(dist = stats::plnorm)
    ),
    "dist"
  )
  expect_error(primary2estim(list(1)), "fitdistdoublecens")
})

test_that("primary2estim works with a real pcd_cmdstan_model fit", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)))
  set.seed(1)
  delays <- primarycensored::rprimarycensored(
    200, rdist = function(n) stats::rgamma(n, shape = 3, rate = 0.8),
    pwindow = 1, swindow = 1, D = Inf
  )
  delay_counts <- as.data.frame(table(delay = delays))
  delay_counts$delay <- as.numeric(as.character(delay_counts$delay))
  names(delay_counts)[2] <- "n"
  delay_counts$delay_upper <- delay_counts$delay + 1
  delay_counts$pwindow <- 1
  delay_counts$relative_obs_time <- Inf
  stan_data <- primarycensored::pcd_as_stan_data(
    delay_counts,
    dist_id = primarycensored::pcd_stan_dist_id("gamma", "delay"),
    primary_id = primarycensored::pcd_stan_dist_id("uniform", "primary"),
    param_bounds = list(lower = c(0, 0), upper = c(Inf, Inf)),
    primary_param_bounds = list(lower = numeric(0), upper = numeric(0)),
    priors = list(location = c(2, 1), scale = c(1, 1)),
    primary_priors = list(location = numeric(0), scale = numeric(0))
  )
  model <- suppressMessages(
    primarycensored::pcd_cmdstan_model(dir = tempdir())
  )
  fit <- suppressMessages(model$sample(
    data = stan_data, chains = 1, iter_warmup = 300, iter_sampling = 200,
    refresh = 0, show_messages = FALSE, seed = 1
  ))
  out <- primary2estim(fit, dist = stats::pgamma, n = 50)
  expect_identical(ncol(out$si_sample), 50L)
  expect_equal(colSums(out$si_sample), rep(1, 50), tolerance = 1e-8)
  k <- seq(0, nrow(out$si_sample) - 1)
  expect_equal(mean(colSums(k * out$si_sample)), 3 / 0.8, tolerance = 0.1)
})
