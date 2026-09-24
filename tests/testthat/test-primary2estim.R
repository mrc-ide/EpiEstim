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
  fit <- primarycensored::fitdistdoublecens(
    censdata, distr = "gamma", start = si_start_values(si_data, "gamma")
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
  pcustom <- function(q, rate) stats::pexp(q, rate)
  fit <- mock_cmdstan_fit(cbind(c(2, 3)))
  expect_error(primary2estim(fit, dist = pcustom), "param_map")
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

test_that("primary2estim errors for a Stan fit without a draws method", {
  fit <- structure(list(), class = "CmdStanMCMC")
  expect_error(primary2estim(fit, dist = stats::pgamma), "draws")
})

test_that("primary2estim supports the exponential distribution", {
  rate <- c(0.2, 0.3, 0.25)
  fit <- mock_cmdstan_fit(cbind(rate))
  out <- primary2estim(fit, dist = stats::pexp)
  k <- seq(0, nrow(out$si_sample) - 1)
  expect_equal(
    out$si_sample,
    expected_si_sample(data.frame(rate = rate), stats::pexp, k),
    tolerance = 1e-10
  )
  expect_identical(out$si_parametric_distr, "exponential")
})

test_that("primary2estim can be extended with param_map and log_params", {
  pcustom <- function(q, rate) stats::pexp(q, rate)
  rate <- c(0.2, 0.3)
  fit <- mock_cmdstan_fit(cbind(rate))
  out <- primary2estim(
    fit, dist = pcustom,
    param_map = function(draws) data.frame(rate = draws[, 1])
  )
  expect_equal(
    out$si_sample,
    primary2estim(data.frame(rate = rate), dist = stats::pexp)$si_sample
  )
  ml_fit <- structure(
    list(estimate = c(rate = 0.25), vcov = matrix(1e-4, 1, 1,
         dimnames = list("rate", "rate"))),
    class = "fitdist"
  )
  custom <- primary2estim(
    ml_fit, dist = pcustom, log_params = "rate", n = 20, seed = 3
  )
  exp_out <- primary2estim(ml_fit, dist = stats::pexp, n = 20, seed = 3)
  expect_equal(custom$si_sample, exp_out$si_sample)
})

test_that("primary2estim accepts fitdistrplus fits", {
  set.seed(2)
  si <- stats::rgamma(200, shape = 4, scale = 1.2)
  fit <- fitdistrplus::fitdist(si, "gamma", start = list(shape = 2, scale = 1))
  out <- primary2estim(fit, dist = stats::pgamma, n = 50, seed = 1)
  expect_identical(ncol(out$si_sample), 50L)
  cens <- data.frame(left = floor(si), right = floor(si) + 1)
  cens_fit <- fitdistrplus::fitdistcens(
    cens, "lnorm", start = list(meanlog = 1, sdlog = 0.5)
  )
  cens_out <- primary2estim(cens_fit, dist = stats::plnorm, n = 50, seed = 1)
  expect_identical(ncol(cens_out$si_sample), 50L)
  expect_equal(colSums(cens_out$si_sample), rep(1, 50), tolerance = 1e-8)
})

test_that("primary2estim checks the contents of fits", {
  bad <- structure(list(estimate = c(a = 1, b = 2)), class = "fitdist")
  expect_error(primary2estim(bad, dist = stats::pgamma), "estimate")
  no_names <- structure(list(estimate = c(1, 2)), class = "fitdistcens")
  expect_error(primary2estim(no_names, dist = stats::pgamma), "named")
  bad_vcov <- structure(
    list(estimate = c(shape = 2, scale = 1), vcov = matrix(1, 2, 2)),
    class = "fitdist"
  )
  expect_error(primary2estim(bad_vcov, dist = stats::pgamma), "vcov")
})

test_that("draws from a maximum likelihood fit match the delta method", {
  estimate <- c(shape = 4, scale = 1.5)
  vcov <- matrix(c(0.5, -0.12, -0.12, 0.04), 2, 2)
  set.seed(10)
  draws <- draw_si_params(estimate, vcov, 2e5, c(TRUE, TRUE))
  log_draws <- log(as.matrix(draws))
  jacobian <- diag(1 / estimate)
  expect_equal(unname(colMeans(log_draws)), unname(log(estimate)),
               tolerance = 0.01)
  expect_equal(unname(stats::cov(log_draws)),
               unname(jacobian %*% vcov %*% jacobian), tolerance = 0.02)
})

test_that("draws from a maximum likelihood fit cover the true serial interval", {
  skip_on_cran()
  mu <- 5
  sigma <- 2.5
  covered <- vapply(seq_len(100), function(i) {
    si_data <- simulate_si_data(50, mu, sigma, obs_time = 200, seed = i)
    censdata <- si_data_to_censdata(si_data, 0)
    fit <- primarycensored::fitdistdoublecens(
      censdata, distr = "gamma", start = list(shape = 2, scale = 2)
    )
    set.seed(i)
    draws <- draw_si_params(
      fit$estimate[c("shape", "scale")],
      fit$vcov[c("shape", "scale"), c("shape", "scale")],
      1000, c(TRUE, TRUE)
    )
    interval <- stats::quantile(draws$shape * draws$scale, c(0.025, 0.975))
    interval[1] <= mu && mu <= interval[2]
  }, logical(1))
  expect_gt(mean(covered), 0.88)
  expect_lt(mean(covered), 0.99)
})
