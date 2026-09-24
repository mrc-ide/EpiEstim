MockRotavirus <- mock_rotavirus()

test_that("si_from_data keeps the output structure", {
  res <- suppressWarnings(run_si_from_data(
    MockRotavirus$si_data, quiet_config(si_parametric_distr = "gamma")
  ))
  expect_s3_class(res, "estimate_R")
  expect_named(res, c(
    "R", "method", "si_distr", "SI.Moments", "dates", "I", "I_local",
    "I_imported", "si_fit_converged", "MCMC_converged"
  ))
  expect_identical(res$method, "si_from_data")
  expect_identical(nrow(res$si_distr), 100L)
  expect_true(all(res$si_distr[, 1] == 0))
  expect_equal(rowSums(res$si_distr), rep(1, 100), tolerance = 1e-8)
  expect_named(res$SI.Moments, c("Mean", "Std"))
  expect_identical(nrow(res$SI.Moments), 100L)
  expect_true(res$si_fit_converged)
  expect_identical(res$MCMC_converged, res$si_fit_converged)
})

test_that("si_from_data works without mcmc_control and uses the config seed", {
  config <- suppressMessages(make_config(
    incid = MockRotavirus$incidence,
    list(si_parametric_distr = "gamma", n1 = 50, n2 = 10, seed = 4)
  ))
  expect_null(config$mcmc_control)
  run <- function() {
    collect <- collect_warnings(
      res <- run_si_from_data(MockRotavirus$si_data, config)
    )
    list(res = res, warnings = collect)
  }
  first <- run()
  second <- run()
  expect_identical(first$res$si_distr, second$res$si_distr)
  expect_false(any(grepl("burnin", first$warnings, fixed = TRUE)))
})

test_that("si_from_data supports the exponential distribution", {
  si_data <- simulate_si_data(300, 5, 5, obs_time = 200)
  config <- suppressMessages(
    si_from_data_config(
      MockRotavirus$incidence, si_parametric_distr = "exponential"
    )
  )
  res <- suppressWarnings(run_si_from_data(si_data, config))
  expect_equal(mean(res$SI.Moments$Mean), 5, tolerance = 0.6, scale = 1)
})

test_that("si_start_values gives moment based starting values", {
  si_data <- simulate_si_data(500, 6, 3, obs_time = 300)
  for (dist in c("gamma", "lognormal", "weibull", "exponential")) {
    start <- si_start_values(si_data, dist)
    expect_type(start, "list")
    expect_true(all(unlist(start[names(start) != "meanlog"]) > 0))
  }
  weibull <- si_start_values(si_data, "weibull")
  mean_weibull <- weibull$scale * gamma(1 + 1 / weibull$shape)
  naive <- (si_data$SR + si_data$SL) / 2 - (si_data$ER + si_data$EL) / 2
  expect_equal(mean_weibull, mean(naive), tolerance = 0.01)
  expect_named(si_start_values(si_data, "lognormal_offset_1"),
               c("meanlog", "sdlog"))
})

test_that("init_mcmc_params is deprecated in favour of si_start_values", {
  si_data <- MockRotavirus$si_data
  expect_warning(old <- init_mcmc_params(si_data, "gamma"), "si_start_values")
  expect_equal(old, unname(unlist(si_start_values(si_data, "gamma"))))
})

test_that("si_from_data is close to coarseDataTools MCMC results", {
  # Reference medians across the sample of serial interval distributions from
  # EpiEstim 3.0.0 using coarseDataTools::dic.fit.mcmc with
  # make_mcmc_control(burnin = 1000, thin = 10, seed = 1), n1 = 500, n2 = 50
  # and seed = 2 on MockRotavirus. Medians are compared because the MCMC
  # posterior for the lognormal has a heavier upper tail than the normal
  # approximation used now.
  reference <- list(
    gamma = c(mean = 2.020, sd = 1.377),
    weibull = c(mean = 2.065, sd = 1.435),
    lognormal = c(mean = 2.332, sd = 2.628)
  )
  for (dist in names(reference)) {
    res <- suppressWarnings(run_si_from_data(
      MockRotavirus$si_data,
      quiet_config(si_parametric_distr = dist, n1 = 500)
    ))
    expect_equal(
      median(res$SI.Moments$Mean), reference[[dist]][["mean"]],
      tolerance = 0.3, scale = 1
    )
    if (dist != "lognormal") {
      expect_equal(
        median(res$SI.Moments$Std), reference[[dist]][["sd"]],
        tolerance = 0.3, scale = 1
      )
    }
    expect_equal(res$R$`Mean(R)`[10], 1.03, tolerance = 0.05, scale = 1)
  }
})

test_that("si_from_data is reproducible with the mcmc_control seed", {
  run <- function(seed) {
    config <- quiet_config(si_parametric_distr = "gamma")
    config$mcmc_control$seed <- seed
    suppressWarnings(run_si_from_data(MockRotavirus$si_data, config))
  }
  expect_identical(run(1)$si_distr, run(1)$si_distr)
  expect_false(identical(run(1)$si_distr, run(3)$si_distr))
})

test_that("si_from_data warns that burnin and thin are ignored", {
  config <- quiet_config(si_parametric_distr = "gamma")
  config$mcmc_control <- make_mcmc_control(burnin = 1000, thin = 5, seed = 1)
  warns <- collect_warnings(run_si_from_data(MockRotavirus$si_data, config))
  expect_true(any(grepl("burnin and thin", warns, fixed = TRUE)))
  config$mcmc_control <- make_mcmc_control(seed = 1)
  warns <- collect_warnings(run_si_from_data(MockRotavirus$si_data, config))
  expect_false(any(grepl("burnin and thin", warns, fixed = TRUE)))
})

test_that("si_from_data uses init_pars as starting values", {
  config <- quiet_config(si_parametric_distr = "gamma")
  default <- suppressWarnings(run_si_from_data(MockRotavirus$si_data, config))
  config$mcmc_control$init_pars <- c(3, 1)
  custom <- suppressWarnings(run_si_from_data(MockRotavirus$si_data, config))
  expect_equal(
    mean(custom$SI.Moments$Mean), mean(default$SI.Moments$Mean),
    tolerance = 1e-3
  )
})

test_that("process_si_data accepts an observation time column", {
  si_data <- MockRotavirus$si_data
  si_data$OT <- 30L
  processed <- process_si_data(si_data)
  expect_identical(processed$OT, si_data$OT)
  si_data$OT <- si_data$SL
  expect_error(process_si_data(si_data), "SL >= OT")
  si_data$OT <- 30
  expect_error(process_si_data(si_data), "OT is non integer")
})

test_that("process_si_data accepts an all NA OT column", {
  si_data <- MockRotavirus$si_data
  si_data$OT <- NA
  expect_identical(process_si_data(si_data)$OT, si_data$OT)
  config <- quiet_config(si_parametric_distr = "gamma")
  with_na <- suppressWarnings(run_si_from_data(si_data, config))
  without <- suppressWarnings(
    run_si_from_data(MockRotavirus$si_data, config)
  )
  expect_identical(with_na$si_distr, without$si_distr)
})

test_that("offset distributions give the offset error before fitting", {
  si_data <- MockRotavirus$si_data
  si_data$SL[1] <- si_data$EL[1]
  si_data$SR[1] <- si_data$EL[1] + 1L
  for (dist in c("gamma_offset_1", "weibull_offset_1", "lognormal_offset_1")) {
    config <- quiet_config(si_parametric_distr = dist)
    expect_error(
      suppressWarnings(run_si_from_data(si_data, config)),
      "offset 1"
    )
  }
})

test_that("process_si_data checks an unnamed fifth column against type", {
  si_data <- MockRotavirus$si_data[, c("EL", "ER", "SL", "SR")]
  si_data$OT <- rep_len(0:2, nrow(si_data))
  names(si_data) <- NULL
  si_data <- as.data.frame(si_data)
  expect_error(suppressWarnings(process_si_data(si_data)), "OT")
  typed <- MockRotavirus$si_data[, c("EL", "ER", "SL", "SR", "type")]
  names(typed) <- NULL
  typed <- as.data.frame(typed)
  expect_warning(processed <- process_si_data(typed), "column names")
  expect_identical(processed$type, MockRotavirus$si_data$type)
})

test_that("process_si_data does not read an unnamed OT column as type", {
  si_data <- MockRotavirus$si_data[, c("EL", "ER", "SL", "SR")]
  si_data$OT <- 30L
  names(si_data) <- NULL
  si_data <- as.data.frame(si_data)
  expect_error(
    suppressWarnings(process_si_data(si_data)),
    "OT"
  )
})

test_that("si_from_data corrects for right truncation with an OT column", {
  mu <- 8
  sigma <- 4
  si_data <- simulate_si_data(
    1000, mu, sigma, obs_time = 40, growth = 0.15, seed = 12
  )
  incid <- MockRotavirus$incidence
  config <- suppressMessages(
    si_from_data_config(incid, si_parametric_distr = "gamma")
  )
  truth <- discr_mean(mu, sigma)

  naive <- suppressWarnings(run_si_from_data(si_data, config, incid))
  ## pairs are observed up to and including day 40
  si_data$OT <- 41L
  truncated <- suppressWarnings(run_si_from_data(si_data, config, incid))
  naive_bias <- mean(naive$SI.Moments$Mean) - truth
  truncated_bias <- mean(truncated$SI.Moments$Mean) - truth
  expect_lt(naive_bias, -0.5)
  expect_lt(abs(truncated_bias), 0.4)
})

test_that("si_from_data supports offset distributions", {
  mu <- 5
  sigma <- 2
  si_data <- simulate_si_data(300, mu - 1, sigma, shift = 1, obs_time = 200)
  incid <- MockRotavirus$incidence
  for (dist in c("gamma_offset_1", "weibull_offset_1", "lognormal_offset_1")) {
    config <- suppressMessages(
      si_from_data_config(incid, si_parametric_distr = dist)
    )
    res <- suppressWarnings(run_si_from_data(si_data, config, incid))
    expect_true(all(res$si_distr[, 1] == 0))
    expect_true(all(res$si_distr[, 2] > 0))
    expect_equal(mean(res$SI.Moments$Mean), mu, tolerance = 0.3, scale = 1)
  }
})

test_that("si_from_data supports other primary event distributions", {
  si_data <- simulate_si_data(300, 5, 2, obs_time = 200)
  incid <- MockRotavirus$incidence
  uniform <- suppressWarnings(run_si_from_data(
    si_data,
    suppressMessages(si_from_data_config(incid, si_parametric_distr = "gamma")),
    incid
  ))
  growth <- suppressWarnings(run_si_from_data(
    si_data,
    suppressMessages(si_from_data_config(
      incid, si_parametric_distr = "gamma",
      si_discr_args = list(
        dprimary = primarycensored::dexpgrowth,
        primary_args = list(r = 0.5)
      )
    )),
    incid
  ))
  expect_false(isTRUE(all.equal(growth$si_distr, uniform$si_distr)))
  expect_equal(rowSums(growth$si_distr), rep(1, 100), tolerance = 1e-8)
})

test_that("si_from_data uses L and D from si_discr_args", {
  config <- quiet_config(
    si_parametric_distr = "gamma", si_discr_args = list(L = 2, D = 5)
  )
  res <- suppressWarnings(run_si_from_data(MockRotavirus$si_data, config))
  k <- seq_len(ncol(res$si_distr)) - 1
  expect_true(all(res$si_distr[, k < 2 | k >= 5] == 0))
  expect_true(all(res$si_distr[, k == 2] > 0))
  expect_equal(rowSums(res$si_distr), rep(1, 100), tolerance = 1e-8)
})

test_that("si_from_data errors if si_discr_args conflicts with the fit", {
  for (args in list(list(dist = stats::plnorm), list(shift = 1))) {
    config <- quiet_config(
      si_parametric_distr = "gamma", si_discr_args = args
    )
    expect_error(
      suppressWarnings(run_si_from_data(MockRotavirus$si_data, config)),
      "si_parametric_distr"
    )
  }
})

test_that("exact dates in si_data map to zero width windows", {
  si_data <- data.frame(
    EL = c(0L, 2L, 5L), ER = c(1L, 2L, 5L),
    SL = c(3L, 6L, 9L), SR = c(4L, 7L, 9L)
  )
  censdata <- si_data_to_censdata(si_data, 0)
  expect_equal(censdata$pwindow, c(1, 0, 0))
  expect_equal(censdata$left, c(3, 4, 4))
  expect_equal(censdata$right, c(4, 5, 4))
  expect_no_message(si_data_to_censdata(si_data, 0))
})

test_that("OT is a continuous upper bound of the observation time", {
  si_data <- data.frame(
    EL = c(0L, 2L), ER = c(1L, 3L), SL = c(3L, 6L), SR = c(4L, 7L),
    OT = c(10L, NA)
  )
  expect_equal(si_data_to_censdata(si_data, 0)$D, c(10, Inf))
  expect_equal(si_data_to_censdata(si_data, 1)$D, c(9, Inf))
})

test_that("si_from_data fits mixed type 0, 1 and 2 rows", {
  si_data <- simulate_si_data(300, 5, 2.5, obs_time = 200, seed = 3)
  exact_primary <- seq(1, 300, by = 3)
  exact_both <- seq(2, 300, by = 3)
  si_data$ER[c(exact_primary, exact_both)] <-
    si_data$EL[c(exact_primary, exact_both)]
  si_data$SR[exact_both] <- si_data$SL[exact_both]
  si_data$type <- NULL
  config <- suppressMessages(
    si_from_data_config(MockRotavirus$incidence, si_parametric_distr = "gamma")
  )
  res <- suppressWarnings(run_si_from_data(si_data, config))
  expect_true(res$si_fit_converged)
  expect_equal(mean(res$SI.Moments$Mean), 5, tolerance = 0.5, scale = 1)
})

test_that("exact rows give the exact and single censored likelihoods", {
  set.seed(8)
  delay <- stats::rgamma(200, shape = 4, scale = 1.2)
  onset <- sample(0:30, 200, replace = TRUE)
  exact <- data.frame(
    EL = onset, ER = onset,
    SL = onset + delay, SR = onset + delay
  )
  fit_exact <- primarycensored::fitdistdoublecens(
    si_data_to_censdata(exact, 0), distr = "gamma",
    start = list(shape = 2, scale = 2)
  )
  reference <- fitdistrplus::fitdist(
    delay, "gamma", start = list(shape = 2, scale = 2)
  )
  expect_equal(fit_exact$estimate, reference$estimate, tolerance = 1e-3)

  single <- data.frame(
    EL = onset, ER = onset,
    SL = onset + floor(delay), SR = onset + floor(delay) + 1
  )
  fit_single <- primarycensored::fitdistdoublecens(
    si_data_to_censdata(single, 0), distr = "gamma",
    start = list(shape = 2, scale = 2)
  )
  reference_single <- fitdistrplus::fitdistcens(
    data.frame(left = floor(delay), right = floor(delay) + 1), "gamma",
    start = list(shape = 2, scale = 2)
  )
  expect_equal(
    fit_single$estimate, reference_single$estimate, tolerance = 1e-3
  )
})

test_that("coarse2estim is deprecated and matches discr_si", {
  samples <- data.frame(shape = c(2, 3), scale = c(1, 1.5))
  expect_warning(
    out <- coarse2estim(dist = "gamma", samples = samples, thin = 1),
    "deprecated"
  )
  k <- seq(0, nrow(out$si_sample) - 1)
  for (i in seq_len(nrow(samples))) {
    expected <- discr_si(
      k, dist = stats::pgamma,
      shape = samples$shape[i], scale = samples$scale[i],
      shift = 0, L = 1
    )
    expect_equal(out$si_sample[, i], expected / sum(expected),
                 tolerance = 1e-10)
  }
  out_off <- suppressWarnings(
    coarse2estim(dist = "gamma_offset_1", samples = samples, thin = 1)
  )
  k <- seq(0, nrow(out_off$si_sample) - 1)
  expected <- discr_si(
    k, dist = stats::pgamma, shape = 2, scale = 1, shift = 1
  )
  expect_equal(out_off$si_sample[, 1], expected / sum(expected),
               tolerance = 1e-10)
})

test_that("check_cdt_samples_convergence is deprecated", {
  set.seed(1)
  samples <- data.frame(a = stats::rnorm(100), b = stats::rnorm(100))
  expect_warning(
    utils::capture.output(check_cdt_samples_convergence(samples)),
    "deprecated"
  )
})

test_that("si_from_data does not warn about right truncation", {
  config <- quiet_config(si_parametric_distr = "gamma")
  warns <- collect_warnings(run_si_from_data(MockRotavirus$si_data, config))
  expect_false(any(grepl("truncation", warns, fixed = TRUE)))
})

test_that("a missing covariance gives an informative error", {
  testthat::local_mocked_bindings(
    fitdistdoublecens = function(...) {
      structure(
        list(estimate = c(shape = 2, scale = 1), vcov = NULL,
             convergence = 0),
        class = "fitdist"
      )
    },
    .package = "primarycensored"
  )
  config <- quiet_config(si_parametric_distr = "gamma")
  expect_error(
    suppressWarnings(run_si_from_data(MockRotavirus$si_data, config)),
    "poor fit"
  )
})

test_that("a non positive definite covariance gives an informative error", {
  expect_error(
    draw_si_params(
      c(shape = 2, scale = 1), matrix(0, 2, 2), 10, c(TRUE, TRUE)
    ),
    "init_pars"
  )
  expect_error(
    draw_si_params(
      c(shape = 2, scale = 1), matrix(NA_real_, 2, 2), 10, c(TRUE, TRUE)
    ),
    "poor fit"
  )
})
