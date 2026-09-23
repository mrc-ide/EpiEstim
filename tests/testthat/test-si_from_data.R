data("MockRotavirus")

quiet_config <- function(...) {
  suppressMessages(si_from_data_config(MockRotavirus$incidence, ...))
}

test_that("si_from_data keeps the output structure", {
  res <- suppressWarnings(run_si_from_data(
    MockRotavirus$si_data, quiet_config(si_parametric_distr = "gamma")
  ))
  expect_s3_class(res, "estimate_R")
  expect_named(res, c(
    "R", "method", "si_distr", "SI.Moments", "dates", "I", "I_local",
    "I_imported", "MCMC_converged"
  ))
  expect_identical(res$method, "si_from_data")
  expect_identical(nrow(res$si_distr), 100L)
  expect_true(all(res$si_distr[, 1] == 0))
  expect_equal(rowSums(res$si_distr), rep(1, 100), tolerance = 1e-8)
  expect_named(res$SI.Moments, c("Mean", "Std"))
  expect_identical(nrow(res$SI.Moments), 100L)
  expect_true(res$MCMC_converged)
})

test_that("si_from_data is close to the previous coarseDataTools results", {
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
  si_data$OT <- si_data$SL - 1L
  expect_error(process_si_data(si_data), "SL > OT")
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

  warns <- collect_warnings(
    naive <- run_si_from_data(si_data, config, incid)
  )
  expect_true(any(grepl("right truncation", warns, fixed = TRUE)))
  si_data$OT <- 40L
  warns <- collect_warnings(
    truncated <- run_si_from_data(si_data, config, incid)
  )
  expect_false(any(grepl("right truncation", warns, fixed = TRUE)))
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

test_that("si_from_data reads exact dates as one day intervals", {
  si_data <- MockRotavirus$si_data
  si_data$ER[1:3] <- si_data$EL[1:3]
  si_data$type <- NULL
  config <- quiet_config(si_parametric_distr = "gamma")
  expect_message(
    suppressWarnings(run_si_from_data(si_data, config)),
    "one day"
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
