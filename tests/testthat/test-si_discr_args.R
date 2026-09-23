data("Flu2009")

flu_t <- nrow(Flu2009$incidence)
lnorm_args <- list(dist = stats::plnorm)
lnorm_si <- discr_si(seq(0, flu_t - 1), 2.6, 1.5, dist = stats::plnorm)

test_that("make_config defaults si_discr_args to an empty list", {
  expect_identical(make_config()$si_discr_args, list())
})

test_that("estimate_R parametric_si is unchanged with default si_discr_args", {
  config <- list(t_start = 2:26, t_end = 8:32, mean_si = 2.6, std_si = 1.5)
  default <- estimate_R(
    Flu2009$incidence, method = "parametric_si", config = make_config(config)
  )
  explicit <- estimate_R(
    Flu2009$incidence, method = "parametric_si",
    config = make_config(c(config, list(si_discr_args = list())))
  )
  expect_identical(default$R, explicit$R)
  expect_identical(default$si_distr, explicit$si_distr)
})

test_that("estimate_R parametric_si passes si_discr_args to discr_si", {
  config <- list(t_start = 2:26, t_end = 8:32)
  parametric <- estimate_R(
    Flu2009$incidence, method = "parametric_si",
    config = make_config(c(
      config, list(mean_si = 2.6, std_si = 1.5, si_discr_args = lnorm_args)
    ))
  )
  non_parametric <- estimate_R(
    Flu2009$incidence, method = "non_parametric_si",
    config = make_config(c(config, list(si_distr = lnorm_si)))
  )
  expect_equal(parametric$R, non_parametric$R)
})

test_that("estimate_R uncertain_si passes si_discr_args to discr_si", {
  config <- list(
    t_start = 2:26, t_end = 8:32,
    mean_si = 2.6, std_mean_si = 1, min_mean_si = 1, max_mean_si = 4.2,
    std_si = 1.5, std_std_si = 0.5, min_std_si = 0.5, max_std_si = 2.5,
    n1 = 20, n2 = 20, seed = 1
  )
  gamma <- estimate_R(
    Flu2009$incidence, method = "uncertain_si", config = make_config(config)
  )
  lnorm <- estimate_R(
    Flu2009$incidence, method = "uncertain_si",
    config = make_config(c(config, list(si_discr_args = lnorm_args)))
  )
  expect_false(isTRUE(all.equal(gamma$si_distr, lnorm$si_distr)))
  expect_true(all(lnorm$si_distr[, 1] == 0))
})

test_that("wallinga_teunis parametric_si passes si_discr_args to discr_si", {
  config <- list(t_start = 2:26, t_end = 8:32, n_sim = 5)
  set.seed(1)
  parametric <- wallinga_teunis(
    Flu2009$incidence, method = "parametric_si",
    config = c(
      config, list(mean_si = 2.6, std_si = 1.5, si_discr_args = lnorm_args)
    )
  )
  set.seed(1)
  non_parametric <- wallinga_teunis(
    Flu2009$incidence, method = "non_parametric_si",
    config = c(config, list(si_distr = lnorm_si))
  )
  expect_equal(parametric$R, non_parametric$R)
})

test_that("estimate_R_agg passes si_discr_args to discr_si", {
  weekly_inc <- c(
    sum(Flu2009$incidence$I[1:7]), sum(Flu2009$incidence$I[8:14]),
    sum(Flu2009$incidence$I[15:21]), sum(Flu2009$incidence$I[22:28])
  )
  run_agg <- function(si_discr_args) {
    suppressWarnings(estimate_R_agg(
      incid = weekly_inc, dt = 7L, dt_out = 7L, iter = 5L,
      config = make_config(list(
        mean_si = 2.6, std_si = 1.5, si_discr_args = si_discr_args
      )),
      method = "parametric_si",
      grid = list(precision = 0.001, min = -1, max = 1)
    ))
  }
  gamma <- run_agg(list())
  lnorm <- run_agg(lnorm_args)
  expect_false(isTRUE(all.equal(gamma$R, lnorm$R)))
})

test_that("si_discr_args is validated", {
  config <- list(t_start = 2:26, t_end = 8:32, mean_si = 2.6, std_si = 1.5)
  run <- function(si_discr_args) {
    estimate_R(
      Flu2009$incidence, method = "parametric_si",
      config = make_config(c(config, list(si_discr_args = si_discr_args)))
    )
  }
  expect_error(run(list(mu = 3)), "si_discr_args")
  expect_error(run(list(foo = 1)), "si_discr_args")
  expect_error(run(list(dist = stats::pweibull)), "si_discr_args")
  expect_error(run(list(dist = "lognormal")), "si_discr_args")
  expect_error(run("lognormal"), "si_discr_args")
  expect_error(run(list(shift = 0)), "serial interval of zero")
  expect_no_error(run(list(shift = 0, L = 1)))
})

test_that("estimate_R_agg truncates the serial interval at the series end", {
  captured <- new.env()
  local_mocked_bindings(
    discr_si_config = function(k, mu, sigma, si_discr_args = NULL) {
      # estimate_R_agg discretises over the whole series, 0 to 28 days
      if (length(k) == 29) {
        captured$args <- si_discr_args
      }
      discr_si(k, mu, sigma)
    }
  )
  weekly_inc <- c(20, 40, 80, 60)
  suppressWarnings(estimate_R_agg(
    incid = weekly_inc, dt = 7L, dt_out = 7L, iter = 2L,
    config = make_config(list(
      mean_si = 2.6, std_si = 1.5, si_discr_args = list(D = 100)
    )),
    method = "parametric_si",
    grid = list(precision = 0.001, min = -1, max = 1)
  ))
  expect_identical(captured$args$D, 29)
})

test_that("si_discr_args with an unnamed element is rejected", {
  # make_config drops unnamed elements, but wallinga_teunis takes a plain list
  expect_error(check_si_discr_args(list(shift = 1, 5)), "named list")
})
