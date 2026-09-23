# Check that R estimates are unchanged from those using the original
# closed form discr_si (old_discr_si in helper-discr_si.R).
data("Flu2009")

flu_t <- nrow(Flu2009$incidence)
old_flu_si <- old_discr_si(seq(0, flu_t - 1), 2.6, 1.5)

test_that("estimate_R parametric_si matches the original discr_si", {
  config <- list(t_start = 2:26, t_end = 8:32)
  new <- estimate_R(
    Flu2009$incidence, method = "parametric_si",
    config = make_config(c(config, list(mean_si = 2.6, std_si = 1.5)))
  )
  old <- estimate_R(
    Flu2009$incidence, method = "non_parametric_si",
    config = make_config(c(config, list(si_distr = old_flu_si)))
  )
  expect_equal(new$R, old$R, tolerance = 1e-10)
})

test_that("estimate_R uncertain_si matches the original discr_si", {
  config <- make_config(list(
    t_start = 2:26, t_end = 8:32,
    mean_si = 2.6, std_mean_si = 1, min_mean_si = 1, max_mean_si = 4.2,
    std_si = 1.5, std_std_si = 0.5, min_std_si = 0.5, max_std_si = 2.5,
    n1 = 20, n2 = 20, seed = 1
  ))
  new <- estimate_R(
    Flu2009$incidence, method = "uncertain_si", config = config
  )
  local_mocked_bindings(discr_si_draws = old_discr_si_draws)
  old <- estimate_R(
    Flu2009$incidence, method = "uncertain_si", config = config
  )
  expect_equal(new$R, old$R, tolerance = 1e-10)
  expect_equal(new$si_distr, old$si_distr, tolerance = 1e-10)
})

test_that("wallinga_teunis parametric_si matches the original discr_si", {
  config <- list(t_start = 2:26, t_end = 8:32, n_sim = 5)
  set.seed(1)
  new <- wallinga_teunis(
    Flu2009$incidence, method = "parametric_si",
    config = c(config, list(mean_si = 2.6, std_si = 1.5))
  )
  set.seed(1)
  old <- wallinga_teunis(
    Flu2009$incidence, method = "non_parametric_si",
    config = c(config, list(si_distr = old_flu_si))
  )
  expect_equal(new$R, old$R, tolerance = 1e-10)
})

test_that("estimate_R_agg parametric_si matches the original discr_si", {
  weekly_inc <- c(
    sum(Flu2009$incidence$I[1:7]), sum(Flu2009$incidence$I[8:14]),
    sum(Flu2009$incidence$I[15:21]), sum(Flu2009$incidence$I[22:28])
  )
  run_agg <- function() {
    suppressWarnings(estimate_R_agg(
      incid = weekly_inc, dt = 7L, dt_out = 7L, iter = 5L,
      config = make_config(list(mean_si = 2.6, std_si = 1.5)),
      method = "parametric_si",
      grid = list(precision = 0.001, min = -1, max = 1)
    ))
  }
  new <- run_agg()
  # The serial interval is truncated at the end of the aggregated series,
  # which normalises it over k
  local_mocked_bindings(
    discr_si_config = function(k, mu, sigma, si_discr_args = NULL) {
      w <- old_discr_si(k, mu, sigma)
      w / sum(w)
    }
  )
  old <- run_agg()
  # The grid search in epitrix::r2R0 amplifies floating point differences
  expect_equal(new$R, old$R, tolerance = 1e-6)
})
