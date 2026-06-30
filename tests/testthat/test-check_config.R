test_that("check_config validates parametric_si inputs", {
  expect_silent(
    EpiEstim:::check_config(
      list(
        t_start = 2:3,
        t_end = 4:5,
        mean_si = 2.5,
        std_si = 1,
        cv_posterior = 0.3
      ),
      "parametric_si",
      n_time_steps = 10
    )
  )

  expect_error(
    EpiEstim:::check_config(
      list(mean_si = 0.9, std_si = 1, cv_posterior = 0.3),
      "parametric_si"
    ),
    "method parametric_si requires a value >1 for config\\$mean_si."
  )

  expect_error(
    EpiEstim:::check_config(
      list(mean_si = 2.5, std_si = 0, cv_posterior = 0.3),
      "parametric_si"
    ),
    "method parametric_si requires a >0 value for config\\$std_si."
  )
})

test_that("check_config validates non_parametric_si inputs and time windows", {
  expect_silent(
    EpiEstim:::check_config(
      list(
        t_start = 2:3,
        t_end = 4:5,
        si_distr = c(0, 1),
        cv_posterior = 0.3
      ),
      "non_parametric_si",
      n_time_steps = 10
    )
  )

  expect_error(
    EpiEstim:::check_config(
      list(si_distr = c(0.2, 0.8), cv_posterior = 0.3),
      "non_parametric_si"
    ),
    "si_distr should be so that si_distr\\[1\\] = 0."
  )

  expect_error(
    EpiEstim:::check_config(
      list(
        t_start = c(2, 4),
        t_end = c(3, 11),
        si_distr = c(0, 1),
        cv_posterior = 0.3
      ),
      "non_parametric_si",
      n_time_steps = 10
    ),
    "t_end must be a vector of integers between 2 and the number of\\s+timesteps in incid."
  )
})

test_that("check_config validates si_from_sample inputs", {
  expect_error(
    EpiEstim:::check_config(
      list(
        t_start = 2:3,
        t_end = 4:5,
        cv_posterior = 0.3
      ),
      "si_from_sample",
      n_time_steps = 10
    ),
    "method si_from_sample requires to specify the config\\$n2 argument."
  )

  expect_silent(
    EpiEstim:::check_config(
      list(
        t_start = 2:3,
        t_end = 4:5,
        n2 = 10,
        cv_posterior = 0.3
      ),
      "si_from_sample",
      n_time_steps = 10
    )
  )
})
