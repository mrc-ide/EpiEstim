data("Flu2009")

test_that("There is an informative error if dates are not ascending", {
  inc_shuffled_dates <- Flu2009$incidence[sample(nrow(Flu2009$incidence)), ]
  expect_error(estimate_R(inc_shuffled_dates, method = "non_parametric_si",
                          config = list(t_start = 2:26, t_end = 8:32,
                                        si_distr = Flu2009$si_distr, 
                                        seed = 1)),
               "dates in incid must be in ascending order")
})


# The following examples are from EpiEstim's documentation.

test_that("Example 1 matches saved output", {
  testthat::local_edition(3)
  out <- estimate_R(Flu2009$incidence, method = "non_parametric_si",
                    config = list(t_start = 2:26, t_end = 8:32,
                                  si_distr = Flu2009$si_distr, 
                                  seed = 1))
  expect_snapshot(out)
})

test_that("Example 2 matches saved output", {
  testthat::local_edition(3)
  data <- c(0, 1, 1, 2, 1, 3, 4, 5, 5, 5, 5, 4, 4, 26, 6, 7, 9)
  location <- c("imported", "local", "imported", "imported", "local",
                "imported", "imported", "imported", "imported",
                "local", "local", "local", "imported", "local",
                "imported", "local", "imported")
  incid <- incidence::incidence(data, groups = location)
  out <- estimate_R(incid, method = "parametric_si",
                    config = list(t_start = 2:21, 
                                  t_end = 8:27,
                                  mean_si = 2.6, 
                                  std_si = 1.5, 
                                  seed = 1))
  expect_snapshot(out)
})

test_that("Example 3 matches saved output", {
  testthat::local_edition(3)
  out <- estimate_R(Flu2009$incidence, method = "parametric_si",
                    config = list(t_start = 2:26, 
                                  t_end = 8:32,
                                  mean_si = 2.6, 
                                  std_si = 1.5, 
                                  seed = 1))
  expect_snapshot(out)
})

test_that("Example 4 matches saved output", {
  testthat::local_edition(3)
  out <- estimate_R(Flu2009$incidence, method = "uncertain_si",
                    config = list(t_start = 2:26, 
                                  t_end = 8:32,
                                  mean_si = 2.6, 
                                  std_mean_si = 1,
                                  min_mean_si = 1,
                                  max_mean_si = 4.2,
                                  std_si = 1.5,
                                  std_std_si = 0.5,
                                  min_std_si = 0.5, 
                                  max_std_si = 2.5,
                                  n1 = 100, 
                                  n2 = 100, 
                                  seed = 1))
  expect_snapshot(out)
})
