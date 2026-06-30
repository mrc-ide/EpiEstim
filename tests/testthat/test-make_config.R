
test_that("make_config triggers error when given method argument", {
  msg <- paste("`method` should be specified as an argument to",
               "`estimate_R`, not `make_config`.")
  expect_error(
    make_config(method = "uncertain_si"),
    msg)
})

test_that("make_config accepts valid direct arguments", {
  expect_silent({
    config <- make_config(
      mean_si = 2.6,
      std_si = 1.5,
      n1 = 250,
      n2 = 75,
      mean_prior = 4,
      std_prior = 6,
      cv_posterior = 0.2
    )
  })

  expect_s3_class(config, "estimate_R_config")
  expect_equal(config$mean_si, 2.6)
  expect_equal(config$std_si, 1.5)
  expect_equal(config$n1, 250)
  expect_equal(config$n2, 75)
  expect_equal(config$mean_prior, 4)
  expect_equal(config$std_prior, 6)
  expect_equal(config$cv_posterior, 0.2)
})

test_that("make_config accepts valid incidence and time windows", {
  incid <- data.frame(
    local = c(0, 1, 2, 1, 0, 1, 2, 1),
    imported = rep(0, 8)
  )

  expect_silent({
    config <- make_config(
      incid = incid,
      t_start = c(2, 4),
      t_end = c(5, 7)
    )
  })

  expect_s3_class(config, "estimate_R_config")
  expect_equal(config$t_start, c(2, 4))
  expect_equal(config$t_end, c(5, 7))
})

test_that("make_config triggers process_I error for invalid incid columns", {
  incid <- data.frame(foo = 1:5, bar = 6:10)
  expect_error(
    make_config(incid = incid),
    "incid must be a vector or a dataframe with either i\\) a column")
})

test_that("make_config triggers process_I warning when incid$local[1] > 0", {
  incid <- data.frame(local = c(3, 1, 0, 2), imported = c(0, 0, 0, 0))
  expect_warning(
    make_config(incid = incid),
    "incid\\$local\\[1\\] is >0 but must be 0")
})

test_that("make_config triggers process_I error for negative counts with dates", {
  incid <- data.frame(
    dates = 1:4,
    local = c(0, -1, 0, 2),
    imported = c(0, 0, 0, 0)
  )
  expect_error(
    make_config(incid = incid),
    "incid must contain only non negative integer values\\.")

  incid_no_dates <- incid[, c("local", "imported")]
  expect_error(
    make_config(incid = incid_no_dates),
    "incid must contain only non negative integer values\\.")
})

test_that("make_config triggers errors through check_times ", {
  ## error for mismatched window lengths
  incid <- data.frame(local = c(0, 1, 2, 1, 0, 1, 2, 1), imported = 0)
  expect_error(
    make_config(incid = incid, t_start = c(2, 3), t_end = c(6)),
    "t_start and t_end must have the same length\\.")


  ## error when t_start exceeds t_end"
  incid <- data.frame(local = c(0, 1, 2, 1, 0, 1, 2, 1), imported = 0)
  expect_error(
    make_config(incid = incid, t_start = c(2, 7), t_end = c(6, 6)),
    "t_start\\[i\\] must be <= t_end\\[i\\] for all i\\.")


  ## error for invalid t_start values"
  incid <- data.frame(local = c(0, 1, 2, 1, 0, 1, 2, 1), imported = 0)
  expect_error(
    make_config(incid = incid, t_start = c(1, 3), t_end = c(6, 7)),
    "t_start must be a vector of integers between 2 and the number of")


  ## error for invalid t_end values", 
  incid <- data.frame(local = c(0, 1, 2, 1, 0, 1, 2, 1), imported = 0)
  expect_error(
    make_config(incid = incid, t_start = c(2, 3), t_end = c(6, 8.5)),
    "t_end must be a vector of integers between 2 and the number of")
})
