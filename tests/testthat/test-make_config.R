
test_that("make_config triggers error when given method argument", {
  msg <- paste("`method` should be specified as an argument to",
               "`estimate_R`, not `make_config`.")
  expect_error(
    make_config(method = "uncertain_si"),
    msg)
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
