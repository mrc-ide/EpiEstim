
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
