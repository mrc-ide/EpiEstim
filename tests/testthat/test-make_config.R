
test_that("make_config triggers error when given method argument", {
  msg <- paste("`method` should be specified as an argument to",
               "`estimate_R`, not `make_config`.")
  expect_error(
    make_config(method = "uncertain_si"),
    msg)
})
