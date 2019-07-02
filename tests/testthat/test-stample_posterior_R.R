
data("Flu2009")

## estimate the reproduction number (method "non_parametric_si")
## when not specifying t_start and t_end in config, they are set to estimate
## the reproduction number on sliding weekly windows                          
res <- estimate_R(incid = Flu2009$incidence, 
                  method = "non_parametric_si",
                  config = make_config(list(si_distr = Flu2009$si_distr)))


test_that("only estimate_R objects are accepted", {
                  
  expect_error(sample_posterior_R(res$R), "R must be generated from the estimate_r function")                  

})

test_that("sampling returns n samples of R", {

  set.seed(2019)            
  R_sample <- sample_posterior_R(res, n = 1001, window = 1L)
  expect_length(R_sample, 1001)

})


test_that("different windows can be specified", {

  set.seed(2019)
  R_sample1 <- sample_posterior_R(res, n = 1001, window = 1L)

  set.seed(2019)
  R_sample2 <- sample_posterior_R(res, n = 1001, window = 1L)

  set.seed(2019)
  R_sample3 <- sample_posterior_R(res, n = 1001, window = 21L)

  # Samples from the same window should be the same
  expect_equal(R_sample1, R_sample2)
  
  # Samples from a different window should be different
  expect_failure(expect_equal(R_sample1, R_sample3))

})
