# devtools::load_all()
require(testthat)

### everything needed for the tests in this file ###

data("SARS2003")
incid <- SARS2003$incidence
dt <- 7L
weekly_inc <- aggregate_inc(incid, dt)

mean_si <- sum(SARS2003$si_distr * seq_along(SARS2003$si_distr))
std_si <- sqrt(sum(SARS2003$si_distr * seq_along(SARS2003$si_distr)^2))
si_distr <- SARS2003$si_distr

### ### ### ### ### ### ### ### ### ### ### ### ### 

test_that("aggregating incidence works", {
  expect_error(aggregate_inc(Flu2009$incidence, 7L),
               "incid should be a vector of integer values")
  expect_error(aggregate_inc(SARS2003$incidence, 7),
               "dt should be an integer >=2")
  expect_error(aggregate_inc(SARS2003$incidence, 1L),
               "dt should be an integer >=2")
  expect_error(aggregate_inc(SARS2003$incidence, -1L),
               "dt should be an integer >=2")
  expect_true(aggregate_inc(SARS2003$incidence, 7L)[1] == sum(SARS2003$incidence[1:7]))
})


test_that("weekly version of estimate_R works in parametric mode", {
  
  ## estimate the reproduction number (method "non_parametric_si")
  ## when not specifying t_start and t_end in config, they are set to estimate
  ## the reproduction number on sliding weekly windows                          
  config <- make_config(list(mean_si = mean_si,
                             std_si = std_si))
  res <- estimate_R_agg(incid = weekly_inc, 
                        dt = 7, # aggregation window of the data
                        dt_out = 7, # desired sliding window length
                        iter = 10,
                        config = config,
                        method = "parametric_si",
                        grid = list(precision = 0.001, min = -1, max = 1))
                        
                          
})
  