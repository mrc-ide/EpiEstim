# devtools::load_all()
require(testthat)

### everything needed for the tests in this file ###

data("Flu2009")
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


test_that("weekly version of estimate_R_agg works in parametric mode", {
  
  ## estimate the reproduction number (method "parametric_si")
  ## when not specifying t_start and t_end in config, they are set to estimate
  ## the reproduction number on sliding weekly windows      
  method <- "parametric_si"
  config <- make_config(list(mean_si = mean_si,
                             std_si = std_si))
  res_weekly <- estimate_R_agg(incid = weekly_inc, 
                        dt = 7, # aggregation window of the data
                        dt_out = 7, # desired sliding window length
                        iter = 10,
                        config = config,
                        method = method,
                        grid = list(precision = 0.001, min = -1, max = 1))
  
  res_daily <- estimate_R(incid = incid, 
                               config = config,
                               method = method)
  
  ######################################################################
  ## test that the weekly incidence matches the aggregated daily one
  ## except for first time window where we impose I = 1 on first day
  ## TODO: explore the issue above of forcing I = 1 on first day
  
  for(i in seq_len(floor(length(res_daily$I) / dt))[-1])
  {
    expect_equal(sum(res_daily$I[i*dt+(1:dt)]), 
                 sum(res_weekly$I[i*dt+(1:dt)]))
  }
  
  ######################################################################
  ## test that the weekly R estimates vaguely match if you use the daily data
  ## or the weekly data
  
  common_t_start <- intersect(res_daily$R$t_start, res_weekly$R$t_start)
  
  relative_error <- abs(res_daily$R$`Mean(R)`[res_daily$R$t_start %in% common_t_start] -   
                          res_weekly$R$`Mean(R)`[res_weekly$R$t_start %in% common_t_start]) / 
    res_daily$R$`Mean(R)`[res_daily$R$t_start %in% common_t_start]

  ## only check after the 10th common_t_start as before hand there is a lot of
  ## day to day variation
  expect_true(all(relative_error[-c(1:20)] < 0.4)) # 0.4 arbitrarily small
  
  ######################################################################
  ## test that the weekly R estimates from weekly data match exactly the weekly
  ## R estimates from the daily reconstructed incidence
  
  res_daily_reconstructed <- estimate_R(incid = res_weekly$I, 
                          config = config,
                          method = method)
  
  common_t_start <- intersect(res_daily_reconstructed$R$t_start, res_weekly$R$t_start)
  
  relative_error <- abs(res_daily_reconstructed$R$`Mean(R)`[res_daily_reconstructed$R$t_start %in% common_t_start] -   
                          res_weekly$R$`Mean(R)`[res_weekly$R$t_start %in% common_t_start]) / 
    res_daily_reconstructed$R$`Mean(R)`[res_daily_reconstructed$R$t_start %in% common_t_start]
  
  expect_true(all(relative_error < 1e-9)) 
                          
})



test_that("weekly version of estimate_R_agg works in  non-parametric mode", {
  
  ## estimate the reproduction number (method "non_parametric_si")
  ## when not specifying t_start and t_end in config, they are set to estimate
  ## the reproduction number on sliding weekly windows     
  method <- "non_parametric_si"
  config <- make_config(list(si_distr = si_distr))
  res_weekly <- estimate_R_agg(incid = weekly_inc, 
                               dt = 7, # aggregation window of the data
                               dt_out = 7, # desired sliding window length
                               iter = 10,
                               config = config,
                               method = method,
                               grid = list(precision = 0.001, min = -1, max = 1))
  
  res_daily <- estimate_R(incid = incid, 
                          config = config,
                          method = method)
  
  ######################################################################
  ## test that the weekly incidence matches the aggregated daily one
  ## except for first time window where we impose I = 1 on first day
  ## TODO: explore the issue above of forcing I = 1 on first day
  
  for(i in seq_len(floor(length(res_daily$I) / dt))[-1])
  {
    expect_equal(sum(res_daily$I[i*dt+(1:dt)]), 
                 sum(res_weekly$I[i*dt+(1:dt)]))
  }
  
  ######################################################################
  ## test that the weekly R estimates vaguely match if you use the daily data
  ## or the weekly data
  
  common_t_start <- intersect(res_daily$R$t_start, res_weekly$R$t_start)
  
  relative_error <- abs(res_daily$R$`Mean(R)`[res_daily$R$t_start %in% common_t_start] -   
                          res_weekly$R$`Mean(R)`[res_weekly$R$t_start %in% common_t_start]) / 
    res_daily$R$`Mean(R)`[res_daily$R$t_start %in% common_t_start]
  
  ## only check after the 10th common_t_start as before hand there is a lot of
  ## day to day variation
  expect_true(all(relative_error[-c(1:20)] < 0.4)) # 0.4 arbitrarily small
  
  ######################################################################
  ## test that the weekly R estimates from weekly data match exactly the weekly
  ## R estimates from the daily reconstructed incidence
  
  res_daily_reconstructed <- estimate_R(incid = res_weekly$I, 
                                        config = config,
                                        method = method)
  
  common_t_start <- intersect(res_daily_reconstructed$R$t_start, res_weekly$R$t_start)
  
  relative_error <- abs(res_daily_reconstructed$R$`Mean(R)`[res_daily_reconstructed$R$t_start %in% common_t_start] -   
                          res_weekly$R$`Mean(R)`[res_weekly$R$t_start %in% common_t_start]) / 
    res_daily_reconstructed$R$`Mean(R)`[res_daily_reconstructed$R$t_start %in% common_t_start]
  
  expect_true(all(relative_error < 1e-9)) 
  
})
  

test_that("weekly version of estimate_R works with aggregated data in parametric mode", {
  
  ## estimate the reproduction number (method "parametric_si")
  ## when not specifying t_start and t_end in config, they are set to estimate
  ## the reproduction number on sliding weekly windows      
  method <- "parametric_si"
  config <- make_config(list(mean_si = mean_si,
                             std_si = std_si))
  res_weekly <- estimate_R(incid = weekly_inc, 
                               dt = 7, # aggregation window of the data
                               dt_out = 7, # desired sliding window length
                               iter = 10,
                               config = config,
                               method = method,
                               grid = list(precision = 0.001, min = -1, max = 1))
  
  res_weekly_agg <- estimate_R_agg(incid = weekly_inc, 
                           dt = 7, # aggregation window of the data
                           dt_out = 7, # desired sliding window length
                           iter = 10,
                           config = config,
                           method = method,
                           grid = list(precision = 0.001, min = -1, max = 1))
  
  res_daily <- estimate_R(incid = incid, 
                          config = config,
                          method = method)
  
  ######################################################################
  ## test that the weekly incidence matches the aggregated daily one
  ## except for first time window where we impose I = 1 on first day
  ## TODO: explore the issue above of forcing I = 1 on first day
  
  for(i in seq_len(floor(length(res_daily$I) / dt))[-1])
  {
    expect_equal(sum(res_daily$I[i*dt+(1:dt)]), 
                 sum(res_weekly$I[i*dt+(1:dt)]))
  }
  
  ######################################################################
  ## test that the weekly R estimates vaguely match if you use the daily data
  ## or the weekly data
  
  common_t_start <- intersect(res_daily$R$t_start, res_weekly$R$t_start)
  
  relative_error <- abs(res_daily$R$`Mean(R)`[res_daily$R$t_start %in% common_t_start] -   
                          res_weekly$R$`Mean(R)`[res_weekly$R$t_start %in% common_t_start]) / 
    res_daily$R$`Mean(R)`[res_daily$R$t_start %in% common_t_start]
  
  ## only check after the 10th common_t_start as before hand there is a lot of
  ## day to day variation
  expect_true(all(relative_error[-c(1:20)] < 0.4)) # 0.4 arbitrarily small
  
  ######################################################################
  ## test that the weekly R estimates from weekly data match exactly the weekly
  ## R estimates from the daily reconstructed incidence
  
  res_daily_reconstructed <- estimate_R(incid = res_weekly$I, 
                                        config = config,
                                        method = method)
  
  common_t_start <- intersect(res_daily_reconstructed$R$t_start, res_weekly$R$t_start)
  
  relative_error <- abs(res_daily_reconstructed$R$`Mean(R)`[res_daily_reconstructed$R$t_start %in% common_t_start] -   
                          res_weekly$R$`Mean(R)`[res_weekly$R$t_start %in% common_t_start]) / 
    res_daily_reconstructed$R$`Mean(R)`[res_daily_reconstructed$R$t_start %in% common_t_start]
  
  expect_true(all(relative_error < 1e-9)) 
  
  ######################################################################
  ## test that the weekly R estimates from estimate_R are the same as from 
  ## estimate_R_agg
  
  expect_true(all(res_weekly$R == res_weekly_agg$R)) 
  
})



test_that("weekly version of estimate_R works with aggregated data in non-parametric mode", {
  
  ## estimate the reproduction number (method "non_parametric_si")
  ## when not specifying t_start and t_end in config, they are set to estimate
  ## the reproduction number on sliding weekly windows      
  method <- "non_parametric_si"
  config <- make_config(list(si_distr = si_distr))
  res_weekly <- estimate_R(incid = weekly_inc, 
                           dt = 7, # aggregation window of the data
                           dt_out = 7, # desired sliding window length
                           iter = 10,
                           config = config,
                           method = method,
                           grid = list(precision = 0.001, min = -1, max = 1))
  
  res_weekly_agg <- estimate_R_agg(incid = weekly_inc, 
                                   dt = 7, # aggregation window of the data
                                   dt_out = 7, # desired sliding window length
                                   iter = 10,
                                   config = config,
                                   method = method,
                                   grid = list(precision = 0.001, min = -1, max = 1))
  
  res_daily <- estimate_R(incid = incid, 
                          config = config,
                          method = method)
  
  ######################################################################
  ## test that the weekly incidence matches the aggregated daily one
  ## except for first time window where we impose I = 1 on first day
  ## TODO: explore the issue above of forcing I = 1 on first day
  
  for(i in seq_len(floor(length(res_daily$I) / dt))[-1])
  {
    expect_equal(sum(res_daily$I[i*dt+(1:dt)]), 
                 sum(res_weekly$I[i*dt+(1:dt)]))
  }
  
  ######################################################################
  ## test that the weekly R estimates vaguely match if you use the daily data
  ## or the weekly data
  
  common_t_start <- intersect(res_daily$R$t_start, res_weekly$R$t_start)
  
  relative_error <- abs(res_daily$R$`Mean(R)`[res_daily$R$t_start %in% common_t_start] -   
                          res_weekly$R$`Mean(R)`[res_weekly$R$t_start %in% common_t_start]) / 
    res_daily$R$`Mean(R)`[res_daily$R$t_start %in% common_t_start]
  
  ## only check after the 10th common_t_start as before hand there is a lot of
  ## day to day variation
  expect_true(all(relative_error[-c(1:20)] < 0.4)) # 0.4 arbitrarily small
  
  ######################################################################
  ## test that the weekly R estimates from weekly data match exactly the weekly
  ## R estimates from the daily reconstructed incidence
  
  res_daily_reconstructed <- estimate_R(incid = res_weekly$I, 
                                        config = config,
                                        method = method)
  
  common_t_start <- intersect(res_daily_reconstructed$R$t_start, res_weekly$R$t_start)
  
  relative_error <- abs(res_daily_reconstructed$R$`Mean(R)`[res_daily_reconstructed$R$t_start %in% common_t_start] -   
                          res_weekly$R$`Mean(R)`[res_weekly$R$t_start %in% common_t_start]) / 
    res_daily_reconstructed$R$`Mean(R)`[res_daily_reconstructed$R$t_start %in% common_t_start]
  
  expect_true(all(relative_error < 1e-9)) 
  
  ######################################################################
  ## test that the weekly R estimates from estimate_R are the same as from 
  ## estimate_R_agg
  
  expect_true(all(res_weekly$R == res_weekly_agg$R)) 
  
})

## TODO: 
## check that if incid is not a vector you can't run with dt > 1 (if it's an incidence object or a matrix / dataframe)
## check that if method is not parametric or non_parametric you can't run with dt > 1
## check that if you don't supply the right arguments it does not work (e.g. no mean_SI for a parametric method)
## do a test with a different dt (3, 10)



