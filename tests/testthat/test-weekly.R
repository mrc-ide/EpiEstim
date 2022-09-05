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

test_that("function to aggregate incidence works", {
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


test_that("estimate_R_agg works in parametric mode", {
  ## estimate the reproduction number (method "parametric_si")
  ## when not specifying t_start and t_end in config, they are set to estimate
  ## the reproduction number on sliding weekly windows      
  method <- "parametric_si"
  config <- make_config(list(mean_si = mean_si,
                             std_si = std_si))
  res_weekly <- suppressWarnings(estimate_R_agg(incid = weekly_inc, 
                        dt = 7L, # aggregation window of the data
                        dt_out = 7L, # desired sliding window length
                        iter = 10L,
                        config = config,
                        method = method,
                        grid = list(precision = 0.001, min = -1, max = 1)))
  
  
  res_daily <- suppressWarnings(estimate_R(incid = incid, 
                               config = config,
                               method = method))
  
  ######################################################################
  ## test that the weekly incidence matches the aggregated daily one
  ## except for first time window where we impose I = 1 on first day
  
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
  
  res_daily_reconstructed <- suppressWarnings(estimate_R(incid = res_weekly$I, 
                          config = config,
                          method = method))
  
  common_t_start <- intersect(res_daily_reconstructed$R$t_start, res_weekly$R$t_start)
  
  relative_error <- abs(res_daily_reconstructed$R$`Mean(R)`[res_daily_reconstructed$R$t_start %in% common_t_start] -   
                          res_weekly$R$`Mean(R)`[res_weekly$R$t_start %in% common_t_start]) / 
    res_daily_reconstructed$R$`Mean(R)`[res_daily_reconstructed$R$t_start %in% common_t_start]
  
  expect_true(all(relative_error < 1e-9)) 
                          
})



test_that("estimate_R_agg works in  non-parametric mode", {
  ## estimate the reproduction number (method "non_parametric_si")
  ## when not specifying t_start and t_end in config, they are set to estimate
  ## the reproduction number on sliding weekly windows     
  method <- "non_parametric_si"
  config <- make_config(list(si_distr = si_distr))
  res_weekly <- suppressWarnings(estimate_R_agg(incid = weekly_inc, 
                               dt = 7L, 
                               dt_out = 7L, 
                               iter = 10L,
                               config = config,
                               method = method,
                               grid = list(precision = 0.001, min = -1, max = 1)))
  
  res_daily <- suppressWarnings(estimate_R(incid = incid, 
                          config = config,
                          method = method))
  
  ######################################################################
  ## test that the weekly incidence matches the aggregated daily one
  ## except for first time window where we impose I = 1 on first day
  
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
  
  res_daily_reconstructed <- suppressWarnings(estimate_R(incid = res_weekly$I, 
                                        config = config,
                                        method = method))
  
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
  res_weekly <- suppressWarnings(estimate_R(incid = weekly_inc, 
                               dt = 7L, 
                               dt_out = 7L, 
                               iter = 10L,
                               config = config,
                               method = method,
                               grid = list(precision = 0.001, min = -1, max = 1)))
  
  res_weekly_agg <- suppressWarnings(estimate_R_agg(incid = weekly_inc, 
                           dt = 7L, 
                           dt_out = 7L, 
                           iter = 10L,
                           config = config,
                           method = method,
                           grid = list(precision = 0.001, min = -1, max = 1)))
  
  res_daily <- suppressWarnings(estimate_R(incid = incid, 
                          config = config,
                          method = method))
  
  ######################################################################
  ## test that the weekly incidence matches the aggregated daily one
  ## except for first time window where we impose I = 1 on first day
  
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
  
  res_daily_reconstructed <- suppressWarnings(estimate_R(incid = res_weekly$I, 
                                        config = config,
                                        method = method))
  
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
  res_weekly <- suppressWarnings(estimate_R(incid = weekly_inc, 
                           dt = 7L, 
                           dt_out = 7L, 
                           iter = 10L,
                           config = config,
                           method = method,
                           grid = list(precision = 0.001, min = -1, max = 1)))
  
  res_weekly_agg <- suppressWarnings(estimate_R_agg(incid = weekly_inc, 
                                   dt = 7L, 
                                   dt_out = 7L, 
                                   iter = 10L,
                                   config = config,
                                   method = method,
                                   grid = list(precision = 0.001, min = -1, max = 1)))
  
  res_daily <- suppressWarnings(estimate_R(incid = incid, 
                          config = config,
                          method = method))
  
  ######################################################################
  ## test that the weekly incidence matches the aggregated daily one
  ## except for first time window where we impose I = 1 on first day
  
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
  
  res_daily_reconstructed <- suppressWarnings(estimate_R(incid = res_weekly$I, 
                                        config = config,
                                        method = method))
  
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



test_that("r grid can be automatically updated with similar results", {
  ## estimate the reproduction number (method "parametric_si")
  ## when not specifying t_start and t_end in config, they are set to estimate
  ## the reproduction number on sliding weekly windows      
  method <- "parametric_si"
  config <- make_config(list(mean_si = mean_si,
                             std_si = std_si))
  
 res_weekly_small_grid <- suppressWarnings(estimate_R_agg(incid = weekly_inc, 
                               dt = 7L, 
                               dt_out = 7L, 
                               iter = 10L,
                               config = config,
                               method = method,
                               grid = list(precision = 0.001, min = 0.01, max = 0.012)))
 
 res_weekly_default_grid <- suppressWarnings(estimate_R_agg(incid = weekly_inc, 
                                         dt = 7L, 
                                         dt_out = 7L, 
                                         iter = 10L,
                                         config = config,
                                         method = method))
 
 expect_true(max(abs(res_weekly_small_grid$R - res_weekly_default_grid$R)) < 1e-9)
  
})



test_that("dt and dt_out in estimate_R_agg are in the correct format", {
  method <- "parametric_si"
  config <- make_config(list(mean_si = mean_si,
                             std_si = std_si))
  
  expect_error(suppressWarnings(estimate_R_agg(incid = weekly_inc, 
                                               dt = 7, 
                                               dt_out = 7L, 
                                               iter = 10L,
                                               config = config,
                                               method = method)),
               "dt and dt_out must be integers e.g. 7L") 
  
  expect_error(suppressWarnings(estimate_R_agg(incid = weekly_inc, 
                                               dt = 7L, 
                                               dt_out = 7, 
                                               iter = 10L,
                                               config = config,
                                               method = method)),
               "dt and dt_out must be integers e.g. 7L") 
  
})



test_that("iter in estimate_R_agg is in the correct format", {
  method <- "parametric_si"
  config <- make_config(list(mean_si = mean_si,
                             std_si = std_si))
  
  expect_error(suppressWarnings(estimate_R_agg(incid = weekly_inc, 
                                               dt = 7L, 
                                               dt_out = 7L, 
                                               iter = 10.5,
                                               config = config,
                                               method = method)),
               "iter must be an integer e.g. 10L") 
  
})


test_that("grid in estimate_R_agg is in the correct format", {
  method <- "parametric_si"
  config <- make_config(list(mean_si = mean_si,
                             std_si = std_si))
  
  expect_error(suppressWarnings(estimate_R_agg(incid = weekly_inc, 
                                               dt = 7L, 
                                               dt_out = 7L, 
                                               iter = 10L,
                                               config = config,
                                               method = method,
                                               grid = c(precision = 0.001, min = -1, max=1))),
               "grid must be a list of 3 elements: precision, min, and max")
  
  expect_error(suppressWarnings(estimate_R_agg(incid = weekly_inc, 
                                               dt = 7L, 
                                               dt_out = 7L, 
                                               iter = 10L,
                                               config = config,
                                               method = method,
                                               grid = list(precision = 0.001, min = -1))),
               "grid must be a list of 3 elements: precision, min, and max")
  
  expect_error(suppressWarnings(estimate_R_agg(incid = weekly_inc, 
                                               dt = 7L, 
                                               dt_out = 7L, 
                                               iter = 10L,
                                               config = config,
                                               method = method,
                                               grid = list(precision = 0.001, min = -1, max = -2))),
               "grid max must be larger than grid min") 
  
  expect_error(suppressWarnings(estimate_R_agg(incid = weekly_inc, 
                                               dt = 7L, 
                                               dt_out = 7L, 
                                               iter = 10L,
                                               config = config,
                                               method = method,
                                               grid = list(precision = 2.1, min = -1, max = 1))),
               "grid precision must be less than grid max - grid min")
  
  expect_error(suppressWarnings(estimate_R_agg(incid = weekly_inc, 
                                               dt = 7L, 
                                               dt_out = 7L, 
                                               iter = 10L,
                                               config = config,
                                               method = method,
                                               grid = list(precision = "0.001", min = -1, max = 1))),
               "grid precision, min, and max, must all be numeric") 
  
})


test_that("you can't run estimate_R_agg unless incidence is supplied as a vector", {
  method <- "parametric_si"
  config <- make_config(list(mean_si = mean_si,
                             std_si = std_si))
  
inc_obj <- as.incidence(weekly_inc)
expect_error(suppressWarnings(estimate_R_agg(incid = inc_obj, 
                                             dt = 7L, 
                                             dt_out = 7L, 
                                             iter = 10L,
                                             config = config,
                                             method = method,
                                             grid = list(precision = 0.001, min = -1, max = 1))))

inc_df <- data.frame(weekly_inc)
expect_error(suppressWarnings(estimate_R_agg(incid = inc_df, 
                                             dt = 7L, 
                                             dt_out = 7L, 
                                             iter = 10L,
                                             config = config,
                                             method = method,
                                             grid = list(precision = 0.001, min = -1, max = 1))))

})


test_that("you can't run estimate_R_agg unless using parametric/non-parametric SI methods", {
  method <- "uncertain_si"
  config <- make_config(list(mean_si = 2.6, std_mean_si = 1,
                             min_mean_si = 1, max_mean_si = 4.2,
                             std_si = 1.5, std_std_si = 0.5,
                             min_std_si = 0.5, max_std_si = 2.5))
  
  expect_error(estimate_R_agg(incid = weekly_inc, 
                 dt = 7L, 
                 dt_out = 7L, 
                 iter = 10L,
                 config = config,
                 method = method,
                 grid = list(precision = 0.001, min = -1, max = 1)),
               "'arg' should be one of “non_parametric_si”, “parametric_si”")
  
})


# Note: need to fix incorrect incidence plotting
test_that("result of weekly version of estimate_R can be plotted without error", {
  method <- "parametric_si"
  config <- make_config(list(mean_si = mean_si,
                             std_si = std_si))
  res_weekly <- suppressWarnings(estimate_R_agg(incid = weekly_inc, 
                                                dt = 7L, 
                                                dt_out = 7L, 
                                                iter = 10L,
                                                config = config,
                                                method = method,
                                                grid = list(precision = 0.001, min = -1, max = 1)))
  
  expect_error(suppressWarnings(plot(res_weekly),NA))
  
}) 

# Note: need to fix issue with other dts
test_that("method works with different dt", {

  inc_three <- aggregate_inc(incid, 3L)
  inc_ten <- aggregate_inc(incid, 10L)
  
  method <- "parametric_si"
  config <- make_config(list(mean_si = mean_si,
                             std_si = std_si))
  
  res_daily <- suppressWarnings(estimate_R(incid = incid,
                                            config = config,
                                            method = method))
  
  
  estimate_R(incid = incid,
             config = make_config(list(mean_si = mean_si,
                                       std_si = std_si, 
                                       t_start = 30,
                                       t_end = 40)) ,
             method = method)
  
  
  res_three <- suppressWarnings(estimate_R(incid = inc_three,
                                                dt = 3L,
                                                dt_out = 7L,
                                                iter = 10L,
                                                config = config,
                                                method = method))
  
  res_weekly <- suppressWarnings(estimate_R(incid = weekly_inc, 
                                               dt = 7L, 
                                               dt_out = 7L, 
                                               iter = 10L,
                                               config = config,
                                               method = method))
  
  res_ten <- suppressWarnings(estimate_R(incid = inc_ten, 
                                                dt = 10L, 
                                                dt_out = 7L, 
                                                iter = 10L,
                                                config = config,
                                                method = method))
  
  ## test that the ten day incidence matches the aggregated daily one
  ## except for first time window where we impose I = 1 on first day
  
  for(i in seq_len(floor(length(res_daily$I) / 10))[-1])
  {
    expect_equal(sum(res_daily$I[i*10+(1:10)]), 
                 sum(res_ten$I[i*10+(1:10)]))
  }
  
  ## test that weekly sliding R estimates are similar when using different 
  ## temporal aggregations of data
  
  common_t_end <- intersect(res_weekly$R$t_end, res_ten$R$t_end)
  
  relative_error <- abs(res_weekly$R$`Mean(R)`[res_weekly$R$t_end %in% common_t_end] -   
                          res_ten$R$`Mean(R)`[res_ten$R$t_end %in% common_t_end]) / 
    res_weekly$R$`Mean(R)`[res_weekly$R$t_end %in% common_t_end]
  
  ## lot of variation early on so excluding first part
  expect_true(all(relative_error[-c(1:20)] < 0.4)) # 0.4 arbitrarily small
  
})
  