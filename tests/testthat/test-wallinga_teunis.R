data(Flu2009)

test_that(
  "wallinga_teunis() outputs the right errors/warnings",
  {

    msg <- "method non_parametric_si requires to specify the config\\$mean_si argument."
    expect_error(
      wallinga_teunis(
        1:10,
        method = "parametric_si",
        config = list(t_start = 3,
                      t_end = 6,
                      n_sim = 10)
      ),
      msg
    )

    msg <- "method non_parametric_si requires to specify the config\\$std_si argument."
    expect_error(
      wallinga_teunis(
        1:10,
        method = "parametric_si",
        config = list(t_start = 3,
                      t_end = 6,
                      n_sim = 10,
                      mean_si = 3)
      ),
      msg
    )
    
    msg <- "method parametric_si requires a value >1 for config\\$mean_si."
    expect_error(
      wallinga_teunis(
        1:10,
        method = "parametric_si",
        config = list(t_start = 3,
                      t_end = 6,
                      n_sim = 10,
                      mean_si = -2,
                      std_si = 3)
      ),
      msg
    )

  }
)




test_that("wallinga_teunis(): results have the right shape/format", {

  i <- flat_incid <- rep(10, 100)
  res <- wallinga_teunis(
    i,
    method = "non_parametric_si",
    config = list(t_start = 2:100, t_end = 2:100,
                  si_distr = Flu2009$si_distr, 
                  seed = 1, 
                  n_sim = 50)
  )
  
  expect_equal(class(res), "wallinga_teunis")
  expect_true(is.list(res))
  exp_names <- c("R", "method", "si_distr", "SI.Moments", "dates", "I", "I_local", 
                 "I_imported")
  expect_identical(names(res), exp_names)
  expect_identical(dim(res$R), c(99L, 6L))
  
}
)





test_that("wallinga_teunis() gives expected results with flat incidence", {
  i <- flat_incid <- rep(10, 100)
  res <- wallinga_teunis(
    i,
    method = "non_parametric_si",
    config = list(t_start = 2:100, t_end = 2:100,
                  si_distr = Flu2009$si_distr, 
                  seed = 1, 
                  n_sim = 50)
  )
  
  ### expect final R estimate to be zero
  expect_equal(tail(res$R$`Mean(R)`, 1), 0)
  
  ### expect R estimates to be ~1 except very early and very late ones
  expect_equal(mean(res$R$`Mean(R)`[10:80]), 1)

})





test_that("wallinga_teunis() gives expected results with increasing incidence", {

  ## Expected results
  ##
  ## The first epicurve grows exponentially, we should find R > 1
  ## The second epicurve grows faster, we should find R_2 > R_1
  i_1 <- rpois(100, lambda = exp(0.0523 * 1:100))
  i_2 <- rpois(100, lambda = exp(0.0765 * 1:100))
  
  res_1 <- wallinga_teunis(
    i_1,
    method = "non_parametric_si",
    config = list(t_start = 10, t_end = 90,
                  si_distr = Flu2009$si_distr, 
                  seed = 1, 
                  n_sim = 50)
  )
  res_2 <- wallinga_teunis(
    i_2,
    method = "non_parametric_si",
    config = list(t_start = 10, t_end = 90,
                  si_distr = Flu2009$si_distr, 
                  seed = 1, 
                  n_sim = 50)
  )
  
  expect_true(res_1$R$`Quantile.0.025(R)` > 1)
  expect_true(res_2$R$`Quantile.0.025(R)` > res_1$R$`Quantile.0.975(R)`)

})



test_that(
  "seed fixing in wallinga_teunis config works as expected", {
    
    # seed = 1 - acts as reference for this test
    out1 <- wallinga_teunis(Flu2009$incidence, method = "non_parametric_si",
                                                config = list(t_start = 2:26, t_end = 8:32,
                                                              si_distr = Flu2009$si_distr, 
                                                              n_sim = 50,
                                                              seed = 1))
    
    # same seed
    out2<- wallinga_teunis(Flu2009$incidence, method = "non_parametric_si",
                            config = list(t_start = 2:26, t_end = 8:32,
                                          si_distr = Flu2009$si_distr, 
                                          n_sim = 50,
                                          seed = 1))
    
    # different seed
    out3 <- wallinga_teunis(Flu2009$incidence, method = "non_parametric_si",
                            config = list(t_start = 2:26, t_end = 8:32,
                                          si_distr = Flu2009$si_distr, 
                                          n_sim = 50,
                                          seed = 2))
    # no seed
    out4 <- wallinga_teunis(Flu2009$incidence, method = "non_parametric_si",
                            config = list(t_start = 2:26, t_end = 8:32,
                                          si_distr = Flu2009$si_distr, 
                                          n_sim = 50,
                                          seed = 2))
    
    # same seed should produce same R estimates
    expect_equal(out1$R, out2$R)
    
    # different seed should produce same mean R estimates but different 
    # uncertainty
    expect_equal(out1$R$`Mean(R)`, out3$R$`Mean(R)`)
    expect_true(any(out1$R$`Quantile.0.025(R)` != out3$R$`Quantile.0.025(R)`))
    expect_true(any(out1$R$`Quantile.0.975(R)` != out3$R$`Quantile.0.975(R)`))
    
    # not fixing the seed should produce same mean R estimates but different 
    # uncertainty
    expect_equal(out1$R$`Mean(R)`, out4$R$`Mean(R)`)
    expect_true(any(out1$R$`Quantile.0.025(R)` != out4$R$`Quantile.0.025(R)`))
    expect_true(any(out1$R$`Quantile.0.975(R)` != out4$R$`Quantile.0.975(R)`))
     }
)

## test_that("Example 1 matches saved output", {
##   out <- wallinga_teunis(Flu2009$incidence, method = "non_parametric_si",
##                     config = list(t_start = 2:26, t_end = 8:32,
##                                   si_distr = Flu2009$si_distr, 
##                                   seed = 1, 
##                                   n_sim = 50))
##   ## TODO: save output (with fixed seed) into am rda as we did for estimate_T
##   #expect_equal_to_reference(out, "../expected_output/Example1.rda", update = FALSE)
## })
