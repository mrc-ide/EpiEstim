set.seed(5000)
data("Flu2009")
## estimate the instantaneous reproduction number
## (method "non_parametric_si")
R_i <- estimate_R(Flu2009$incidence,
                  method = "non_parametric_si",
                  config = list(t_start = seq(2, 26), 
                                t_end = seq(8, 32), 
                                si_distr = Flu2009$si_distr
                               )
                 )

R_c <- wallinga_teunis(Flu2009$incidence, 
                      method = "non_parametric_si",
                      config = list(t_start = seq(2, 26), 
                                    t_end = seq(8, 32), 
                                    si_distr = Flu2009$si_distr,
                                    n_sim = 10L
                                   )
                     )

test_that("plot.estimate_R doesn't have to include the legend", {
  skip("vdiffr tests not working - see #106")
  skip_if_not_installed("vdiffr")
  vdiffr::expect_doppelganger("Flu2009-instantaneous-no-legend", 
                              plot(R_i, legend = FALSE))

})

test_that("incidence can be plotted separately with imported cases", {
  skip("vdiffr tests not working - see #106")
  skip_if_not_installed("vdiffr")
  vdiffr::expect_doppelganger("Flu2009-incidence-import", 
                              plot(R_i, "incid", add_imported_cases=TRUE))

})

test_that("serial interval distribution can be plotted separately", {
  skip("vdiffr tests not working - see #106")
  skip_if_not_installed("vdiffr")
  vdiffr::expect_doppelganger("Flu2009-SI", 
                              plot(R_i, "SI"))

})

test_that("Reproduction numbers can be plotted separately", {
  skip("vdiffr tests not working - see #106")
  skip_if_not_installed("vdiffr")
  vdiffr::expect_doppelganger("Flu2009-Ri", 
                              plot(R_i, "R", options_R = list(ylim = c(0, 4))))
  vdiffr::expect_doppelganger("Flu2009-Rc",
                              plot(R_c, "R", options_R = list(ylim = c(0, 4))))

})
