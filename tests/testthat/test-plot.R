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

R_c <- wallinga_teunis(
  Flu2009$incidence,
  method = "non_parametric_si",
  config = list(t_start = seq(2, 26), 
                t_end = seq(8, 32), 
                si_distr = Flu2009$si_distr,
                n_sim = 10L
                )
)

# Create snapshots by system, prevent fragility across developers
system <- Sys.info()[["sysname"]]

test_that("plot.estimate_R doesn't have to include the legend", {
  epi_expect_doppelganger(
    "Flu2009-instantaneous-no-legend", 
    test = "plot",
    plot(R_i, legend = FALSE),
    variant = system
  ) |> suppressWarnings() # TODO: Fix use of aes_string in incidence to avoid this warning

  epi_expect_doppelganger(
    "Flu2009-Rc-instantaneous-no-legend",
    test = "plot",
    plot(R_c, legend = FALSE),
    variant = system
  ) |> suppressWarnings() # TODO: Fix use of aes_string in incidence to avoid this warning
})

test_that("incidence can be plotted separately with imported cases", {
  epi_expect_doppelganger(
    "Flu2009-incidence-import", 
    test = "plot",
    plot(R_i, "incid", add_imported_cases = TRUE),
    variant = system
  ) |> suppressMessages()

  epi_expect_doppelganger(
    "Flu2009-Rc-incidence-import", 
    test = "plot",
    plot(R_c, "incid", add_imported_cases = TRUE),
    variant = system
  ) |> suppressMessages()
})

test_that("serial interval distribution can be plotted separately", {
  epi_expect_doppelganger(
    "Flu2009-SI", 
    test = "plot",
    plot(R_i, "SI"),
    variant = system
  )
  epi_expect_doppelganger(
    "Flu2009-Rc-SI", 
    test = "plot",
    plot(R_c, "SI"),
    variant = system
  )
})

test_that("Reproduction numbers can be plotted separately", {
  epi_expect_doppelganger(
    "Flu2009-Ri", 
    test = "plot",
    plot(R_i, "R", options_R = list(ylim = c(0, 4))),
    variant = system
  )
  epi_expect_doppelganger(
    "Flu2009-Rc",
    test = "plot",
    plot(R_c, "R", options_R = list(ylim = c(0, 4))),
    variant = system
  )
})

