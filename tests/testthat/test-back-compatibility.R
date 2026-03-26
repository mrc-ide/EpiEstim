test_that("These things still calculate the same after refactoring", {
  local_edition(3)
  data("Flu2009")
  
  x1 <- estimate_R(
    Flu2009$incidence,
    method = "non_parametric_si",
    config = make_config(list(si_distr = Flu2009$si_distr))
  )
  expect_snapshot_value(x1, style = "serialize")
  
  x2 <- estimate_R(
    Flu2009$incidence$I,
    method = "non_parametric_si",
    config = make_config(list(si_distr = Flu2009$si_distr))
  )
  expect_snapshot_value(x2, style = "serialize")
  
  x3 <- estimate_R(
    data.frame(Flu2009$incidence$I),
    method = "non_parametric_si",
    config = make_config(list(si_distr = Flu2009$si_distr))
  )
  expect_snapshot_value(x3, style = "serialize")
  
  
})