# Shared setup ----------------------------------------------------------------
n_v <- 2   # 2 variants
n_loc <- 3 # 3 locations
time_stps <- 100   # 100 time steps
priors <- default_priors()

w_v <- c(0, 0.2, 0.5, 0.3)
si_distr <- cbind(w_v, w_v)
epsilon <- 1

incid0 <- array(10, dim = c(time_stps, n_loc, n_v))
incid_mv0 <- process_I_multivariant(incid0)
lambda0 <- compute_lambda(incid_mv0, si_distr)

#R <- matrix(1, nrow = T, ncol = n_loc)
#R[1, ] <- NA

test_that("defaults", {
  expect_silent(default_priors()) |>
    expect_equal(
      list(epsilon = list(shape = 1, scale = 1), R = list(shape = 0.04, scale = 25))    
    )
  expect_silent(default_mcmc_controls()) |>
    expect_equal(
      list(n_iter = 1100, burnin = 10, thin = 10)
    )
})


test_that("process_I_multivariant()", {
  incid <- array(10, dim = c(100, 1, 2))

  epi_snapshot_value(process_I_multivariant(incid))
})

test_that("get_shape_R_flat()", {
  incid <- array(10, dim = c(100, 3, 2))

  expect_silent(get_shape_R_flat(incid, default_priors())) |>
    expect_all_equal(20.04)
})

test_that("compute_lambda()", {
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v)
  incid <- array(10, dim = c(100, 3, 2))
  incid_mv <- process_I_multivariant(incid)

  expect_silent(lambda <- compute_lambda(incid_mv, si_distr))
  epi_snapshot_value(lambda)
})

test_that("Epsilon: get_shape_epsilon(), draw_epsilon()", {
  incid <- array(10, dim = c(100, 4, 3))
  incid_mv <- process_I_multivariant(incid)
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v, w_v)
  lambda <- compute_lambda(incid_mv, si_distr)

  expect_silent(get_shape_epsilon(incid_mv$local, lambda, default_priors())) |>
    expect_equal(c(3961, 3961))
  
  R <- matrix(1, nrow = 100, ncol = 4)
  R[1, ] <- NA # no estimates of R on first time step

  expect_silent(
    draw_epsilon(R, incid_mv$local, lambda, default_priors(), seed = 1)
  ) |>
    expect_equal(c(1.001066, 1.032584), tolerance = 0.000001)
})

test_that("draw_R()", {
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v)
  incid <- array(10, dim = c(100, 3, 2))
  incid_mv <- process_I_multivariant(incid)
  lambda <- compute_lambda(incid_mv, si_distr)

  expect_silent(result <- draw_R(1, incid_mv$local, lambda, default_priors(),
                                 seed = 1, t_min = 2L))
  epi_snapshot_value(result)
})

test_that("first_nonzero_incid()", {
  incid <- array(10, dim = c(100, 3, 2))
  incid_mv <- process_I_multivariant(incid)

  expect_silent(first_nonzero_incid(incid_mv$local)) |>
    expect_equal(2)
})

test_that("compute_si_cutoff snapshot", {
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v)
  
  expect_silent(compute_si_cutoff(si_distr)) |>
    expect_equal(4)
})

test_that("compute_t_min snapshot", {
  incid <- array(10, dim = c(100, 3, 2))
  incid_mv <- process_I_multivariant(incid)

  expect_silent(compute_t_min(incid_mv$local, si_distr)) |>
    expect_equal(6)
})

test_that("estimate_advantage()", {
  incid <- array(10, dim = c(100, 3, 2))
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v)
  
  expect_message(
    est <- estimate_advantage(incid, si_distr, default_priors(), seed = 1)
  )

  expect_s3_class(est$diag[[1]], "gelman.diag")
  epi_snapshot_value(est)
})

