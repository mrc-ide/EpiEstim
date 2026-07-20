context("Gibbs samplers")

test_that("draw_epsilon produces expected results (2 variants, 4 locations)", {
  n_v <- 2 # 2 variants
  n_loc <- 4 # 4 locations
  T <- 100 # 100 time steps

  priors <- default_priors()

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))
  incid <- process_I_multivariant(incid)

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v)
  lambda <- compute_lambda(incid, si_distr)

  # Constant reproduction number of 1
  R <- matrix(1, nrow = T, ncol = n_loc)
  R[1, ] <- NA # no estimates of R on first time step

  set.seed(1)
  x <- sapply(1:1000, function(e) draw_epsilon(R, incid$local, lambda, priors))

  ## epsilon should be approximately 1
  ## not exactly 1 because of the first few timesteps & because of priors
  expect_equal(mean(x), 1, tolerance = 0.05)
})


test_that("draw_epsilon produces expected results (2 variants, 1 location)", {
  n_v <- 2 # 2 variants
  n_loc <- 1 # 1 location
  T <- 100 # 100 time steps

  priors <- default_priors()

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))
  incid <- process_I_multivariant(incid)

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v)
  lambda <- compute_lambda(incid, si_distr)

  # Constant reproduction number of 1
  R <- matrix(1, nrow = T, ncol = n_loc)
  R[1, ] <- NA # no estimates of R on first time step

  set.seed(1)
  x <- sapply(1:1000, function(e) draw_epsilon(R, incid$local, lambda, priors))

  ## epsilon should be approximately 1
  ## not exactly 1 because of the first few timesteps & because of priors
  expect_equal(mean(x), 1, tolerance = 0.05)
})


test_that("draw_epsilon produces expected results (>2 variants, 4 locations)", {
  n_v <- 3 # 3 variants
  n_loc <- 4 # 4 locations
  T <- 100 # 100 time steps

  priors <- default_priors()

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))
  incid <- process_I_multivariant(incid)

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v, w_v)
  lambda <- compute_lambda(incid, si_distr)

  # Constant reproduction number of 1
  R <- matrix(1, nrow = T, ncol = n_loc)
  R[1, ] <- NA # no estimates of R on first time step

  set.seed(1)
  x <- sapply(1:1000, function(e) draw_epsilon(R, incid$local, lambda, priors))

  ## epsilon should be approximately 1
  ## not exactly 1 because of the first few timesteps & because of priors
  expect_equal(
    rowMeans(x), c(1, 1), tolerance = 0.05
  )
})


test_that("draw_epsilon produces expected results (>2 variants, 1 location)", {
  n_v <- 3 # 3 variants
  n_loc <- 1 # 1 location
  T <- 100 # 100 time steps

  priors <- default_priors()

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))
  incid <- process_I_multivariant(incid)

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v, w_v)
  lambda <- compute_lambda(incid, si_distr)

  # Constant reproduction number of 1
  R <- matrix(1, nrow = T, ncol = n_loc)
  R[1, ] <- NA # no estimates of R on first time step

  set.seed(1)
  x <- sapply(1:1000, function(e) draw_epsilon(R, incid$local, lambda, priors))

  ## epsilon should be approximately 1
  ## not exactly 1 because of the first few timesteps & because of priors
  expect_equal(rowMeans(x), c(1, 1), tolerance = 0.05)
})


test_that("draw_R produces expected results (2 variants, 4 locations)", {
  n_v <- 2 # 2 variants
  n_loc <- 4 # 4 locations
  T <- 100 # 100 time steps

  priors <- default_priors()

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))
  incid <- process_I_multivariant(incid)

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v)
  lambda <- compute_lambda(incid, si_distr)

  # Epsilon = 1 i.e. no transmission advantage
  epsilon <- 1

  set.seed(1)
  x <- lapply(1:1000, function(e) draw_R(epsilon, incid$local, lambda, priors, t_min = 2L))
  x_mean <- Reduce("+", x) / length(x)

  ## R should be approximately 1
  ## not exactly 1 because of the first few timesteps & because of priors
  ## so ignore fisrt timesteps
  expect_lt(max(abs(x_mean[-c(1, 2, 3), ] - 1)), 0.05)
})


test_that("draw_R produces expected results (2 variants, 1 location)", {
  n_v <- 2 # 2 variants
  n_loc <- 1 # 1 locations
  T <- 100 # 100 time steps

  priors <- default_priors()

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))
  incid <- process_I_multivariant(incid)

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v)
  lambda <- compute_lambda(incid, si_distr)

  # Epsilon = 1 i.e. no transmission advantage
  epsilon <- 1

  set.seed(1)
  x <- lapply(1:1000, function(e) draw_R(epsilon, incid$local, lambda, priors, t_min = 2L))
  x_mean <- Reduce("+", x) / length(x)

  ## R should be approximately 1
  ## not exactly 1 because of the first few timesteps & because of priors
  ## so ignore fisrt timesteps
  expect_lt(max(abs(x_mean[-c(1, 2, 3), ] - 1)), 0.05)
})


test_that("draw_R produces expected results (>2 variants, 4 locations)", {
  n_v <- 3 # 3 variants
  n_loc <- 4 # 4 locations
  T <- 100 # 100 time steps

  priors <- default_priors()

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))
  incid <- process_I_multivariant(incid)

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v, w_v)
  lambda <- compute_lambda(incid, si_distr)

  # Epsilon = 1 i.e. no transmission advantage
  epsilon <- c(1, 1)

  set.seed(1)
  x <- lapply(1:1000, function(e) draw_R(epsilon, incid$local, lambda, priors, t_min = 2L))
  x_mean <- Reduce("+", x) / length(x)

  ## R should be approximately 1
  ## not exactly 1 because of the first few timesteps & because of priors
  ## so ignore fisrt timesteps
  expect_lt(max(abs(x_mean[-c(1, 2, 3), ] - 1)), 0.05)
})


test_that("draw_R produces expected results (>2 variants, 1 location)", {
  n_v <- 3 # 3 variants
  n_loc <- 1 # 1 locations
  T <- 100 # 100 time steps

  priors <- default_priors()

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))
  incid <- process_I_multivariant(incid)

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v, w_v)
  lambda <- compute_lambda(incid, si_distr)

  # Epsilon = 1 i.e. no transmission advantage
  epsilon <- c(1, 1)

  set.seed(1)
  x <- lapply(1:1000, function(e) draw_R(epsilon, incid$local, lambda, priors, t_min = 2L))
  x_mean <- Reduce("+", x) / length(x)

  ## R should be approximately 1
  ## not exactly 1 because of the first few timesteps & because of priors
  ## so ignore fisrt timesteps
  expect_lt(max(abs(x_mean[-c(1, 2, 3), ] - 1)), 0.05)
})


test_that("estimate_advantage produces expected results (2 variants 3 locations)", {
  n_v <- 2 # 2 variants
  n_loc <- 3 # 3 locations
  T <- 100 # 100 time steps

  priors <- default_priors()

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v)

  x <- estimate_advantage(incid, si_distr, priors, seed = 1, t_min = 2L)

  ## epsilon should be approximately 1
  expect_equal(mean(x$epsilon), 1, tolerance = 0.05)

  ## R should be approximately 1
  ## not exactly 1 because of the first few timesteps & because of priors
  ## so ignore fisrt timesteps
  mean_R <- rowMeans(x$R, dims = 2)
  expect_lt(max(abs(mean_R[-c(1, 2, 3), ] - 1)), 0.1)
})


test_that("estimate_advantage produces expected results (2 variants 1 location)", {
  n_v <- 2 # 2 variants
  n_loc <- 1 # 1 locations
  T <- 100 # 100 time steps

  priors <- default_priors()

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v)

  x <- estimate_advantage(incid, si_distr, priors, seed = 1, t_min = 2L)

  ## epsilon should be approximately 1
  expect_equal(mean(x$epsilon), 1, tolerance = 0.05)

  ## R should be approximately 1
  ## not exactly 1 because of the first few timesteps & because of priors
  ## so ignore fisrt timesteps
  mean_R <- rowMeans(x$R, dims = 2)
  expect_lt(max(abs(mean_R[-c(1, 2, 3), ] - 1)), 0.1)
})


test_that("estimate_advantage produces expected results (>2 variants 4 locs)", {
  n_v <- 3 # 3 variants
  n_loc <- 4 # 4 locations
  T <- 100 # 100 time steps

  priors <- default_priors()

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v, w_v)

  x <- estimate_advantage(incid, si_distr, priors, seed = 1, t_min = 2L)

  ## epsilon should be approximately 1
  ## FIXME this should be apply(x$epsilon, 1, mean)
  ## as 2 epsilons are returned here.

  expect_equal(
    rowMeans(x$epsilon), c(1, 1), tolerance = 0.05
  )

  ## R should be approximately 1
  ## not exactly 1 because of the first few timesteps & because of priors
  ## so ignore fisrt timesteps
  mean_R <- rowMeans(x$R, dims = 2)
  expect_lt(max(abs(mean_R[-c(1, 2, 3), ] - 1)), 0.1)
})


test_that("estimate_advantage produces expected results (>2 variants 1 loc)", {
  n_v <- 3 # 3 variants
  n_loc <- 1 # 1 location
  T <- 100 # 100 time steps

  priors <- default_priors()

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v, w_v)

  x <- estimate_advantage(incid, si_distr, priors, seed = 1, t_min = 2L)

  ## epsilon should be approximately 1
  expect_equal(
    rowMeans(x$epsilon), c(1, 1), tolerance = 0.05
  )

  ## R should be approximately 1
  ## not exactly 1 because of the first few timesteps & because of priors
  ## so ignore fisrt timesteps
  mean_R <- rowMeans(x$R, dims = 2)
  expect_lt(max(abs(mean_R[-c(1, 2, 3), ] - 1)), 0.1)
})




test_that("process_I_multivariant rejects wrong inputs", {
  n_v <- 3 # 3 variants
  n_loc <- 1 # 1 location
  T <- 100 # 100 time steps

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))
  incid_imported <- array(1, dim = c(T, n_loc, n_v))

  expect_error(process_I_multivariant(incid, incid_imported[-1, , ]),
               "'incid' and 'incid_imported' have incompatible dimensions")
})


test_that("process_I_multivariant works as expected", {
  n_v <- 3 # 3 variants
  n_loc <- 1 # 1 location
  T <- 100 # 100 time steps

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))
  incid_imported <- array(1, dim = c(T, n_loc, n_v))

  ## specifying imported cases
  incid_processed <- process_I_multivariant(incid, incid_imported)
  expect_equal(incid_processed$local + incid_processed$imported, incid)
  expect_equal(incid_processed$imported, incid_imported)

  ## with the default
  incid_processed <- process_I_multivariant(incid)
  expect_equal(incid_processed$local + incid_processed$imported, incid)
  expect_true(all(incid_processed$imported[-1, , ] == 0))
  expect_true(all(incid_processed$imported[1, , ] == incid[1, , ]))
  expect_true(all(incid_processed$local[-1, , ] == incid[-1, , ]))
  expect_true(all(incid_processed$local[1, , ] == 0))
})


test_that("compute_lambda rejects invalid incid inputs", {
  n_v <- 3 # 3 variants
  n_loc <- 1 # 1 location
  T <- 100 # 100 time steps

  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v)

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))

  expect_error(compute_lambda(incid, si_distr),
      "'incid 'should be an 'incid_multivariant' object.")
})


test_that("estimate_advantage produces expected results (>2var, 1loc, imports)", {
  n_v <- 3 # 3 variants
  n_loc <- 1 # 1 location
  T <- 100 # 100 time steps

  priors <- default_priors()

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))
  incid_imported <- array(0, dim = c(T, n_loc, n_v))
  # make all the cases of non reference variant imported
  incid_imported[, , 2:n_v] <- incid[, , 2:n_v]
  # all cases at first time step must be imported
  incid_imported[1, , ] <- incid[1, , ]

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v, w_v)

  x <- estimate_advantage(incid, si_distr, priors, seed = 1,
                      incid_imported = incid_imported, t_min = 2L)

  ## epsilon should be approximately 0
  expect_equal(
    rowMeans(x$epsilon), c(0, 0), tolerance = 0.05
  )

  ## R should be approximately 1
  ## not exactly 1 because of the first few timesteps & because of priors
  ## so ignore fisrt timesteps
  mean_R <- rowMeans(x$R, dims = 2)
  expect_lt(max(abs(mean_R[-c(1, 2, 3), ] - 1)), 0.1)
})


test_that("estimate_advantage produces expected results (>2var, 4loc, imports)", {
  n_v <- 3 # 3 variants
  n_loc <- 4 # 4 locations
  T <- 100 # 100 time steps

  priors <- default_priors()

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))
  incid_imported <- array(0, dim = c(T, n_loc, n_v))
  # make all the cases of non reference variant imported
  incid_imported[, , 2:n_v] <- incid[, , 2:n_v]
  # all cases at first time step must be imported
  incid_imported[1, , ] <- incid[1, , ]

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v, w_v)
  ## Need to run for longer after changing priors
  x <- estimate_advantage(
    incid, si_distr, priors, seed = 1,
    incid_imported = incid_imported,
    mcmc_control = list(n_iter = 2000L, burnin = 100L, thin = 10L), t_min = 2L
  )

  ## epsilon should be approximately 0
  expect_equal(
    rowMeans(x$epsilon), c(0, 0), tolerance = 0.05
  )

  ## R should be approximately 1
  ## not exactly 1 because of the first few timesteps & because of priors
  ## so ignore fisrt timesteps
  mean_R <- rowMeans(x$R, dims = 2)
  expect_lt(max(abs(mean_R[-c(1, 2, 3), ] - 1)), 0.1)
})


test_that("estimate_advantage produces expected results (2var, 1loc, imports)", {
  n_v <- 2 # 2 variants
  n_loc <- 1 # 1 location
  T <- 100 # 100 time steps

  priors <- default_priors()

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))
  incid_imported <- array(0, dim = c(T, n_loc, n_v))
  # make all the cases of second variant and third imported
  incid_imported[, , 2] <- incid[, , 2]
  # all cases at first time step must be imported
  incid_imported[1, , ] <- incid[1, , ]

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v, w_v)

  x <- estimate_advantage(incid, si_distr, priors, seed = 1,
                      incid_imported = incid_imported, t_min = 2L)

  ## epsilon should be approximately 0
  expect_equal(mean(x$epsilon), 0, tolerance = 0.05)

  ## R should be approximately 1
  ## not exactly 1 because of the first few timesteps & because of priors
  ## so ignore fisrt timesteps
  mean_R <- rowMeans(x$R, dims = 2)
  expect_lt(max(abs(mean_R[-c(1, 2, 3), ] - 1)), 0.1)
})


test_that("estimate_advantage produces expected results (2var, 4loc, imports)", {
  n_v <- 2 # 2 variants
  n_loc <- 4 # 4 locations
  T <- 100 # 100 time steps

  priors <- default_priors()

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))
  incid_imported <- array(0, dim = c(T, n_loc, n_v))
  # make all the cases of second variant and third imported
  incid_imported[, , 2] <- incid[, , 2]
  # all cases at first time step must be imported
  incid_imported[1, , ] <- incid[1, , ]

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v)

  x <- estimate_advantage(
    incid, si_distr, priors, seed = 1,
    incid_imported = incid_imported,
    mcmc_control = list(n_iter = 2000L, burnin = 100L, thin = 10L), t_min = 2L
  )

  ## epsilon should be approximately 0
  expect_equal(mean(x$epsilon), 0, tolerance = 0.05)

  ## R should be approximately 1
  ## not exactly 1 because of the first few timesteps & because of priors
  ## so ignore fisrt timesteps
  mean_R <- rowMeans(x$R, dims = 2)
  expect_lt(max(abs(mean_R[-c(1, 2, 3), ] - 1)), 0.1)
})


test_that("estimate_advantage produces expected results (2 var, 2 loc, R_loc1 = 1.1, R_loc2 = 1.5)", {
  skip_if_not_installed("projections")
  n_v <- 2 # 2 variants
  n_loc <- 2 # 2 locations
  T <- 100 # 100 time steps

  transm_adv <- 1.5 # Var 2 has TA of 1.5
  R_loc1 <- 1.1
  R_loc2 <- 1.5

  R_L1V1 <- R_loc1
  R_L1V2 <- R_loc1*transm_adv
  R_L2V1 <- R_loc2
  R_L2V2 <- R_loc2*transm_adv

  R <- array(NA, dim = c(T, n_loc, n_v))
  R[,1,1] <- rep(R_L1V1, each=T)
  R[,2,1] <- rep(R_L2V1, each=T)
  R[,1,2] <- rep(R_L1V2, each=T)
  R[,2,2] <- rep(R_L2V2, each=T)

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v)

  # simulate incidence
  incid_init <- incidence::incidence(rep(1, 20))
  incid <- array(NA, dim = c(T, n_loc, n_v))

  for (loc in seq_len(n_loc)) {
    for (v in seq_len(n_v)) {
      incid[, loc, v] <- rbind(
        incid_init$counts,
        as.matrix( #
          projections::project(
            incid_init,
            ## R in the future so removing time of seeding
            R = R[-1, loc, v],
            si = si_distr[-1, v],
            n_sim = 1,
            n_days = T - 1,
            time_change = seq_len(
              length(R[, loc, v]) - 2
            ) - 1
          )
        )
      )
    }
  }


  priors <- default_priors()
  x <- estimate_advantage(
    incid, si_distr, priors, seed = 1, t_min = 2L
  )

  ## R should be approx 1.1 for loc1 and 1.5 for loc2
  expect_equal(mean(x$R[,1,], na.rm = TRUE), 1.1, tolerance = 0.5)
  expect_equal(mean(x$R[,2,], na.rm = TRUE), 1.5, tolerance = 0.5)

})

test_that("estimate_advantage faster with precompute (2 variants 3 locations)", {
  n_v <- 2 # 2 variants
  n_loc <- 3 # 3 locations
  T <- 100 # 100 time steps

  priors <- default_priors()

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v)

  t1 <- system.time(
    x1 <- estimate_advantage(incid, si_distr, priors, seed = 1, precompute = TRUE, t_min = 2L)
  )

  t2 <- system.time(
    x2 <- estimate_advantage(incid, si_distr, priors, seed = 1, precompute = FALSE, t_min = 2L)
  )

  ## t1 should be < t2
  expect_lt(t1[["elapsed"]], t2[["elapsed"]])

  ## epsilon should be approximately 1 in both cases
  expect_equal(mean(x1$epsilon), 1, tolerance = 0.05)
  expect_equal(mean(x2$epsilon), 1, tolerance = 0.05)

  ## R should be approximately 1 in both cases
  ## not exactly 1 because of the first few timesteps & because of priors
  ## so ignore fisrt timesteps
  mean_R1 <- rowMeans(x1$R, dims = 2)
  expect_lt(max(abs(mean_R1[-c(1, 2, 3), ] - 1)), 0.1)
  mean_R2 <- rowMeans(x2$R, dims = 2)
  expect_lt(max(abs(mean_R2[-c(1, 2, 3), ] - 1)), 0.1)
})


test_that("estimate_advantage faster with precompute (2 variants 1 location)", {
  n_v <- 2 # 2 variants
  n_loc <- 1 # 1 locations
  T <- 100 # 100 time steps

  priors <- default_priors()

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v)

  t1 <- system.time(
    x1 <- estimate_advantage(incid, si_distr, priors, seed = 1, precompute = TRUE, t_min = 2L)
  )

  t2 <- system.time(
    x2 <- estimate_advantage(incid, si_distr, priors, seed = 1, precompute = FALSE, t_min = 2L)
  )

  ## t1 should be < t2
  expect_lt(t1[["elapsed"]], t2[["elapsed"]])

  ## epsilon should be approximately 1 in both cases
  expect_equal(mean(x1$epsilon), 1, tolerance = 0.05)
  expect_equal(mean(x2$epsilon), 1, tolerance = 0.05)

  ## R should be approximately 1 in both cases
  ## not exactly 1 because of the first few timesteps & because of priors
  ## so ignore fisrt timesteps
  mean_R1 <- rowMeans(x1$R, dims = 2)
  expect_lt(max(abs(mean_R1[-c(1, 2, 3), ] - 1)), 0.1)
  mean_R2 <- rowMeans(x2$R, dims = 2)
  expect_lt(max(abs(mean_R2[-c(1, 2, 3), ] - 1)), 0.1)
})


test_that("estimate_advantage faster with precompute (3 variants 4 locations)", {
  n_v <- 3 # 3 variants
  n_loc <- 4 # 4 locations
  T <- 100 # 100 time steps

  priors <- default_priors()

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v, w_v)

  t1 <- system.time(
    x1 <- estimate_advantage(incid, si_distr, priors, seed = 1, precompute = TRUE, t_min = 2L)
  )

  t2 <- system.time(
    x2 <- estimate_advantage(incid, si_distr, priors, seed = 1, precompute = FALSE, t_min = 2L)
  )

  ## t1 should be < t2
  expect_lt(t1[["elapsed"]], t2[["elapsed"]])

    ## epsilon should be approximately 1 in both cases
  expect_equal(mean(x1$epsilon), 1, tolerance = 0.05)
  expect_equal(mean(x2$epsilon), 1, tolerance = 0.05)

  ## R should be approximately 1 in both cases
  ## not exactly 1 because of the first few timesteps & because of priors
  ## so ignore fisrt timesteps
  mean_R1 <- rowMeans(x1$R, dims = 2)
  expect_lt(max(abs(mean_R1[-c(1, 2, 3), ] - 1)), 0.1)
  mean_R2 <- rowMeans(x2$R, dims = 2)
  expect_lt(max(abs(mean_R2[-c(1, 2, 3), ] - 1)), 0.1)
})


test_that("estimate_advantage faster with precompute (3 variants 1 location)", {
  n_v <- 3 # 3 variants
  n_loc <- 1 # 1 locations
  T <- 100 # 100 time steps

  priors <- default_priors()

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v, w_v)
  
  mcmc_control <- default_mcmc_controls()
  mcmc_control$n_iter <- 2100L

  t1 <- system.time(
    x1 <- estimate_advantage(incid, si_distr, priors, seed = 1, 
                             precompute = TRUE, t_min = 2L, 
                             mcmc_control = mcmc_control)
  )

  t2 <- system.time(
    x2 <- estimate_advantage(incid, si_distr, priors, seed = 1, 
                             precompute = FALSE, t_min = 2L, 
                             mcmc_control = mcmc_control)
  )

  ## t1 should be < t2
  expect_lt(t1[["elapsed"]], t2[["elapsed"]])

  ## epsilon should be approximately 1 in both cases
  expect_equal(mean(x1$epsilon), 1, tolerance = 0.05)
  expect_equal(mean(x2$epsilon), 1, tolerance = 0.05)

  ## R should be approximately 1 in both cases
  ## not exactly 1 because of the first few timesteps & because of priors
  ## so ignore fisrt timesteps
  mean_R1 <- rowMeans(x1$R, dims = 2)
  expect_lt(max(abs(mean_R1[-c(1, 2, 3), ] - 1)), 0.1)
  mean_R2 <- rowMeans(x2$R, dims = 2)
  expect_lt(max(abs(mean_R2[-c(1, 2, 3), ] - 1)), 0.1)
})


test_that("compute_si_cutoff returns the correct value", {
  si_1 <- c(0.1, 0.2, 0.3, 0.2, 0.1, 0.03, 0.02,
            0.01, 0.01, 0.01, 0.02)
  si <- cbind(si_1, si_1)
  ## With default cut-off we expect it to be 9
  expect_equal(compute_si_cutoff(si), 7L)
  ## Change miss_at_most to 0.5. We now should get
  ## 3
  expect_equal(compute_si_cutoff(si, 0.5), 3L)
  ## With different SI distributions for the 2
  ## variants
  si_2 <- c(1, rep(0, 10))
  si <- cbind(si_1, si_2)
  expect_equal(compute_si_cutoff(si), 7L)
  expect_equal(compute_si_cutoff(si, 0.5), 3L)
})

test_that("first_nonzero_incid returns correct value", {
  incid <- array(0, dim = c(10, 2, 2))
  ## Put the non-zero incidence in different
  ## places so that we can be sure we are getting
  ## the right index back.
  incid[2, 1, 1] <- 1
  incid[3, 2, 1] <- 1
  incid[4, 1, 2] <- 1
  incid[5, 2, 2] <- 1
  expect_equal(first_nonzero_incid(incid), 5L)

  ## Edge cases
  incid <- array(0, dim = c(10, 2, 2))
  incid[1, , ] <- 1
  expect_equal(first_nonzero_incid(incid), 1L)

  incid <- array(0, dim = c(10, 2, 2))
  incid[10, , ] <- 1
  expect_equal(first_nonzero_incid(incid), 10L)
})

test_that("compute_t_min works correctly", {
  incid <- array(0, dim = c(10, 2, 2))
  incid[2, 1, 1] <- 1
  incid[3, 2, 1] <- 1
  incid[4, 1, 2] <- 1
  incid[5, 2, 2] <- 1
  si_1 <- c(0.1, 0.2, 0.3, 0.2, 0.1, 0.03, 0.02,
            0.01, 0.01, 0.01, 0.02)
  si <- cbind(si_1, si_1)
  ## first_nonzero_incid will return 5 and
  ## compute_si_cutoff  will be 7 so that here we
  ## expect 12.
  expect_equal(compute_t_min(incid, si), 12L)
})

test_that("estimate_advantage uses the correct t_min", {
  n_v <- 2 # 2 variants
  n_loc <- 3 # 3 locations
  T <- 100 # 100 time steps

  priors <- default_priors()

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v)

  x <- estimate_advantage(incid, si_distr, priors, t_min = 2L, seed = 1)
  ## If t_min is 2, the first row if x$R will be NA
  expect_true(all(is.na(x$R[1, , ])))
  ## and not after that.
  expect_false(anyNA(x$R[seq(2, dim(x$R)[1]), , ]))
  ## if t_min is NULL, t_min would be set to
  ## compute_t_min.
  t_min <- compute_t_min(incid, si_distr)
  x <- estimate_advantage(incid, si_distr, priors, seed = 1)
  expect_true(all(is.na(x$R[seq(1, t_min - 1, 1), , ])))
  expect_false(anyNA(x$R[seq(t_min, dim(x$R)[1]), , ]))
})

test_that("estimate_advantage parametric_si matches an equivalent discrete SI matrix", {
  n_v <- 2 # 2 variants
  n_loc <- 2 # 2 locations
  T <- 40 # 40 time steps

  priors <- default_priors()
  mcmc_control <- list(n_iter = 100L, burnin = 20L, thin = 5L)

  incid <- array(10, dim = c(T, n_loc, n_v))
  mean_si <- c(2.6, 3.1)
  std_si <- c(1.5, 1.2)

  si_parametric <- matrix(NA_real_, nrow = T, ncol = n_v)
  for (k in seq_len(n_v)) {
    col <- discr_si(seq(0, T - 1), mean_si[k], std_si[k])
    si_parametric[, k] <- col / sum(col)
  }

  parametric_out <- estimate_advantage(
    incid,
    priors = priors,
    mcmc_control = mcmc_control,
    seed = 1,
    t_min = 2L,
    method = "parametric_si",
    mean_si = mean_si,
    std_si = std_si
  )

  discrete_out <- estimate_advantage(
    incid,
    si_parametric,
    priors = priors,
    mcmc_control = mcmc_control,
    seed = 1,
    t_min = 2L
  )

  expect_equal(parametric_out$epsilon, discrete_out$epsilon)
  expect_equal(parametric_out$R, discrete_out$R)
  expect_equal(parametric_out$convergence, discrete_out$convergence)
  expect_equal(parametric_out$diag, discrete_out$diag)
  expect_equal(parametric_out$si_distr, discrete_out$si_distr)
  expect_equal(parametric_out$SI.Moments, discrete_out$SI.Moments)
})

test_that("estimate_advantage parametric_si validates SI inputs", {
  incid <- array(10, dim = c(20, 2, 2))
  priors <- default_priors()

  expect_error(
    estimate_advantage(incid, priors = priors, method = "parametric_si", std_si = 1.5),
    "mean_si must be supplied"
  )
  expect_error(
    estimate_advantage(incid, priors = priors, method = "parametric_si", mean_si = 2.6),
    "std_si must be supplied"
  )
  expect_error(
    estimate_advantage(incid, priors = priors, method = "parametric_si", mean_si = 1, std_si = 1.5),
    "mean_si"
  )
  expect_error(
    estimate_advantage(incid, priors = priors, method = "parametric_si", mean_si = 2.6, std_si = 0),
    "std_si"
  )
  expect_error(
    estimate_advantage(
      incid,
      priors = priors,
      method = "parametric_si",
      mean_si = c(2.6, 3.1, 3.4),
      std_si = c(1.5, 1.2, 1.1)
    ),
    "length 1 or n_v"
  )
})

test_that("estimate_advantage parametric_si uses the computed t_min when omitted", {
  n_v <- 2
  n_loc <- 2
  T <- 50
  incid <- array(10, dim = c(T, n_loc, n_v))
  mean_si <- 2.6
  std_si <- 1.5

  si_parametric <- matrix(NA_real_, nrow = T, ncol = n_v)
  for (k in seq_len(n_v)) {
    col <- discr_si(seq(0, T - 1), mean_si, std_si)
    si_parametric[, k] <- col / sum(col)
  }
  expected_t_min <- compute_t_min(incid, si_parametric)

  out <- estimate_advantage(
    incid,
    priors = default_priors(),
    mcmc_control = list(n_iter = 15L, burnin = 5L, thin = 5L),
    seed = 1,
    method = "parametric_si",
    mean_si = mean_si,
    std_si = std_si
  )

  expect_true(all(is.na(out$R[seq(1, expected_t_min - 1), , seq_len(dim(out$R)[3])])))
  expect_false(anyNA(out$R[seq(expected_t_min, dim(out$R)[1]), , seq_len(dim(out$R)[3])]))
})



test_that("estimate_advantage convergence checks work with >2 variants", {
  n_v <- 3 # 3 variants
  n_loc <- 4 # 4 locations
  T <- 100 # 100 time steps

  priors <- default_priors()

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v, w_v)
  low_iter <- list(n_iter = 60L, burnin = 10L, thin = 1L)
  x <- estimate_advantage(
    incid, si_distr, priors, seed = 1, t_min = 2L, mcmc_control = low_iter
  )
  ## convergence should be a list of length 2.
  ## not checking whether chains have converged or not.
  ## that is tested in a different set of tests.
  expect_length(x$convergence, 2)
  expect_length(x$diag, 2)
})

test_that("estimate_advantage produces expected warning message", {
  n_v <- 2 # 2 variants
  n_loc <- 1 # 1 locations
  T <- 100 # 100 time steps
  ## Try to use non-default priors
  priors <-   list(epsilon = list(shape = 10, scale = 10),
                   R = list(shape = 0.04, scale = 25))

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v)

  expect_warning(
    estimate_advantage(incid, si_distr, priors, seed = 1, t_min = 2L),
    "Priors where the mean of epsilon is different from 1 are not currently supported."
  )

})
