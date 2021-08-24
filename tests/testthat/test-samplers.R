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
    apply(x, 1, mean), c(1, 1), tolerance = 0.05
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
  expect_equal(apply(x, 1, mean), c(1, 1), tolerance = 0.05)
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
  expect_true(max(abs(x_mean[-c(1, 2, 3), ] - 1)) < 0.05)
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
  expect_true(max(abs(x_mean[-c(1, 2, 3), ] - 1)) < 0.05)
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
  expect_true(max(abs(x_mean[-c(1, 2, 3), ] - 1)) < 0.05)
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
  expect_true(max(abs(x_mean[-c(1, 2, 3), ] - 1)) < 0.05)
})


test_that("estimate_joint produces expected results (2 variants 3 locations)", {
  n_v <- 2 # 2 variants
  n_loc <- 3 # 3 locations
  T <- 100 # 100 time steps

  priors <- default_priors()

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v)

  x <- estimate_joint(incid, si_distr, priors, seed = 1, t_min = 2L)

  ## epsilon should be approximately 1
  expect_equal(mean(x$epsilon), 1, tolerance = 0.05)

  ## R should be approximately 1
  ## not exactly 1 because of the first few timesteps & because of priors
  ## so ignore fisrt timesteps
  mean_R <- apply(x$R, c(1, 2), mean)
  expect_true(max(abs(mean_R[-c(1, 2, 3), ] - 1)) < 0.1)
})


test_that("estimate_joint produces expected results (2 variants 1 location)", {
  n_v <- 2 # 2 variants
  n_loc <- 1 # 1 locations
  T <- 100 # 100 time steps

  priors <- default_priors()

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v)

  x <- estimate_joint(incid, si_distr, priors, seed = 1, t_min = 2L)

  ## epsilon should be approximately 1
  expect_equal(mean(x$epsilon), 1, tolerance = 0.05)

  ## R should be approximately 1
  ## not exactly 1 because of the first few timesteps & because of priors
  ## so ignore fisrt timesteps
  mean_R <- apply(x$R, c(1, 2), mean)
  expect_true(max(abs(mean_R[-c(1, 2, 3), ] - 1)) < 0.1)
})


test_that("estimate_joint produces expected results (>2 variants 4 locs)", {
  n_v <- 3 # 3 variants
  n_loc <- 4 # 4 locations
  T <- 100 # 100 time steps

  priors <- default_priors()

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v, w_v)

  x <- estimate_joint(incid, si_distr, priors, seed = 1, t_min = 2L)

  ## epsilon should be approximately 1
  ## FIXME this should be apply(x$epsilon, 1, mean)
  ## as 2 epsilons are returned here.

  expect_equal(
    apply(x$epsilon, 1, mean), c(1, 1), tolerance = 0.05
  )

  ## R should be approximately 1
  ## not exactly 1 because of the first few timesteps & because of priors
  ## so ignore fisrt timesteps
  mean_R <- apply(x$R, c(1, 2), mean)
  expect_true(max(abs(mean_R[-c(1, 2, 3), ] - 1)) < 0.1)
})


test_that("estimate_joint produces expected results (>2 variants 1 loc)", {
  n_v <- 3 # 3 variants
  n_loc <- 1 # 1 location
  T <- 100 # 100 time steps

  priors <- default_priors()

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v, w_v)

  x <- estimate_joint(incid, si_distr, priors, seed = 1, t_min = 2L)

  ## epsilon should be approximately 1
  expect_equal(
    apply(x$epsilon, 1, mean), c(1, 1), tolerance = 0.05
  )

  ## R should be approximately 1
  ## not exactly 1 because of the first few timesteps & because of priors
  ## so ignore fisrt timesteps
  mean_R <- apply(x$R, c(1, 2), mean)
  expect_true(max(abs(mean_R[-c(1, 2, 3), ] - 1)) < 0.1)
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


test_that("estimate_joint produces expected results (>2var, 1loc, imports)", {
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

  x <- estimate_joint(incid, si_distr, priors, seed = 1,
                      incid_imported = incid_imported, t_min = 2L)

  ## epsilon should be approximately 0
  expect_equal(
    apply(x$epsilon, 1, mean), c(0, 0), tolerance = 0.05
  )

  ## R should be approximately 1
  ## not exactly 1 because of the first few timesteps & because of priors
  ## so ignore fisrt timesteps
  mean_R <- apply(x$R, c(1, 2), mean)
  expect_true(max(abs(mean_R[-c(1, 2, 3), ] - 1)) < 0.1)
})


test_that("estimate_joint produces expected results (>2var, 4loc, imports)", {
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
  x <- estimate_joint(
    incid, si_distr, priors, seed = 1,
    incid_imported = incid_imported,
    mcmc_control = list(n_iter = 2000L, burnin = 100L, thin = 10L), t_min = 2L
  )

  ## epsilon should be approximately 0
  expect_equal(
    apply(x$epsilon, 1, mean), c(0, 0), tolerance = 0.05
  )

  ## R should be approximately 1
  ## not exactly 1 because of the first few timesteps & because of priors
  ## so ignore fisrt timesteps
  mean_R <- apply(x$R, c(1, 2), mean)
  expect_true(max(abs(mean_R[-c(1, 2, 3), ] - 1)) < 0.1)
})


test_that("estimate_joint produces expected results (2var, 1loc, imports)", {
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

  x <- estimate_joint(incid, si_distr, priors, seed = 1,
                      incid_imported = incid_imported, t_min = 2L)

  ## epsilon should be approximately 0
  expect_equal(mean(x$epsilon), 0, tolerance = 0.05)

  ## R should be approximately 1
  ## not exactly 1 because of the first few timesteps & because of priors
  ## so ignore fisrt timesteps
  mean_R <- apply(x$R, c(1, 2), mean)
  expect_true(max(abs(mean_R[-c(1, 2, 3), ] - 1)) < 0.1)
})


test_that("estimate_joint produces expected results (2var, 4loc, imports)", {
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

  x <- estimate_joint(
    incid, si_distr, priors, seed = 1,
    incid_imported = incid_imported,
    mcmc_control = list(n_iter = 2000L, burnin = 100L, thin = 10L), t_min = 2L
  )

  ## epsilon should be approximately 0
  expect_equal(mean(x$epsilon), 0, tolerance = 0.05)

  ## R should be approximately 1
  ## not exactly 1 because of the first few timesteps & because of priors
  ## so ignore fisrt timesteps
  mean_R <- apply(x$R, c(1, 2), mean)
  expect_true(max(abs(mean_R[-c(1, 2, 3), ] - 1)) < 0.1)
})


test_that("estimate_joint produces expected results (2 var, 2 loc, R_loc1 = 1.1, R_loc2 = 1.5)", {
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
  x <- estimate_joint(
    incid, si_distr, priors, seed = 1, t_min = 2L
  )

  ## R should be approx 1.1 for loc1 and 1.5 for loc2
  expect_equal(mean(x$R[,1,], na.rm=T), 1.1, tolerance = 0.5)
  expect_equal(mean(x$R[,2,], na.rm=T), 1.5, tolerance = 0.5)

})

test_that("estimate_joint faster with precompute (2 variants 3 locations)", {
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
    x1 <- estimate_joint(incid, si_distr, priors, seed = 1, precompute = TRUE, t_min = 2L)
  )

  t2 <- system.time(
    x2 <- estimate_joint(incid, si_distr, priors, seed = 1, precompute = FALSE, t_min = 2L)
  )

  ## t1 should be < t2
  expect_true(t1[["elapsed"]] < t2[["elapsed"]])

  ## epsilon should be approximately 1 in both cases
  expect_equal(mean(x1$epsilon), 1, tolerance = 0.05)
  expect_equal(mean(x2$epsilon), 1, tolerance = 0.05)

  ## R should be approximately 1 in both cases
  ## not exactly 1 because of the first few timesteps & because of priors
  ## so ignore fisrt timesteps
  mean_R1 <- apply(x1$R, c(1, 2), mean)
  expect_true(max(abs(mean_R1[-c(1, 2, 3), ] - 1)) < 0.1)
  mean_R2 <- apply(x2$R, c(1, 2), mean)
  expect_true(max(abs(mean_R2[-c(1, 2, 3), ] - 1)) < 0.1)
})


test_that("estimate_joint faster with precompute (2 variants 1 location)", {
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
    x1 <- estimate_joint(incid, si_distr, priors, seed = 1, precompute = TRUE, t_min = 2L)
  )

  t2 <- system.time(
    x2 <- estimate_joint(incid, si_distr, priors, seed = 1, precompute = FALSE, t_min = 2L)
  )

  ## t1 should be < t2
  expect_true(t1[["elapsed"]] < t2[["elapsed"]])

  ## epsilon should be approximately 1 in both cases
  expect_equal(mean(x1$epsilon), 1, tolerance = 0.05)
  expect_equal(mean(x2$epsilon), 1, tolerance = 0.05)

  ## R should be approximately 1 in both cases
  ## not exactly 1 because of the first few timesteps & because of priors
  ## so ignore fisrt timesteps
  mean_R1 <- apply(x1$R, c(1, 2), mean)
  expect_true(max(abs(mean_R1[-c(1, 2, 3), ] - 1)) < 0.1)
  mean_R2 <- apply(x2$R, c(1, 2), mean)
  expect_true(max(abs(mean_R2[-c(1, 2, 3), ] - 1)) < 0.1)
})


test_that("estimate_joint faster with precompute (3 variants 4 locations)", {
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
    x1 <- estimate_joint(incid, si_distr, priors, seed = 1, precompute = TRUE, t_min = 2L)
  )

  t2 <- system.time(
    x2 <- estimate_joint(incid, si_distr, priors, seed = 1, precompute = FALSE, t_min = 2L)
  )

  ## t1 should be < t2
  expect_true(t1[["elapsed"]] < t2[["elapsed"]])

    ## epsilon should be approximately 1 in both cases
  expect_equal(mean(x1$epsilon), 1, tolerance = 0.05)
  expect_equal(mean(x2$epsilon), 1, tolerance = 0.05)

  ## R should be approximately 1 in both cases
  ## not exactly 1 because of the first few timesteps & because of priors
  ## so ignore fisrt timesteps
  mean_R1 <- apply(x1$R, c(1, 2), mean)
  expect_true(max(abs(mean_R1[-c(1, 2, 3), ] - 1)) < 0.1)
  mean_R2 <- apply(x2$R, c(1, 2), mean)
  expect_true(max(abs(mean_R2[-c(1, 2, 3), ] - 1)) < 0.1)
})


test_that("estimate_joint faster with precompute (3 variants 1 location)", {
  n_v <- 3 # 3 variants
  n_loc <- 1 # 1 locations
  T <- 100 # 100 time steps

  priors <- default_priors()

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v, w_v)

  t1 <- system.time(
    x1 <- estimate_joint(incid, si_distr, priors, seed = 1, precompute = TRUE, t_min = 2L)
  )

  t2 <- system.time(
    x2 <- estimate_joint(incid, si_distr, priors, seed = 1, precompute = FALSE, t_min = 2L)
  )

  ## t1 should be < t2
  expect_true(t1[["elapsed"]] < t2[["elapsed"]])

  ## epsilon should be approximately 1 in both cases
  expect_equal(mean(x1$epsilon), 1, tolerance = 0.05)
  expect_equal(mean(x2$epsilon), 1, tolerance = 0.05)

  ## R should be approximately 1 in both cases
  ## not exactly 1 because of the first few timesteps & because of priors
  ## so ignore fisrt timesteps
  mean_R1 <- apply(x1$R, c(1, 2), mean)
  expect_true(max(abs(mean_R1[-c(1, 2, 3), ] - 1)) < 0.1)
  mean_R2 <- apply(x2$R, c(1, 2), mean)
  expect_true(max(abs(mean_R2[-c(1, 2, 3), ] - 1)) < 0.1)
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

test_that("estimate_joint uses the correct t_min", {
  n_v <- 2 # 2 variants
  n_loc <- 3 # 3 locations
  T <- 100 # 100 time steps

  priors <- default_priors()

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v)

  x <- estimate_joint(incid, si_distr, priors, t_min = 2L, seed = 1)
  ## If t_min is 2, the first row if x$R will be NA
  expect_true(all(is.na(x$R[1, , ])))
  ## and not after that.
  expect_true(! any(is.na(x$R[seq(2, dim(x$R)[1]), , ])))
  ## if t_min is NULL, t_min would be set to
  ## compute_t_min.
  t_min <- compute_t_min(incid, si_distr)
  x <- estimate_joint(incid, si_distr, priors, seed = 1)
  expect_true(all(is.na(x$R[seq(1, t_min - 1, 1), , ])))
  expect_true(! any(is.na(x$R[seq(t_min, dim(x$R)[1]), , ])))
})



test_that("estimate_joint convergence checks work with >2 variants", {
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
  x <- estimate_joint(
    incid, si_distr, priors, seed = 1, t_min = 2L, mcmc_control = low_iter
  )
  ## convergence should be a list of length 2.
  ## not checking whether chains have converged or not.
  ## that is tested in a different set of tests.
  expect_equal(length(x$convergence), 2)
  expect_equal(length(x$diag), 2)
})
