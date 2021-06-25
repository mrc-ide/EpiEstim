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
  x <- lapply(1:1000, function(e) draw_R(epsilon, incid$local, lambda, priors))
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
  x <- lapply(1:1000, function(e) draw_R(epsilon, incid$local, lambda, priors))
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
  x <- lapply(1:1000, function(e) draw_R(epsilon, incid$local, lambda, priors))
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
  x <- lapply(1:1000, function(e) draw_R(epsilon, incid$local, lambda, priors))
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

  x <- estimate_joint(incid, si_distr, priors, seed = 1)

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

  x <- estimate_joint(incid, si_distr, priors, seed = 1)

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

  x <- estimate_joint(incid, si_distr, priors, seed = 1)

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

  x <- estimate_joint(incid, si_distr, priors, seed = 1)

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
                      incid_imported = incid_imported)

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
    mcmc_control = list(n_iter = 2000L, burnin = 100L, thin = 10L)
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
                      incid_imported = incid_imported)

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
    mcmc_control = list(n_iter = 2000L, burnin = 100L, thin = 10L)
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
            si = si_distr[, v],
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
    incid, si_distr, priors, seed = 1
  )

  ## R should be approx 1.1 for loc1 and 1.5 for loc2
  expect_equal(mean(x$R[,1,], na.rm=T), 1.1, tolerance = 0.5)
  expect_equal(mean(x$R[,2,], na.rm=T), 1.5, tolerance = 0.5)

})


test_that("estimate_joint produces expected results (2 variants 1 location, t_min NULL)", {
  n_v <- 2 # 2 variants
  n_loc <- 1 # 1 locations
  T <- 100 # 100 time steps

  priors <- default_priors()

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))

  ## Reasonable serial interval so that
  ## t_min is sensibly calculated
  w_v1 <- discr_si(0:30, 7, 5)
  ## Different mean for 2nd variant,
  w_v2 <- discr_si(0:30, 14, 5)
  si_distr <- cbind(w_v1, w_v2)

  x <- estimate_joint(incid, si_distr, priors, seed = 1, t_min = NULL)

  ## epsilon should be approximately 1
  expect_equal(mean(x$epsilon), 1, tolerance = 0.05)

  ## R should be approximately 1
  ## not exactly 1 because of the first few timesteps & because of priors
  ## so ignore fisrt timesteps
  mean_R <- apply(x$R, c(1, 2), mean)
  ## if t_min is bigger than 2 the first t_min - 1
  ## rows will be NA
  t_min <- compute_si_cutoff(si_distr)
  expect_true(max(abs(mean_R[-seq(1, t_min - 1, 1), ] - 1)) < 0.1)
})
