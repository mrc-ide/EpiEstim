context("Gibbs samplers")

test_that("draw_epsilon produces expected results (2 variants, 4 locations)", {
  n_v <- 2 # 2 variants
  n_loc <- 4 # 4 locations
  T <- 100 # 100 time steps

  priors <- EpiEstim:::default_priors()

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))
  incid <- EpiEstim:::process_I_multivariant(incid)

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v)
  lambda <- EpiEstim:::compute_lambda(incid, si_distr)

  # Constant reproduction number of 1
  R <- matrix(1, nrow = T, ncol = n_loc)
  R[1, ] <- NA # no estimates of R on first time step

  set.seed(1)
  x <- sapply(1:1000, function(e) EpiEstim:::draw_epsilon(R, incid$local, lambda, priors))

  ## epsilon should be approximately 1
  ## not exactly 1 because of the first few timesteps & because of priors
  expect_equal(mean(x), 1, tolerance = 0.05)
})


test_that("draw_epsilon produces expected results (2 variants, 1 location)", {
  n_v <- 2 # 2 variants
  n_loc <- 1 # 1 location
  T <- 100 # 100 time steps
  
  priors <- EpiEstim:::default_priors()
  
  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))
  incid <- EpiEstim:::process_I_multivariant(incid)
  
  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v)
  lambda <- EpiEstim:::compute_lambda(incid, si_distr)
  
  # Constant reproduction number of 1
  R <- matrix(1, nrow = T, ncol = n_loc)
  R[1, ] <- NA # no estimates of R on first time step
  
  set.seed(1)
  x <- sapply(1:1000, function(e) EpiEstim:::draw_epsilon(R, incid$local, lambda, priors))
  
  ## epsilon should be approximately 1
  ## not exactly 1 because of the first few timesteps & because of priors
  expect_equal(mean(x), 1, tolerance = 0.05)
})


test_that("draw_epsilon produces expected results (>2 variants, 4 locations)", {
  n_v <- 3 # 3 variants
  n_loc <- 4 # 4 locations
  T <- 100 # 100 time steps

  priors <- EpiEstim:::default_priors()

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))
  incid <- EpiEstim:::process_I_multivariant(incid)

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v, w_v)
  lambda <- EpiEstim:::compute_lambda(incid, si_distr)

  # Constant reproduction number of 1
  R <- matrix(1, nrow = T, ncol = n_loc)
  R[1, ] <- NA # no estimates of R on first time step

  set.seed(1)
  x <- sapply(1:1000, function(e) EpiEstim:::draw_epsilon(R, incid$local, lambda, priors))

  ## epsilon should be approximately 1
  ## not exactly 1 because of the first few timesteps & because of priors
  expect_equal(mean(x), 1, tolerance = 0.05)
})


test_that("draw_epsilon produces expected results (>2 variants, 1 location)", {
  n_v <- 3 # 3 variants
  n_loc <- 1 # 1 location
  T <- 100 # 100 time steps
  
  priors <- EpiEstim:::default_priors()
  
  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))
  incid <- EpiEstim:::process_I_multivariant(incid)
  
  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v, w_v)
  lambda <- EpiEstim:::compute_lambda(incid, si_distr)
  
  # Constant reproduction number of 1
  R <- matrix(1, nrow = T, ncol = n_loc)
  R[1, ] <- NA # no estimates of R on first time step
  
  set.seed(1)
  x <- sapply(1:1000, function(e) EpiEstim:::draw_epsilon(R, incid$local, lambda, priors))
  
  ## epsilon should be approximately 1
  ## not exactly 1 because of the first few timesteps & because of priors
  expect_equal(mean(x), 1, tolerance = 0.05)
})


test_that("draw_R produces expected results (2 variants, 4 locations)", {
  n_v <- 2 # 2 variants
  n_loc <- 4 # 4 locations
  T <- 100 # 100 time steps

  priors <- EpiEstim:::default_priors()

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))
  incid <- EpiEstim:::process_I_multivariant(incid)

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v)
  lambda <- EpiEstim:::compute_lambda(incid, si_distr)

  # Epsilon = 1 i.e. no transmission advantage
  epsilon <- 1

  set.seed(1)
  x <- lapply(1:1000, function(e) EpiEstim:::draw_R(epsilon, incid$local, lambda, priors))
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
  
  priors <- EpiEstim:::default_priors()
  
  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))
  incid <- EpiEstim:::process_I_multivariant(incid)
  
  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v)
  lambda <- EpiEstim:::compute_lambda(incid, si_distr)
  
  # Epsilon = 1 i.e. no transmission advantage
  epsilon <- 1
  
  set.seed(1)
  x <- lapply(1:1000, function(e) EpiEstim:::draw_R(epsilon, incid$local, lambda, priors))
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

  priors <- EpiEstim:::default_priors()

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))
  incid <- EpiEstim:::process_I_multivariant(incid)

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v, w_v)
  lambda <- EpiEstim:::compute_lambda(incid, si_distr)

  # Epsilon = 1 i.e. no transmission advantage
  epsilon <- c(1, 1)

  set.seed(1)
  x <- lapply(1:1000, function(e) EpiEstim:::draw_R(epsilon, incid$local, lambda, priors))
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
  
  priors <- EpiEstim:::default_priors()
  
  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))
  incid <- EpiEstim:::process_I_multivariant(incid)
  
  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v, w_v)
  lambda <- EpiEstim:::compute_lambda(incid, si_distr)
  
  # Epsilon = 1 i.e. no transmission advantage
  epsilon <- c(1, 1)
  
  set.seed(1)
  x <- lapply(1:1000, function(e) EpiEstim:::draw_R(epsilon, incid$local, lambda, priors))
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

  priors <- EpiEstim:::default_priors()

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v)

  x <- EpiEstim:::estimate_joint(incid, si_distr, priors, seed = 1)

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
  
  priors <- EpiEstim:::default_priors()
  
  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))
  
  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v)
  
  x <- EpiEstim:::estimate_joint(incid, si_distr, priors, seed = 1)
  
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

  priors <- EpiEstim:::default_priors()

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v, w_v)

  x <- EpiEstim:::estimate_joint(incid, si_distr, priors, seed = 1)

  ## epsilon should be approximately 1
  expect_equal(mean(x$epsilon), 1, tolerance = 0.05)

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
  
  priors <- EpiEstim:::default_priors()
  
  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))
  
  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v, w_v)
  
  x <- EpiEstim:::estimate_joint(incid, si_distr, priors, seed = 1)
  
  ## epsilon should be approximately 1
  expect_equal(mean(x$epsilon), 1, tolerance = 0.05)
  
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
  
  expect_error(EpiEstim:::process_I_multivariant(incid, incid_imported[-1, , ]),
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
  incid_processed <- EpiEstim:::process_I_multivariant(incid, incid_imported)
  expect_equal(incid_processed$local + incid_processed$imported, incid)
  expect_equal(incid_processed$imported, incid_imported)
  
  ## with the default
  incid_processed <- EpiEstim:::process_I_multivariant(incid)
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
  
  expect_error(EpiEstim:::compute_lambda(incid, si_distr),
      "'incid 'should be an 'incid_multivariant' object.")
})


test_that("estimate_joint produces expected results (>2var, 1loc, imports)", {
  n_v <- 3 # 3 variants
  n_loc <- 1 # 1 location
  T <- 100 # 100 time steps
  
  priors <- EpiEstim:::default_priors()
  
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
  
  x <- EpiEstim:::estimate_joint(incid, si_distr, priors, seed = 1,
                      incid_imported = incid_imported)
  
  ## epsilon should be approximately 0
  expect_equal(mean(1/x$epsilon), 0, tolerance = 0.05)
  
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
  
  priors <- EpiEstim:::default_priors()
  
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
  
  x <- EpiEstim:::estimate_joint(incid, si_distr, priors, seed = 1,
                      incid_imported = incid_imported)
  
  ## epsilon should be approximately 0
  expect_equal(mean(x$epsilon), 0, tolerance = 0.05)
  
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
  
  priors <- EpiEstim:::default_priors()
  
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
  
  x <- EpiEstim:::estimate_joint(incid, si_distr, priors, seed = 1,
                      incid_imported = incid_imported)
  
  ## epsilon should be approximately 0
  expect_equal(mean(1/x$epsilon), 0, tolerance = 0.05)
  
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
  
  priors <- EpiEstim:::default_priors()
  
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
  
  x <- EpiEstim:::estimate_joint(incid, si_distr, priors, seed = 1,
                      incid_imported = incid_imported)
  
  ## epsilon should be approximately 0
  expect_equal(mean(1/x$epsilon), 0, tolerance = 0.05)
  
  ## R should be approximately 1
  ## not exactly 1 because of the first few timesteps & because of priors
  ## so ignore fisrt timesteps
  mean_R <- apply(x$R, c(1, 2), mean)
  expect_true(max(abs(mean_R[-c(1, 2, 3), ] - 1)) < 0.1)
})


test_that("estimate_joint produces expected results (2 var, 2 loc, R_L1 = 1, RL2 = 2)", {
  n_v <- 2 # 2 variants
  n_loc <- 3 # 3 locations
  T <- 100 # 100 time steps
  
  priors <- EpiEstim:::default_priors()
  
  # incidence for R in location 1 = 1
  
  # incidence for R in location 2 = 2
  
  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))
  
  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v)
  
  x <- EpiEstim:::estimate_joint(incid, si_distr, priors, seed = 1)
  
  ## epsilon should be approximately 1
  expect_equal(mean(x$epsilon), 1, tolerance = 0.05)
  
  ## R should be approximately 1
  ## not exactly 1 because of the first few timesteps & because of priors
  ## so ignore fisrt timesteps
  mean_R <- apply(x$R, c(1, 2), mean)
  expect_true(max(abs(mean_R[-c(1, 2, 3), ] - 1)) < 0.1)
})
