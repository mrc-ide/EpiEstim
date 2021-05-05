context("Gibbs samplers")

test_that("draw_epsilon produces expected results (2 variants, 4 locations)", {
  n_v <- 2 # 2 variants
  n_loc <- 4 # 4 locations
  T <- 100 # 100 time steps

  priors <- default_priors()

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v)
  lambda <- compute_lambda(incid, si_distr)

  # Constant reproduction number of 1
  R <- matrix(1, nrow = T, ncol = n_loc)
  R[1, ] <- NA # no estimates of R on first time step

  set.seed(1)
  x <- sapply(1:1000, function(e) draw_epsilon(R, incid, lambda, priors))

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
  
  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v)
  lambda <- compute_lambda(incid, si_distr)
  
  # Constant reproduction number of 1
  R <- matrix(1, nrow = T, ncol = n_loc)
  R[1, ] <- NA # no estimates of R on first time step
  
  set.seed(1)
  x <- sapply(1:1000, function(e) draw_epsilon(R, incid, lambda, priors))
  
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

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v, w_v)
  lambda <- compute_lambda(incid, si_distr)

  # Constant reproduction number of 1
  R <- matrix(1, nrow = T, ncol = n_loc)
  R[1, ] <- NA # no estimates of R on first time step

  set.seed(1)
  x <- sapply(1:1000, function(e) draw_epsilon(R, incid, lambda, priors))

  ## epsilon should be approximately 1
  ## not exactly 1 because of the first few timesteps & because of priors
  expect_equal(mean(x), 1, tolerance = 0.05)
})


test_that("draw_epsilon produces expected results (>2 variants, 1 location)", {
  n_v <- 3 # 3 variants
  n_loc <- 1 # 1 location
  T <- 100 # 100 time steps
  
  priors <- default_priors()
  
  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))
  
  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v, w_v)
  lambda <- compute_lambda(incid, si_distr)
  
  # Constant reproduction number of 1
  R <- matrix(1, nrow = T, ncol = n_loc)
  R[1, ] <- NA # no estimates of R on first time step
  
  set.seed(1)
  x <- sapply(1:1000, function(e) draw_epsilon(R, incid, lambda, priors))
  
  ## epsilon should be approximately 1
  ## not exactly 1 because of the first few timesteps & because of priors
  expect_equal(mean(x), 1, tolerance = 0.05)
})


test_that("draw_R produces expected results (2 variants, 4 locations)", {
  n_v <- 2 # 2 variants
  n_loc <- 4 # 4 locations
  T <- 100 # 100 time steps

  priors <- default_priors()

  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v)
  lambda <- compute_lambda(incid, si_distr)

  # Epsilon = 1 i.e. no transmission advantage
  epsilon <- 1

  set.seed(1)
  x <- lapply(1:1000, function(e) draw_R(epsilon, incid, lambda, priors))
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
  
  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v)
  lambda <- compute_lambda(incid, si_distr)
  
  # Epsilon = 1 i.e. no transmission advantage
  epsilon <- 1
  
  set.seed(1)
  x <- lapply(1:1000, function(e) draw_R(epsilon, incid, lambda, priors))
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

  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v, w_v)
  lambda <- compute_lambda(incid, si_distr)

  # Epsilon = 1 i.e. no transmission advantage
  epsilon <- c(1, 1)

  set.seed(1)
  x <- lapply(1:1000, function(e) draw_R(epsilon, incid, lambda, priors))
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
  
  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v, w_v)
  lambda <- compute_lambda(incid, si_distr)
  
  # Epsilon = 1 i.e. no transmission advantage
  epsilon <- c(1, 1)
  
  set.seed(1)
  x <- lapply(1:1000, function(e) draw_R(epsilon, incid, lambda, priors))
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
  
  priors <- default_priors()
  
  # constant incidence 10 per day everywhere
  incid <- array(10, dim = c(T, n_loc, n_v))
  
  # arbitrary serial interval
  w_v <- c(0, 0.2, 0.5, 0.3)
  si_distr <- cbind(w_v, w_v, w_v)
  
  x <- estimate_joint(incid, si_distr, priors, seed = 1)
  
  ## epsilon should be approximately 1
  expect_equal(mean(x$epsilon), 1, tolerance = 0.05)
  
  ## R should be approximately 1
  ## not exactly 1 because of the first few timesteps & because of priors
  ## so ignore fisrt timesteps
  mean_R <- apply(x$R, c(1, 2), mean)
  expect_true(max(abs(mean_R[-c(1, 2, 3), ] - 1)) < 0.1)
})
