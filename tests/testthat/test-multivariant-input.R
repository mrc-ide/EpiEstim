# devtools::load_all()
require(testthat)

################################
## Tests for compute_lambda() ##
################################

n_v <- 2 # 2 variants
n_loc <- 3 # 3 locations
T <- 100 # 100 time steps
incid <- array(10, dim = c(T, n_loc, n_v)) # constant incidence 10 per day everywhere
incid_processed <- process_I_multivariant(incid)

  # si distr tests

sidistr_1=cbind(c(0.1,0.1,0.5,0.3),c(0,0.2,0.5,0.3)) # first SI value not 0
sidistr_2=cbind(c(0,0.1,0.5,0.3),c(0,0.2,0.5,0.3)) # sum of first SI doesn't equal 1
sidistr_3=cbind(c(0,-0.1,0.7,0.4),c(0,0.2,0.5,0.3)) # negative value

test_that("si_distr is specified correctly", {
  expect_error(compute_lambda(incid=incid_processed, si_distr=sidistr_1),
               "Values in the first row of si_distr must be 0")
  expect_error(compute_lambda(incid=incid_processed, si_distr=sidistr_2),
               "The sum of each column in si_distr should be equal to 1")
  expect_error(compute_lambda(incid=incid_processed, si_distr=sidistr_3),
               "si_distr must be >=0")
})

##############################
## Tests for draw_epsilon() ##
##############################

R <- matrix(1, nrow = T, ncol = n_loc)
w_v <- c(0, 0.2, 0.5, 0.3)
si_distr <- cbind(w_v, w_v, w_v)
lambda <- compute_lambda(incid_processed, si_distr)
priors <- default_priors()

  # tmin/tmax tests

tmin1 <- 2.5 # not an integer
tmin2 <- 1L # less than 2
tmin3 <- as.integer(nrow(incid)+1) # greater than nrow(incid)

test_that("tmin and tmax are specified correctly", {
  expect_error(draw_epsilon(R=R, incid=incid, lambda=lambda, priors=priors,
                              t_min = tmin1, t_max = nrow(incid),
                              seed = NULL),
               "t_min and t_max must be integers")
  expect_error(draw_epsilon(R=R, incid=incid, lambda=lambda, priors=priors,
                            t_min = tmin2, t_max = nrow(incid),
                            seed = NULL),
               "t_min and t_max must be >=2")
  expect_error(draw_epsilon(R=R, incid=incid, lambda=lambda, priors=priors,
                            t_min = tmin3, t_max = nrow(incid),
                            seed = NULL),
               "t_min and t_max must be <= nrow(incid)", fixed=TRUE)
})

  # R test

test_that("R is specified correctly",{
  Rneg <- R
  Rneg[1,1] <- -1
  expect_error(draw_epsilon(R=Rneg, incid=incid, lambda=lambda, priors=priors,
                            t_min = 2L, t_max = nrow(incid),
                            seed = NULL),
               "R must be >= 0")
})

  # seed test

seed <- "a"

test_that("seed is specified correctly",{
  expect_error(draw_epsilon(R=R, incid=incid, lambda=lambda, priors=priors,
                            t_min = 2L, t_max = nrow(incid),
                            seed = "a"),
               "seed must be numeric")
})

########################
## Tests for draw_R() ##
########################

epsilon <- 1
incid <- array(10, dim = c(T, n_loc, n_v))
incid_processed <- process_I_multivariant(incid)
w_v <- c(0, 0.2, 0.5, 0.3)
si_distr <- cbind(w_v, w_v)
lambda <- compute_lambda(incid_processed, si_distr)
priors <- default_priors()

  # tmin/max tests

tmin1 <- 2.5 # not an integer
tmin2 <- 1L # less than 2
tmin3 <- as.integer(nrow(incid)+1) # greater than nrow(incid)

test_that("tmin and tmax are specified correctly", {
  expect_error(draw_R(epsilon=epsilon, incid=incid, lambda=lambda, priors=priors,
                            t_min = tmin1, t_max = nrow(incid),
                            seed = NULL),
               "t_min and t_max must be integers")
  expect_error(draw_R(epsilon=epsilon, incid=incid, lambda=lambda, priors=priors,
                            t_min = tmin2, t_max = nrow(incid),
                            seed = NULL),
               "t_min and t_max must be >=2")
  expect_error(draw_R(epsilon=epsilon, incid=incid, lambda=lambda, priors=priors,
                            t_min = tmin3, t_max = nrow(incid),
                            seed = NULL),
               "t_min and t_max must be <= nrow(incid)", fixed=TRUE)
})


  # seed test

seed <- "a"

test_that("seed is specified correctly",{
  expect_error(draw_R(epsilon=epsilon, incid=incid, lambda=lambda, priors=priors,
                            t_min = 2L, t_max = nrow(incid),
                            seed = seed),
               "seed must be numeric")
})


  # epsilon test

epsilon <- -1

test_that("epsilon is specified correctly",{
  expect_error(draw_R(epsilon=epsilon, incid=incid, lambda=lambda, priors=priors,
                      t_min = 2L, t_max = nrow(incid),
                      seed = NULL),
               "epsilon must be > 0")
})


################################
## Tests for estimate_advantage() ##
################################

incid <- array(10, dim = c(T, n_loc, n_v))
incid_processed <- process_I_multivariant(incid)
w_v <- c(0, 0.2, 0.5, 0.3)
si_distr <- cbind(w_v, w_v)
priors <- default_priors()

  # tmin/max tests

tmin1 <- 2.5 # not an integer
tmin2 <- 1L # less than 2
tmin3 <- as.integer(nrow(incid)+1) # greater than nrow(incid)

test_that("tmin and tmax are specified correctly", {
  expect_error(estimate_advantage(incid=incid, si_distr=si_distr, priors=priors,
                              mcmc_control = default_mcmc_controls(),
                              t_min = tmin1, t_max = nrow(incid),
                              seed = NULL),
               "t_min and t_max must be integers")
  expect_error(estimate_advantage(incid=incid, si_distr=si_distr, priors=priors,
                              mcmc_control = default_mcmc_controls(),
                              t_min = tmin2, t_max = nrow(incid),
                              seed = NULL),
               "t_min and t_max must be >=2")
  expect_error(estimate_advantage(incid=incid, si_distr=si_distr, priors=priors,
                              mcmc_control = default_mcmc_controls(),
                              t_min = tmin3, t_max = nrow(incid),
                              seed = NULL),
               "t_min and t_max must be <= nrow(incid)", fixed=TRUE)
})


  # si_distr tests

sidistr_1=cbind(c(0.1,0.1,0.5,0.3),c(0,0.2,0.5,0.3)) # first SI value not 0
sidistr_2=cbind(c(0,0.1,0.5,0.3),c(0,0.2,0.5,0.3)) # sum of first SI doesn't equal 1
sidistr_3=cbind(c(0,-0.1,0.7,0.4),c(0,0.2,0.5,0.3)) # negative value

test_that("si_distr is specified correctly", {
  expect_error(estimate_advantage(incid=incid, si_distr=sidistr_1, priors=priors,
                              mcmc_control = default_mcmc_controls(),
                              t_min = 2L, t_max = nrow(incid),
                              seed = NULL),
               "Values in the first row of si_distr must be 0")
  expect_error(estimate_advantage(incid=incid, si_distr=sidistr_2, priors=priors,
                              mcmc_control = default_mcmc_controls(),
                              t_min = 2L, t_max = nrow(incid),
                              seed = NULL),
               "The sum of each column in si_distr should be equal to 1")
  expect_error(estimate_advantage(incid=incid, si_distr=sidistr_3, priors=priors,
                              mcmc_control = default_mcmc_controls(),
                              t_min = 2L, t_max = nrow(incid),
                              seed = NULL),
               "si_distr must be >=0")
})


  # mcmc_control tests

mcmc_control1 <- function() {
  list(n_iter = -1,
         burnin = 10L,
         thin = 10L)
}

mcmc_control2 <- function() {
  list(n_iter = 1100L,
       burnin = -1,
       thin = 10L)
}

mcmc_control3 <- function() {
  list(n_iter = 1100L,
       burnin = 10L,
       thin = -1)
}

test_that("mcmc_control is specified correctly", {
  expect_error(estimate_advantage(incid=incid, si_distr=si_distr, priors=priors,
                              mcmc_control = mcmc_control1(),
                              t_min = 2L, t_max = nrow(incid),
                              seed = NULL),
               "n_iter in mcmc_control must be a positive integer")
  expect_error(estimate_advantage(incid=incid, si_distr=si_distr, priors=priors,
                              mcmc_control = mcmc_control2(),
                              t_min = 2L, t_max = nrow(incid),
                              seed = NULL),
               "burnin in mcmc_control must be a positive integer")
  expect_error(estimate_advantage(incid=incid, si_distr=si_distr, priors=priors,
                              mcmc_control = mcmc_control3(),
                              t_min = 2L, t_max = nrow(incid),
                              seed = NULL),
               "thin in mcmc_control must be a positive integer")
})


  # seed test

seed <- "a"

test_that("seed is specified correctly",{
  expect_error(estimate_advantage(incid=incid, si_distr=si_distr, priors=priors,
                              mcmc_control = default_mcmc_controls(),
                              t_min = 2L, t_max = nrow(incid),
                              seed = seed),
               "seed must be numeric")
})

  # convergence check
  # Note: At least ~70 iterations, otherwise mcmc chains too short to run
  # gelman diagnostic
  # with set.seed(42) the low_iter MCMC should fail to converge, while the
  # high_iter MCMC should succeed in converging

set.seed(42)

low_iter <- function() {
  list(n_iter = 70L,
       burnin = 10L,
       thin = 10L)
}

high_iter <- function() {
  list(n_iter = 5000L,
       burnin = 2500L,
       thin = 10L)
}

test_that("convergence check returns FALSE for non-converging MCMC", {
  out <- estimate_advantage(incid=incid, si_distr=si_distr, priors=priors,
                        mcmc_control = low_iter(),
                        t_min = 2L, t_max = nrow(incid),
                        seed = NULL)

  expect_false(out$convergence[1])
})

test_that("convergence check returns TRUE for converging MCMC", {
  out <- estimate_advantage(incid=incid, si_distr=si_distr, priors=priors,
                        mcmc_control = high_iter(),
                        t_min = 2L, t_max = nrow(incid),
                        seed = NULL)

  expect_true(out$convergence[1])
})


  # reordering incidence tests

  # use dummy incidence data based on Rt_ref = 1.6, epsilon = 2

incid <- readRDS("../incidence_files/incid.rds")
incid <- as.array(incid[[1]])
incid <- incid[[1]]
si_for_est <- readRDS("../incidence_files/si_for_est.rds")
si_for_est <- si_for_est[[1]]
priors <- default_priors()

test_that("estimate of epsilon is the same with or without incidence reordering", {
  
  out <- estimate_advantage(incid=incid, si_distr=si_for_est, priors=priors,
                            mcmc_control = default_mcmc_controls(),
                            t_min = 2L, t_max = nrow(incid),
                            seed = NULL, reorder_incid = TRUE)
  
  out2 <- estimate_advantage(incid=incid, si_distr=si_for_est, priors=priors,
                            mcmc_control = default_mcmc_controls(),
                            t_min = 2L, t_max = nrow(incid),
                            seed = NULL, reorder_incid = FALSE)
  
  expect_equal(mean(out$epsilon), mean(out2$epsilon), tolerance = 0.001)
})


test_that("the Rt of the reference variant is returned both with and without incidence reordering", {
  
  out <- estimate_advantage(incid=incid, si_distr=si_for_est, priors=priors,
                            mcmc_control = default_mcmc_controls(),
                            t_min = 2L, t_max = nrow(incid),
                            seed = NULL, reorder_incid = TRUE)
  
  out2 <- estimate_advantage(incid=incid, si_distr=si_for_est, priors=priors,
                             mcmc_control = default_mcmc_controls(),
                             t_min = 2L, t_max = nrow(incid),
                             seed = NULL, reorder_incid = FALSE)
  
  expect_equal(median(out$R, na.rm = T), median(out2$R, na.rm = T), tolerance = 0.001)
})