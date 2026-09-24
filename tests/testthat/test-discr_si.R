si_cases <- list(
  c(2.6, 1.5),  # Flu
  c(4.7, 2.9),  # COVID
  c(8.4, 3.8),  # Ebola-like
  c(14.9, 3.9), # Measles-like
  c(22, 8),     # Smallpox-like
  c(1.2, 0.1),  # Short and narrow
  c(3, 6)       # Long tailed
)

test_that("discr_si matches the original implementation by default", {
  for (x in si_cases) {
    k <- seq(0, 60)
    expect_equal(
      discr_si(k, x[1], x[2]), old_discr_si(k, x[1], x[2]),
      tolerance = 1e-10
    )
  }
})

test_that("discr_si matches the original implementation for unsorted k", {
  k <- c(5, 0, 2, 2, 30, 1)
  expect_equal(
    discr_si(k, 4.7, 2.9), old_discr_si(k, 4.7, 2.9),
    tolerance = 1e-10
  )
})

test_that("discr_si reproduces the documented Flu example", {
  expect_equal(
    round(discr_si(0:4, mu = 2.6, sigma = 1.5), 4),
    round(old_discr_si(0:4, mu = 2.6, sigma = 1.5), 4)
  )
  expect_identical(discr_si(0, 2.6, 1.5), 0)
})

test_that("discr_si keeps the original input checks", {
  expect_error(discr_si(0:5, 2, -1), "sigma must be >=0")
  expect_error(discr_si(0:5, 1, 1), "mu must be >1")
  expect_error(discr_si(-1:5, 2, 1), "all values in k must be >=0")
  expect_error(discr_si(0:5, 4.7, 2.9, shift = -1), "shift must be >=0")
  expect_error(discr_si(0:5, mu = 4.7), "Both mu and sigma")
  expect_error(discr_si(0:5, sigma = 2.9), "Both mu and sigma")
  expect_error(discr_si(0:5, 4.7, 2.9, L = 5, D = 2), "L must be smaller than D")
})

test_that("DiscrSI still returns discr_si output", {
  expect_equal(
    suppressWarnings(DiscrSI(0:10, 2.6, 1.5)), discr_si(0:10, 2.6, 1.5)
  )
})

test_that("discr_si preserves the mean for Gamma and Lognormal mu and sigma", {
  # With a uniform primary event over one day, the discretised delay has the
  # same mean as the continuous delay.
  # A long support is needed for the long tailed case
  k <- seq(0, 2000)
  for (dist in list(stats::pgamma, stats::plnorm)) {
    for (x in si_cases) {
      w <- discr_si(k, x[1], x[2], dist = dist)
      expect_equal(sum(w), 1, tolerance = 1e-6)
      expect_equal(pmf_mean(k, w), x[1], tolerance = 1e-4)
    }
  }
})

test_that("discr_si supports shift = 0", {
  k <- seq(0, 100)
  w <- discr_si(k, 4.7, 2.9, shift = 0)
  expect_gt(w[1], 0)
  expect_equal(sum(w), 1, tolerance = 1e-6)
  expect_equal(pmf_mean(k, w), 4.7, tolerance = 1e-4)
})

test_that("discr_si with no shift matches primarycensored directly", {
  k <- seq(0, 30)
  mu <- 4.7
  sigma <- 2.9
  expected <- primarycensored::dprimarycensored(
    k, stats::pgamma,
    shape = (mu / sigma)^2, scale = sigma^2 / mu,
    pwindow = 1
  )
  expect_equal(discr_si(k, mu, sigma, shift = 0), expected, tolerance = 1e-10)
})

test_that("discr_si left truncation removes mass below L", {
  k <- seq(0, 100)
  w <- discr_si(k, 4.7, 2.9, shift = 0, L = 1)
  expect_identical(w[1], 0)
  expect_equal(sum(w), 1, tolerance = 1e-6)
  untruncated <- discr_si(k, 4.7, 2.9, shift = 0)
  expect_equal(w[-1], untruncated[-1] / sum(untruncated[-1]), tolerance = 1e-8)
})

test_that("discr_si right truncation removes mass at and above D", {
  k <- seq(0, 30)
  w <- discr_si(k, 8.4, 3.8, D = 10)
  expect_true(all(w[k >= 10] == 0))
  expect_equal(sum(w), 1, tolerance = 1e-8)
  untruncated <- discr_si(k, 8.4, 3.8)
  expect_equal(
    w[k < 10], untruncated[k < 10] / sum(untruncated[k < 10]),
    tolerance = 1e-8
  )
})

test_that("discr_si combines shift and both truncations", {
  k <- seq(0, 30)
  w <- discr_si(k, 8.4, 3.8, shift = 1, L = 3, D = 12)
  expect_true(all(w[k < 3 | k >= 12] == 0))
  expect_equal(sum(w), 1, tolerance = 1e-8)
})

test_that("discr_si accepts a CDF function with its own parameters", {
  k <- seq(0, 40)
  mu <- 4.7
  sigma <- 2.9
  shape <- ((mu - 1) / sigma)^2
  scale <- sigma^2 / (mu - 1)
  expect_equal(
    discr_si(k, dist = stats::pgamma, shape = shape, scale = scale),
    discr_si(k, mu, sigma),
    tolerance = 1e-10
  )
})

test_that("discr_si errors when both mu/sigma and a CDF function are given", {
  expect_error(
    discr_si(0:5, 2, 1, dist = stats::pgamma, shape = 2, scale = 1),
    "not both"
  )
})

test_that("discr_si only takes mu and sigma for Gamma and Lognormal", {
  expect_error(discr_si(0:5, 2, 1, dist = stats::pweibull), "mu and sigma")
  expect_error(discr_si(0:5, 2, 1, dist = "gamma"), "dist")
})

test_that("discr_si passes native parameters to other distributions", {
  k <- seq(0, 200)
  w <- discr_si(k, dist = stats::pweibull, shape = 1.5, scale = 4)
  expect_equal(sum(w), 1, tolerance = 1e-6)
  expect_equal(pmf_mean(k, w), 1 + 4 * gamma(1 + 1 / 1.5), tolerance = 1e-4)
})

test_that("discr_si supports other primary event distributions", {
  k <- seq(0, 100)
  w_unif <- discr_si(k, 4.7, 2.9)
  w_growth <- discr_si(
    k, 4.7, 2.9,
    dprimary = primarycensored::dexpgrowth, primary_args = list(r = 0.5)
  )
  expect_equal(sum(w_growth), 1, tolerance = 1e-6)
  # Growth puts more primary events late in the day, so secondary events
  # fall on later days
  expect_gt(pmf_mean(k, w_growth), pmf_mean(k, w_unif))
})

test_that("discr_si agrees with Monte Carlo simulation", {
  set.seed(123)
  n <- 1e5
  mu <- 4.7
  sigma <- 2.9
  sdlog <- sqrt(log(1 + sigma^2 / (mu - 1)^2))
  meanlog <- log(mu - 1) - sdlog^2 / 2
  samples <- 1 + primarycensored::rprimarycensored(
    n,
    rdist = function(n) stats::rlnorm(n, meanlog, sdlog),
    rprimary = primarycensored::rexpgrowth,
    rprimary_args = list(r = 0.2),
    pwindow = 1, swindow = 1, D = Inf
  )
  k <- seq(0, 60)
  empirical <- tabulate(samples + 1, nbins = length(k)) / n
  w <- discr_si(
    k, mu, sigma, dist = stats::plnorm,
    dprimary = primarycensored::dexpgrowth, primary_args = list(r = 0.2)
  )
  expect_lt(max(abs(w - empirical)), 0.005)
})
