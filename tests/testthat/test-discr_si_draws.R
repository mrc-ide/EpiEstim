draw_mu <- c(2.6, 4.7, 8.4, 3)
draw_sigma <- c(1.5, 2.9, 3.8, 6)

expect_draws_match <- function(k, si_discr_args = list()) {
  draws <- discr_si_draws(k, draw_mu, draw_sigma, si_discr_args)
  expected <- t(vapply(
    seq_along(draw_mu),
    function(i) {
      do.call(
        discr_si,
        c(list(k = k, mu = draw_mu[i], sigma = draw_sigma[i]), si_discr_args)
      )
    },
    numeric(length(k))
  ))
  testthat::expect_equal(draws, expected, tolerance = 1e-12)
}

test_that("discr_si_draws matches discr_si for each draw by default", {
  expect_draws_match(seq(0, 40))
})

test_that("discr_si_draws matches discr_si for Gamma and Lognormal", {
  for (dist in list(stats::pgamma, stats::plnorm)) {
    expect_draws_match(seq(0, 40), list(dist = dist))
  }
})

test_that("discr_si_draws matches discr_si with shift and truncation", {
  expect_draws_match(seq(0, 40), list(shift = 0, L = 1))
  expect_draws_match(seq(0, 40), list(D = 12))
  expect_draws_match(seq(0, 40), list(shift = 1, L = 3, D = 15))
  expect_draws_match(c(5, 0, 2, 2, 30, 1), list(D = 20))
})

test_that("discr_si_draws matches discr_si with other primary distributions", {
  expect_draws_match(
    seq(0, 30),
    list(dprimary = primarycensored::dexpgrowth, primary_args = list(r = 0.2))
  )
})

test_that("discr_si_draws returns one row per draw", {
  draws <- discr_si_draws(seq(0, 10), draw_mu, draw_sigma)
  expect_identical(dim(draws), c(length(draw_mu), 11L))
})

test_that("discr_si_draws checks its inputs", {
  expect_error(
    discr_si_draws(0:5, c(2, 1), c(1, 1)), "mu must be >1"
  )
  expect_error(
    discr_si_draws(0:5, c(2, 3), c(1, -1)), "sigma must be >=0"
  )
  expect_error(
    discr_si_draws(0:5, 2, 1, list(foo = 1)), "si_discr_args"
  )
})

test_that("discr_si_param_draws matches discr_si for native parameters", {
  k <- seq(0, 30)
  params <- list(list(shape = 1.5, scale = 3), list(shape = 3, scale = 2))
  si_args <- si_discr_defaults(list(dist = stats::pweibull, D = 20))
  draws <- discr_si_param_draws(k, params, si_args)
  expected <- t(vapply(
    params,
    function(p) {
      do.call(discr_si, c(list(k = k, dist = stats::pweibull, D = 20), p))
    },
    numeric(length(k))
  ))
  expect_equal(draws, expected, tolerance = 1e-12)
})
