# Simulate infector/infected pairs with daily onset dates.
# Infector onsets fall on day EL (ER = EL + 1) with day index weighted by
# exp(growth * day), and the serial interval is Gamma distributed with mean
# mu and sd sigma, plus shift. Pairs whose infected onset falls after
# obs_time are not observed.
simulate_si_data <- function(n, mu, sigma, shift = 0, obs_time = 40,
                             growth = 0, seed = 1) {
  set.seed(seed)
  days <- seq(0, obs_time)
  n_sim <- 20 * n
  EL <- sample(days, n_sim, replace = TRUE, prob = exp(growth * days))
  si <- shift + stats::rgamma(
    n_sim, shape = (mu / sigma)^2, scale = sigma^2 / mu
  )
  SL <- floor(EL + stats::runif(n_sim) + si)
  keep <- which(SL <= obs_time)[seq_len(n)]
  data.frame(
    EL = as.integer(EL[keep]), ER = as.integer(EL[keep] + 1),
    SL = as.integer(SL[keep]), SR = as.integer(SL[keep] + 1),
    type = 0L
  )
}

si_from_data_config <- function(incid, ...) {
  make_config(incid = incid, list(
    n1 = 100, n2 = 10, seed = 2,
    mcmc_control = make_mcmc_control(seed = 1),
    ...
  ))
}

# The MockRotavirus dataset
mock_rotavirus <- function() {
  env <- new.env()
  utils::data("MockRotavirus", package = "EpiEstim", envir = env)
  env$MockRotavirus
}

# si_from_data_config for the MockRotavirus incidence without messages
quiet_config <- function(...) {
  suppressMessages(si_from_data_config(mock_rotavirus()$incidence, ...))
}

# Run estimate_R with method "si_from_data" quietly, keeping warnings
run_si_from_data <- function(si_data, config, incid = NULL) {
  if (is.null(incid)) {
    incid <- mock_rotavirus()$incidence
  }
  utils::capture.output(
    res <- estimate_R(
      incid, method = "si_from_data", si_data = si_data, config = config
    )
  )
  res
}

# Mean of the discretised serial interval used for non-offset distributions
discr_mean <- function(mu, sigma, shift = 0, L = 1) {
  k <- seq(0, 200)
  sum(k * discr_si(k, mu, sigma, shift = shift, L = L))
}

# Evaluate expr and return the messages of all warnings raised
collect_warnings <- function(expr) {
  seen <- new.env()
  seen$messages <- character()
  withCallingHandlers(
    suppressMessages(expr),
    warning = function(w) {
      seen$messages <- c(seen$messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  seen$messages
}
