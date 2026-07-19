n_v <- 2
n_loc <- 3 # 3 locations
T <- 100 # 100 time steps
priors <- default_priors()
# constant incidence 10 per day everywhere
incid <- array(10, dim = c(T, n_loc, n_v))
incid <- process_I_multivariant(incid)
# arbitrary serial interval, same for both variants
w_v <- c(0, 0.2, 0.5, 0.3)
si_distr <- cbind(w_v, w_v)
lambda <- compute_lambda(incid, si_distr)
