
#' @title Estimated Instantaneous Reproduction Number from coarsely aggregated data
#'
#' @param incid aggregated incidence data, supplied as a vector
#' @param dt length of temporal aggregation of the data (numeric, 7 time units (typically days) by default)
#' @param dt_out length of the sliding windows for R estimates (numeric, 7 time units (typically days) by default)
#' @param iter number of iterations of the EM algorithm (numeric, 10 by default)
#' @param config An object of class \code{estimate_R_config}, as returned by 
#' function \code{make_config}. 
#' @param method One of "non_parametric_si" or "parametric_si" (see details).
#' @param grid named list containing "precision", "min", and "max" which are used to
#' define a grid of growth rate parameters that are used inside the EM algorithm (see details). 
#' We recommend using the default values. 
#'
#' @return {
#' an object of class \code{estimate_R}, with components:
#' \itemize{
#'
#' \item{R}{: a dataframe containing:
#' the times of start and end of each time window considered ;
#' the posterior mean, std, and 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975
#' quantiles of the reproduction number for each time window.}
#'
#' \item{method}{: the method used to estimate R, one of "non_parametric_si",
#' "parametric_si", "uncertain_si", "si_from_data" or "si_from_sample"}
#'
#' \item{si_distr}{: a vector or dataframe (depending on the method) containing
#'  the discrete serial interval distribution(s) used for estimation}
#'
#' \item{SI.Moments}{: a vector or dataframe (depending on the method)
#' containing the mean and std of the discrete serial interval distribution(s)
#' used for estimation}
#'
#' \item{I}{: the time series of total incidence}
#'
#' \item{I_local}{: the time series of incidence of local cases (so that
#' \code{I_local + I_imported = I})}
#'
#' \item{I_imported}{: the time series of incidence of imported cases (so that
#' \code{I_local + I_imported = I})}
#'
#' \item{dates}{: a vector of dates corresponding to the incidence time series}

#' }
#' }
#' @importFrom epitrix gamma_mucv2shapescale r2R0
#' @importFrom distcrete distcrete
#' @export
#'
#' @details TODO: add details
#' 
#' @examples 
#' ## load data on SARS in 2003
#' data("SARS2003")
#' 
#' ## aggregate the data to weekly counts of cases
#' incid <- SARS2003$incidence
#' dt <- 7L
#' weekly_inc <- aggregate_inc(incid, dt)
#' si_distr <- SARS2003$si_distr
#' 
#' ## estimate Rt using the default parameters (method "non_parametric_si")
#' method <- "non_parametric_si"
#' config <- make_config(list(si_distr = si_distr))
#' res_weekly <- estimate_R_agg(incid = weekly_inc, 
#'                             dt = 7, # aggregation window of the data
#'                             dt_out = 7, # desired sliding window length
#'                             iter = 10,
#'                             config = config,
#'                             method = method,
#'                             grid = list(precision = 0.001, min = -1, max = 1))
#' 
#' ## plot the result
#' plot(res_weekly)
estimate_R_agg <- function(incid,
                           dt = 7L, # aggregation window of the data
                           dt_out = 7L, # desired sliding window length
                           iter = 10L,
                           config = make_config(), 
                           method = c("non_parametric_si", "parametric_si"),
                           grid = list(precision = 0.001, min = -1, max = 1)){ 
  
  if (!is.integer(dt) | !is.integer(dt_out)) {
    stop("dt and dt_out must be integers e.g. 7L")
  }
  if (!is.integer(iter)) {
    stop("iter must be an integer e.g. 10L")
  }
  if (!is.list(grid) | !length(grid) == 3){
    stop("grid must be a list of 3 elements: precision, min, and max")
  }
  if (grid$max < grid$min){
    stop("grid max must be larger than grid min")
  }
  if (grid$precision > grid$max-grid$min){
    stop("grid precision must be less than grid max - grid min")
  }
  if (!is.numeric(grid$precision) | !is.numeric(grid$min) | !is.numeric(grid$max)){
    stop("grid precision, min, and max, must all be numeric")
  }
  
  # Two configs:
  # 'config' for the R estimates used to reconstruct the incidence (internal to the
  # EM algorithm). These use a fixed window length matched to dt (aggregation window)
  # 'config_out' for the final estimated R using sliding windows
  
  method <- match.arg(method) # potentially add an error message but maybe automatic
  config <- process_config(config)
  check_config(config, method)
  config_out <- config 
  
  if(is.null(config_out$t_start)) {
    config_out$t_start <- seq(from = dt + 1,
                              to = (length(incid) * dt) -
                                (dt_out - 1), 1)
  }
  if(is.null(config_out$t_end)) {
    config_out$t_end <- config_out$t_start + (dt_out - 1)
  }
  
  all_T <- length(incid)*dt
  
  # config$t_start and config$t_end used for the reconstruction. 
  # Rt estimation starts on 1st day of second dt. 
  # Width of fixed time windows match aggregations (dt):
  config$t_start <- seq(from = dt + 1, to = all_T - (dt - 1), dt)
  config$t_end <- seq(from = min(config$t_start) + (dt - 1),to = all_T,dt)
  
  niter <- seq(1,iter,1) 
  sim_inc <- matrix(NA, nrow = all_T, ncol = iter)
  
  for(i in 1:length(niter)){
    if(niter[i]==1){
      
      # Initialisation of EM. Aggregated incidence split evenly:
      dis <- incid/dt
      dis_inc <- rep(dis, each = dt)
      
      # Remember that the EpiEstim SI needs to start at 0
      # vs. projections which needs to start at 1
      R <- estimate_R(dis_inc, 
                      method = method,
                      config = config)
      
      config$t_start <- R$R$t_start
      config$t_end <- R$R$t_end
      
      print(paste0("Estimated R for iteration: ", i))
      Mean_R <- R$R$`Mean(R)`
      
      
      # Translate R to growth rate
      get_r_from_R <- function(R, gt_mean, gt_sd, 
                               gt_distr,
                               grid) {
        r_grid <- seq(grid$min, grid$max, grid$precision)
        if(is.null(gt_distr)) {
          gt_pars <- gamma_mucv2shapescale(mu = gt_mean, cv = gt_sd / gt_mean)
          gt_distr <- distcrete("gamma", interval = 1,
                                           shape = gt_pars$shape,
                                           scale = gt_pars$scale, w = 0.5)
        }
        # using a grid of r values translate that into R using r2R0
        R_grid <- r2R0(r = r_grid, w = gt_distr)
        
        # find location of the value in the R grid which has the smallest 
        # difference to the input of R the user provided e.g. R_grid[idx_r]:
        idx_r <- vapply(R, function(e) which.min(abs(R_grid - e)), numeric(1L)) 
        
        while (any(idx_r == 1) || any(idx_r == length(r_grid))) {
          # rerunning get_r_from_R with a wider r_grid
          grid_multiplier <- 5
          if(grid$max > 0) grid$max <- grid_multiplier * grid$max else grid$max <- - grid$max
          if(grid$min < 0) grid$min <- grid_multiplier * grid$min else grid$min <- - grid$min
          r_grid <- seq(grid$min, grid$max, grid$precision)
          R_grid <- r2R0(r = r_grid, w = gt_distr)
          idx_r <- vapply(R, function(e) which.min(abs(R_grid - e)), numeric(1L))
        }
        r <- vapply(idx_r, function(e) r_grid[e], numeric(1L))
      }
      
      gr <- get_r_from_R(R = Mean_R, 
                         gt_mean = config$mean_si, gt_sd = config$std_si, 
                         gt_distr = config$si_distr,
                         grid = grid)
      
      
      
      # so that the incidence for the first dt can be 
      # reconstructed, making the assumption that the gr
      # for the first week would be the same as the second:
      
      gr <- c(gr[1],gr) 
      
      # 5. Estimate incidence using growth rate
      
      # Assume that It is a constant (k) multiplied by exp(gr[for that dt]*t)
      ngroups <- length(incid)
      d <- numeric(length=ngroups)
      k <- numeric(length=ngroups)
      
      for (w in 1:length(k)){
        d[w] <- sum(exp(gr[w] * seq(1, dt - 1, 1)))
        k[w] <- incid[w] / (exp(gr[w]) * (1 + d[w]))
      }
      
      k_seq <- rep(k, each = dt)
      gr_seq <- rep(gr, each = dt)
      
      est_inc <- numeric(length = all_T)
      days <- seq(1, all_T, 1)
      w_day <- rep(seq(1, dt), ngroups)
      
      # Reconstruct daily incidence
      for(t in seq_along(days)){
        est_inc[t] <- k_seq[t] * exp(gr_seq[t] * w_day[t])
      }
      
      sim_inc[,i] <- est_inc
      
      print(paste0("Reconstructed incidence for iteration: ", i))
      
      
    } else {
      
      # Use new adjusted incidence as starting point
      new_inc <- sim_inc[,i-1]
      
      # Re-Estimate R
      
      R2 <- estimate_R(new_inc,
                       method = method,
                       config = config)
      
      print(paste0("Estimated R for iteration: ", i))
      
      
      Mean_R_2 <- R2$R$`Mean(R)`
      
      # Translate R to growth rate again
      gr2 <- get_r_from_R(R = Mean_R_2, 
                         gt_mean = config$mean_si, gt_sd = config$std_si, 
                         gt_distr = config$si_distr,
                         grid = grid)
      
      # again, making the assumption that the gr for the first 
      # dt would be the same as the second dt:
      gr2 <- c(gr2[1], gr2) 
      
      # Estimate incidence 
      k2 <- numeric(length=ngroups)
      d2 <- numeric(length=ngroups)
      
      for (w in seq_along(k2)){
        d2[w] <- sum(exp(gr2[w] * seq(1, dt - 1, 1)))
        k2[w] <- incid[w] / (exp(gr2[w]) * (1 + d2[w]))
      }
      
      k2_seq <- rep(k2, each = dt)
      gr2_seq <- rep(gr2, each = dt)
      
      est_inc2 <- numeric(length = all_T)
      days <- seq(1, all_T, 1)
      w_day <- rep(seq(1, dt), ngroups)
      
      for(t in 1:length(days)){
        est_inc2[t] <- k2_seq[t] * exp(gr2_seq[t] * w_day[t])
      }
      
      sim_inc[,i] <- est_inc2
      
      # monitor progress:
      print(paste0("Reconstructed incidence for iteration: ", i))
      
      if(niter[i] == max(niter)){
        R_out <- estimate_R(sim_inc[,i],
                            method = method,
                            config = config_out)
      }
      
    }
  }
  
  R_out
  
}
