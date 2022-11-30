
#' @title Estimated Instantaneous Reproduction Number from coarsely aggregated data
#'
#' @param incid aggregated incidence data, supplied as a vector
#' 
#' @param dt length of temporal aggregations of the data (integer or vector of 
#' integers, integer of 7 time units (typically days) by default - see details).
#' 
#' @param dt_out length of the sliding windows for R estimates (integer, 7 time 
#' units (typically days) by default).
#' 
#' @param iter number of iterations of the EM algorithm (integer, 10 by default)
#' 
#' @param config An object of class \code{estimate_R_config}, as returned by 
#' function \code{make_config}. 
#' 
#' @param method One of "non_parametric_si" or "parametric_si" (see details).
#' 
#' @param grid named list containing "precision", "min", and "max" which are used to
#' define a grid of growth rate parameters that are used inside the EM algorithm 
#' (see details). We recommend using the default values. 
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
#' \item{I}{: the time series of daily incidence reconstructed by the EM algorithm.
#' For the initial incidence that cannot be reconstructed (e.g. the first aggregation
#' window and aggregation windows where incidence is too low to estimate Rt - see 
#' details) then the incidence returned will be the naive disaggregation of the 
#' incidence data used to initialise the EM algorithm. }
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
#' @details 
#' Estimation of the time-varying reproduction number from temporally-aggregated
#' incidence data. For full details about how Rt is estimated within this 
#' bayesian framework, see details in \code{\link{estimate_R}}.
#' 
#' Here, an expectation maximisation (EM) algorithm is used to reconstruct daily 
#' incidence from data that is provided on another timescale. In addition to the 
#' usual parameters required in \code{estimate_R}, the \code{estimate_R_agg()} 
#' function requires that the user specify two additional parameters: dt and 
#' dt_out. There are two other parameters that the user may modify, iter and 
#' grid, however we recommend that the default values are used.
#' 
#' \itemize{
#' \item{dt is the length of the temporal aggregations of the incidence data in 
#' days. This can be a single integer if the aggregations do not change. E.g. 
#' for weekly data, the user would supply dt = 7L. If the aggregation windows 
#' vary in length, a vector of integers can be supplied. If the aggregations 
#' have a repeating pattern, for instance, if incidence is always reported three 
#' times per week on the same day of the week, you can supply a vector such as: 
#' dt = c(2L,2L,3L). If the aggregations change over time, or do not have a 
#' repeating pattern, the user can supply a full vector of aggregations matching 
#' the length of the incidence data supplied.}
#' \item{dt_out is the length of the sliding window used to estimate Rt from the 
#' reconstructed daily data. By default, Rt estimation uses weekly sliding time 
#' windows, however, we recommend that the sliding window is at least equal to 
#' the length of the longest aggregation window (dt) in the data.}
#' \item{iter is the number of iterations of the EM algorithm used to reconstruct
#' the daily incidence data. By default, iter = 10L, which has been demonstrated
#' to exceed the number of iterations necessary to reach convergence in simulation 
#' studies and analysis of real-world data by the package authors (manuscript in 
#' prep).}
#' \item{grid is a named list containing "precision", "min", and "max" which are 
#' used to define a grid of growth rate parameters used inside the EM algorithm.
#' The grid is used to convert reproduction number estimates for each aggregation
#' of incidence data into growth rates, which are then used to reconstruct the 
#' daily incidence data assuming exponential growth. The grid will auto-adjust 
#' if it is not large enough, so we recommend using the default values.}
#' }
#' 
#' There are three stages of the EM algorithm:
#' \itemize{
#' \item{Initialisation. The EM algorithm is initialised with a naive disaggregation 
#' of the incidence data. For example, if there were 70 cases over the course of a 
#' week, this would be naively split into 10 cases per day.}
#' \item{Expectation. The reproduction number is estimated for each aggregation 
#' window, except for the first aggregation window (as there is no past incidence 
#' data). This means that the earliest the incidence reconstruction can start is
#' at least the first day of the second aggregation window. Additionally, if the 
#' disaggregated incidence in subsequent aggregation windows is too low to estimate 
#' the reproduction number, this will mean that the reconstruction will not start 
#' until case numbers are sufficiently high.}
#' \item{Maximisation. The reproduction number estimates are then translated 
#' into growth rates for each aggregation window (Wallinga & Lipsitch, 2007) and 
#' used to reconstruct daily incidence data assuming exponential growth. The daily
#' incidence is adjusted by a constant to ensure that if the daily incidence were
#' to be re-aggregated, it would still sum to the original aggregated totals. The
#' expectation and maximisation steps repeat iteratively until convergence.}
#' }
#' 
#' The daily incidence that is reconstructed after the final iteration of the EM 
#' algorithm is then used to estimate Rt using the same process as the original 
#' \code{estimate_R} function, with sliding weekly time windows used as the default.
#' 
#' 
#' @seealso \itemize{
#'  \item{\code{\link{estimate_R}}}{ for details of the core function}
#'  }
#'  
#' @author Rebecca Nash \email{r.nash@imperial.ac.uk} and Anne Cori 
#' \email{a.cori@imperial.ac.uk}
#'  @references {
#' TODO: add Nash, R et al. 
#' Wallinga & Lipsitch. How generation intervals shape the relationship between 
#' growth rates and reproductive numbers (Proc Biol Sci 2007).
#' 
#' 
#' @examples 
#' ## Example for constant aggregation windows e.g. weekly reporting
#' 
#' # Load data on SARS in 2003
#' data("SARS2003")
#' 
#' # this is daily data, but for this example we will aggregate it to weekly 
#' # counts using the `aggregate_inc()` function
#' incid <- SARS2003$incidence
#' dt <- 7L
#' weekly_incid <- aggregate_inc(incid, dt)
#' si_distr <- SARS2003$si_distr
#' 
#' # estimate Rt using the default parameters (method "non_parametric_si")
#' method <- "non_parametric_si"
#' config <- make_config(list(si_distr = si_distr))
#' res_weekly <- estimate_R_agg(incid = weekly_incid, 
#'                             dt = 7L, # aggregation window of the data
#'                             dt_out = 7L, # desired sliding window length
#'                             iter = 10L,
#'                             config = config,
#'                             method = method,
#'                             grid = list(precision = 0.001, min = -1, max = 1))
#' 
#' # Plot the result
#' plot(res_weekly)
#' 
#' 
#' ## Example using repeating vector of aggregation windows e.g. for consistent 
#' ## reporting 3 times a week
#' 
#' # Using the SARS data again
#' data("SARS2003")
#' 
#' # For this example we will pretend data is being reported three times a week,
#' # in 2-day, 2-day, 3-day aggregations
#' incid <- SARS2003$incidence
#' dt <- c(2L,2L,3L)
#' agg_incid <- aggregate_inc(incid, dt)
#' si_distr <- SARS2003$si_distr 
#' 
#' # estimate Rt using the default parameters (method "non_parametric_si")
#' method <- "non_parametric_si"
#' config <- make_config(list(si_distr = si_distr))
#' res_agg <- estimate_R_agg(incid = agg_incid, 
#'                             dt = c(2L,2L,3L), # aggregation windows of the data
#'                             dt_out = 7L, # desired sliding window length
#'                             iter = 10L,
#'                             config = config,
#'                             method = method,
#'                             grid = list(precision = 0.001, min = -1, max = 1))
#' 
#' # Plot the result
#' plot(res_agg)
#'
#' ## Example using full vector of aggregation windows
#'
#' dt <- rep(c(2L,2L,3L), length.out=length(agg_incid))
#' si_distr <- SARS2003$si_distr 
#' 
#' # estimate Rt using the default parameters (method "non_parametric_si")
#' method <- "non_parametric_si"
#' config <- make_config(list(si_distr = si_distr))
#' res_agg <- estimate_R_agg(incid = agg_incid, 
#'                             dt = dt, # aggregation windows of the data
#'                             dt_out = 7L, # desired sliding window length
#'                             iter = 10L,
#'                             config = config,
#'                             method = method,
#'                             grid = list(precision = 0.001, min = -1, max = 1))
#' 
#' # Plot the result
#' plot(res_agg)
#'
estimate_R_agg <- function(incid,
                           dt = 7L, # aggregation window of the data
                           dt_out = 7L, # desired sliding window length
                           iter = 10L,
                           config = make_config(), 
                           method = c("non_parametric_si", "parametric_si"),
                           grid = list(precision = 0.001, min = -1, max = 1)){ 
  
  if (!is.integer(dt) | !is.integer(dt_out)) {
    stop ("dt and dt_out must be integers e.g. 7L")
  }
  if (!is.integer(iter)) {
    stop ("iter must be an integer e.g. 10L")
  }
  if (!is.list(grid) | !length(grid) == 3){
    stop ("grid must be a list of 3 elements: precision, min, and max")
  }
  if (grid$max < grid$min){
    stop ("grid max must be larger than grid min")
  }
  if (grid$precision > grid$max-grid$min){
    stop ("grid precision must be less than grid max - grid min")
  }
  if (!is.numeric(grid$precision) | !is.numeric(grid$min) | !is.numeric(grid$max)){
    stop ("grid precision, min, and max, must all be numeric")
  }
  if (!method == "parametric_si" && !method == "non_parametric_si"){
    stop ("'arg' should be one of 'non_parametric_si' and 'parametric_si'")
  }
  
  # Two configs:
  # 'config' for the R estimates used to reconstruct the incidence (internal to the
  # EM algorithm). These use a fixed window length matched to dt (aggregation window)
  # 'config_out' for the final estimated R using sliding windows (supplied by user)
  
  method <- match.arg(method) # potentially add an error message but maybe automatic
  config <- process_config(config)
  check_config(config, method)
  config_out <- config 
  
  # config$t_start and config$t_end used for the reconstruction. 
  # Rt estimation starts on 1st day of second dt. 
  # Width of fixed time windows match aggregations (dt):
  
  if (length(dt) == 1){
    T <- length(incid) * dt
    config$t_start <- seq(from = dt + 1, to = T - (dt - 1), dt)
    config$t_end <- seq(from = min(config$t_start) + (dt - 1),to = T, dt)
    } else if (length(dt) == length(incid)){
      T <- sum(dt)
      config$t_start <- cumsum(c(dt[1] + 1, dt[2:length(dt[-1])]))
      config$t_end <- cumsum(c(config$t_start[1] + dt[2] - 1, dt[3:length(dt)]))
      } else { # vector of repeating aggregations
        T <- sum(rep(dt, length.out = length(incid)))
        n_dt <- length(incid) # number of aggregations
        # reorder dt as R estimation starts on second aggregation window
        reo_dt_start <- c(dt[2:length(dt)], dt[1])
        reo_dt_end <- c(reo_dt_start[2:length(reo_dt_start)], reo_dt_start[1])
        config$t_start <- cumsum(c(dt[1] + 1, 
                               rep(reo_dt_start, length.out = n_dt - dt[1])))
        config$t_end <- cumsum(c(config$t_start[1] + reo_dt_start[1] - 1, 
                             rep(reo_dt_end, length.out = n_dt - dt[1])))
  }
  
  niter <- seq(1, iter, 1) 
  sim_inc <- matrix(NA, nrow = T, ncol = iter)
  
  for (i in 1:length(niter)){
    if (niter[i] == 1){
      # Initialisation of EM. Aggregated incidence split evenly:
      if (length(dt) == 1){
      dis <- incid / dt
      dis_inc <- rep(dis, each = dt)
      full_dt <- rep(dt, length(incid))
      } else if (length(dt) == length(incid)){
        dis <- incid / dt
        dis_inc <- rep(dis, times = dt)
        full_dt <- dt
        } else {
          full_dt <- rep(dt, length.out = length(incid))
          dis <- incid / full_dt
          dis_inc <- rep(dis, times = full_dt)
          }
      
      # Estimate R
      R <- estimate_R(dis_inc, 
                      method = method,
                      config = config)
      
      message("Estimated R for iteration: ", i)
      Mean_R <- R$R$`Mean(R)`
      
      if (any(is.na(Mean_R))){
        idx_na <- which(is.na(Mean_R))
        idx_reconstruct <- seq(min(R$R$t_start[-idx_na]), length(dis_inc))
        Mean_R <- Mean_R[!is.na(Mean_R)]
      } else {
        idx_reconstruct <- seq(min(R$R$t_start), length(dis_inc))
      }
      
      ## Incidence that can't be reconstructed (can't estimate R for
      ## first agg window or if incidence is too low)
      incid_not_to_reconstruct <- dis_inc[-idx_reconstruct]
      
      ## Incidence that can be reconstructed
      aggs_to_reconstruct <- seq(ceiling(idx_reconstruct[1] / dt[1]), length(incid))
      incid_to_reconstruct <- incid[aggs_to_reconstruct]
      
      # Translate R to growth rate
      get_r_from_R <- function(R, gt_mean, gt_sd, 
                               gt_distr,
                               grid) {
        r_grid <- seq(grid$min, grid$max, grid$precision)
        if (is.null(gt_distr)) {
          gt_pars <- gamma_mucv2shapescale(mu = gt_mean, cv = gt_sd / gt_mean)
          gt_distr <- distcrete("gamma", interval = 1,
                                           shape = gt_pars$shape,
                                           scale = gt_pars$scale, w = 0.5)
        }
        # using a grid of r values translate that into R using r2R0
        R_grid <- epitrix::r2R0(r = r_grid, w = gt_distr)
        
        # find location of the value in the R grid which has the smallest 
        # difference to the input of R the user provided e.g. R_grid[idx_r]:
        idx_r <- vapply(R, function(e) which.min(abs(R_grid - e)), numeric(1L)) 
        
        while (any(idx_r == 1) || any(idx_r == length(r_grid))) {
          # if necessary rerun get_r_from_R with a wider r_grid
          grid_multiplier <- 5
          if (grid$max > 0) grid$max <- grid_multiplier * grid$max else grid$max <- - grid$max
          if (grid$min < 0) grid$min <- grid_multiplier * grid$min else grid$min <- - grid$min
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
      
      # Estimate incidence using growth rate
      
      # Assume that It is a constant (k) multiplied by exp(gr[for that dt]*t)
      d <- numeric(length(aggs_to_reconstruct))
      k <- numeric(length(aggs_to_reconstruct))
      dt_seq <- full_dt[aggs_to_reconstruct]
      
      for (w in 1:length(k)){
        d[w] <- sum(exp(gr[w] * seq(1, dt_seq[w] - 1, 1)))
        k[w] <- incid_to_reconstruct[w] / (exp(gr[w]) * (1 + d[w]))
      }
      
      recon_df <- data.frame(k = k, gr = gr, dt = dt_seq)
      k_ls <- list()
      gr_ls <- list()
      
      for (f in 1:length(incid_to_reconstruct)){
        k_ls[[f]] <- rep(recon_df$k[f], recon_df$dt[f])
        gr_ls[[f]] <- rep(recon_df$gr[f], recon_df$dt[f])
      }
      k_seq <- unlist(k_ls)
      gr_seq <- unlist(gr_ls)
      
      recon_days <- seq(1, sum(dt_seq))
      w_day <- list()

      for (x in 1:length(dt_seq)){
        w_day[[x]] <- seq(1,dt_seq[x])
      }

      w_day <- unlist(w_day)

      # Reconstruct daily incidence
      est_inc <- rep(NA, length(recon_days))
      
      for (t in seq_along(recon_days)){
        est_inc[t] <- k_seq[t] * exp(gr_seq[t] * w_day[t])
      }
      
      ## For the incidence that can't be reconstructed (can't estimate R for
      ## first agg window or if incidence is too low), using the initial dis_inc
      sim_inc[,i] <- c(incid_not_to_reconstruct, est_inc)
      
      message("Reconstructed incidence for iteration: ", i)
      
      
    } else {
      
      # Use new adjusted incidence as starting point
      new_inc <- sim_inc[,i-1]
      
      # Re-Estimate R
      R <- estimate_R(new_inc,
                      method = method,
                      config = config)
      
      message("Estimated R for iteration: ", i)
      
      Mean_R <- R$R$`Mean(R)`
      
      if (any(is.na(Mean_R))){
        idx_na <- which(is.na(Mean_R))
        idx_reconstruct <- seq(min(R$R$t_start[-idx_na]), length(dis_inc))
        Mean_R <- Mean_R[!is.na(Mean_R)]
      } else {
        idx_reconstruct <- seq(min(R$R$t_start), length(dis_inc))
      }
      
      incid_not_to_reconstruct <- dis_inc[-idx_reconstruct]
      
      aggs_to_reconstruct <- seq(ceiling(idx_reconstruct[1] / dt[1]), length(incid))
      incid_to_reconstruct <- incid[aggs_to_reconstruct]
      
      # Translate R to growth rate again
      gr <- get_r_from_R(R = Mean_R, 
                         gt_mean = config$mean_si, gt_sd = config$std_si, 
                         gt_distr = config$si_distr,
                         grid = grid)
      
      # Estimate incidence
      d <- numeric(length(aggs_to_reconstruct))
      k <- numeric(length(aggs_to_reconstruct))
      dt_seq <- full_dt[aggs_to_reconstruct]
      
      for (w in 1:length(k)){
        d[w] <- sum(exp(gr[w] * seq(1, dt_seq[w] - 1, 1)))
        k[w] <- incid_to_reconstruct[w] / (exp(gr[w]) * (1 + d[w]))
      }
      
      recon_df <- data.frame(k = k, gr = gr, dt = dt_seq)
      k_ls <- list()
      gr_ls <- list()
      
      for (f in 1:length(incid_to_reconstruct)){
        k_ls[[f]] <- rep(recon_df$k[f], recon_df$dt[f])
        gr_ls[[f]] <- rep(recon_df$gr[f], recon_df$dt[f])
      }
      k_seq <- unlist(k_ls)
      gr_seq <- unlist(gr_ls)
      
      w_day <- list()
      for (x in 1:length(dt_seq)){
        w_day[[x]] <- seq(1 , dt_seq[x])
      }
      w_day <- unlist(w_day)

      recon_days <- seq(1, sum(dt_seq))    
      est_inc <- rep(NA, length(recon_days))
      for (t in seq_along(recon_days)){
        est_inc[t] <- k_seq[t] * exp(gr_seq[t] * w_day[t])
        }
        
      sim_inc[,i] <- c(incid_not_to_reconstruct, est_inc)
      
      # monitor progress:
      message("Reconstructed incidence for iteration: ", i)
      
      # Final estimate R starting on the first aggregation window 
      # that incidence was able to be reconstructed over
      
      if (is.null(config_out$t_start)) {
        config_out$t_start <- seq(from = min(R$R$t_start), 
                                  to = T - (dt_out - 1), 1)
      }
      if (is.null(config_out$t_end)) {
        config_out$t_end <- config_out$t_start + (dt_out - 1)
      }
      
      if (niter[i] == max(niter)){
        R_out <- estimate_R(sim_inc[,i],
                            method = method,
                            config = config_out)
        message("R estimation starts on day ", R_out$R$t_start[1])
      }
      
    }
  }
  
  R_out
  
}
