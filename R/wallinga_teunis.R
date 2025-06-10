##########################################################################
## wallinga_teunis function to estimate Rc the case reproduction number ##
##########################################################################

#' Estimation of the case reproduction number using the Wallinga and Teunis
#' method
#'
#' \code{wallinga_teunis} estimates the case reproduction number of an epidemic,
#' given the incidence time series and the serial interval distribution.
#'
#' @param incid One of the following
#'   * Vector (or a dataframe with
#'   a column named 'incid') of non-negative integers containing an incidence
#'   time series. If the dataframe contains a column \code{incid$dates}, this is
#'   used for plotting. \code{incid$dates} must contains only dates in a row.
#'
#'   * An object of class [incidence]

#'
#' @param method the method used to estimate R, one of "non_parametric_si",
#'   "parametric_si", "uncertain_si", "si_from_data" or "si_from_sample"
#'
#' @param config a list with the following elements: 
#'   * t_start:
#'   Vector of positive integers giving the starting times of each window over
#'   which the reproduction number will be estimated. These must be in ascending
#'   order, and so that for all \code{i}, \code{t_start[i]<=t_end[i]}.
#'   t_start[1] should be strictly after the first day with non null incidence.
#'   * t_end: Vector of positive integers giving the ending times of each
#'   window over which the reproduction number will be estimated. These must be
#'   in ascending order, and so that for all \code{i},
#'   \code{t_start[i]<=t_end[i]}.
#'   * method: One of "non_parametric_si" or
#'   "parametric_si" (see details).
#'   * mean_si: For method "parametric_si" ;
#'   positive real giving the mean serial interval.
#'   * std_si: For method
#'   "parametric_si" ; non negative real giving the standard deviation of the
#'   serial interval.
#'   * si_distr: For method "non_parametric_si" ; vector
#'   of probabilities giving the discrete distribution of the serial interval,
#'   starting with \code{si_distr[1]} (probability that the serial interval is
#'   zero), which should be zero.
#'   * n_sim: A positive integer giving the
#'   number of simulated epidemic trees used for computation of the confidence
#'   intervals of the case reproduction number (see details).
#' @return { a list with components: 
#'   * R: a dataframe
#'   containing: the times of start and end of each time window considered ; the
#'   estimated mean, std, and 0.025 and 0.975 quantiles of the reproduction
#'   number for each time window.
#'   * si_distr: a vector containing the
#'   discrete serial interval distribution used for estimation
#'   * SI.Moments: a vector containing the mean and std of the discrete
#'   serial interval distribution(s) used for estimation
#'   * I: the time
#'   series of total incidence 
#'   * I_local: the time series of incidence of
#'   local cases (so that \code{I_local + I_imported = I})
#'   * I_imported:
#'   the time series of incidence of imported cases (so that \code{I_local +
#'   I_imported = I}) 
#'   * dates: a vector of dates corresponding to the
#'   incidence time series }
#'
#' @details Estimates of the case reproduction number for an epidemic over
#' predefined time windows can be obtained, for a given discrete distribution of
#' the serial interval, as proposed by Wallinga and Teunis (AJE, 2004).
#' Confidence intervals are obtained by simulating a number (config$n_sim) of
#' possible transmission trees (only done if config$n_sim > 0).
#'
#' The methods vary in the way the serial interval distribution is specified.
#'
#' ----------------------- \code{method "non_parametric_si"}
#' -----------------------
#'
#' The discrete distribution of the serial interval is directly specified in the
#' argument \code{config$si_distr}.
#'
#'
#' ----------------------- \code{method "parametric_si"} -----------------------
#'
#' The mean and standard deviation of the continuous distribution of the serial
#' interval are given in the arguments \code{config$mean_si} and
#' \code{config$std_si}. The discrete distribution of the serial interval is
#' derived automatically using \code{\link{discr_si}}.
#'
#'
#' @seealso \code{\link{discr_si}}, \code{\link{estimate_R}}
#' @author Anne Cori \email{a.cori@imperial.ac.uk}
#' @references { Cori, A. et al. A new framework and software to estimate
#'   time-varying reproduction numbers during epidemics (AJE 2013). Wallinga, J.
#'   and P. Teunis. Different epidemic curves for severe acute respiratory
#'   syndrome reveal similar impacts of control measures (AJE 2004). }
#' @export
#' @import reshape2 gridExtra
#' @importFrom ggplot2 last_plot ggplot aes geom_step ggtitle geom_ribbon
#'   geom_line xlab ylab xlim geom_hline ylim geom_histogram
#' @examples
#' ## load data on pandemic flu in a school in 2009
#' data("Flu2009")
#'
#' ## estimate the case reproduction number (method "non_parametric_si")
#' res <- wallinga_teunis(Flu2009$incidence,
#'    method="non_parametric_si",
#'    config = list(t_start = seq(2, 26), t_end = seq(8, 32),
#'                  si_distr = Flu2009$si_distr,
#'                  n_sim = 100))
#' plot(res)
#' ## the second plot produced shows, at each each day,
#' ## the estimate of the case reproduction number over the 7-day window
#' ## finishing on that day.
#'
#' ## estimate the case reproduction number (method "parametric_si")
#' res <- wallinga_teunis(Flu2009$incidence, method="parametric_si",
#'    config = list(t_start = seq(2, 26), t_end = seq(8, 32),
#'                  mean_si = 2.6, std_si = 1.5,
#'                  n_sim = 100))
#' plot(res)
#' ## the second plot produced shows, at each each day,
#' ## the estimate of the case reproduction number over the 7-day window
#' ## finishing on that day.
wallinga_teunis <- function(incid,
                            method = c("non_parametric_si", "parametric_si"),
                            config) {
  
  ### Functions ###
  
  #########################################################
  # Draws a possile transmission tree                     #
  #########################################################
  
  draw_one_set_of_ancestries <- function() {
    res <- vector()
    for (t in seq_len(T))
    {
      if (length(which(Onset == t)) > 0) {
        if (length(possible_ances_time[[t]]) > 0) {
          prob <- config$si_distr[t - possible_ances_time[[t]] + 1] *
            incid[possible_ances_time[[t]]]
          ot <- which(Onset == t)
          if (any(prob > 0)) {
            res[ot] <-
              possible_ances_time[[t]][which(rmultinom(length(ot),
                                                       size = 1, prob = prob)
                                             == TRUE, arr.ind = TRUE)[, 1]]
          } else {
            res[ot] <- NA
          }
        } else {
          res[which(Onset == t)] <- NA
        }
      }
    }
    return(res)
  }
  
  ### Error messages ###
  
  method <- match.arg(method)
  
  incid <- process_I(incid)
  if (!is.null(incid$dates)) {
    dates <- check_dates(incid)
    incid <- process_I_vector(rowSums(incid[, c("local", "imported")]))
    T <- length(incid)
  } else {
    incid <- process_I_vector(rowSums(incid[, c("local", "imported")]))
    T <- length(incid)
    dates <- seq_len(T)
  }
  
  
  ### Adjusting t_start and t_end so that at least an incident case has been
  ### observed before t_start[1] ###
  
  i <- 1
  while (sum(incid[seq_len(config$t_start[i] - 1)]) == 0) {
    i <- i + 1
  }
  temp <- which(config$t_start < i)
  if (length(temp > 0)) {
    config$t_start <- config$t_start[-temp]
    config$t_end <- config$t_end[-temp]
  }
  
  check_times(config$t_start, config$t_end, T)
  nb_time_periods <- length(config$t_start)
  
  if (is.null(config$n_sim)) {
    config$n_sim <- 10
    warning("setting config$n_sim to 10 as config$n_sim was not specified. ")
  }
  
  if (method == "non_parametric_si") {
    check_si_distr(config$si_distr)
    config$si_distr <- c(config$si_distr, 0)
  }
  
  if (method == "parametric_si") {
    if (is.null(config$mean_si)) {
      stop("method non_parametric_si requires to specify the config$mean_si
           argument.")
    }
    if (is.null(config$std_si)) {
      stop("method non_parametric_si requires to specify the config$std_si
           argument.")
    }
    if (config$mean_si < 1) {
      stop("method parametric_si requires a value >1 for config$mean_si.")
    }
    if (config$std_si < 0) {
      stop("method parametric_si requires a >0 value for config$std_si.")
    }
  }
  
  if (!is.numeric(config$n_sim)) {
    stop("config$n_sim must be a positive integer.")
  }
  if (config$n_sim < 0) {
    stop("config$n_sim must be a positive integer.")
  }
  
  ### What does each method do ###
  
  if (method == "non_parametric_si") {
    parametric_si <- "N"
  }
  if (method == "parametric_si") {
    parametric_si <- "Y"
  }
  
  if (parametric_si == "Y") {
    config$si_distr <- discr_si(seq(0,T - 1), config$mean_si, config$std_si)
  }
  if (length(config$si_distr) < T + 1) {
    config$si_distr[seq(length(config$si_distr) + 1, T + 1)] <- 0
  }
  
  final_mean_si <- sum(config$si_distr * (seq(0, length(config$si_distr) - 1)))
  final_std_si <- sqrt(sum(config$si_distr * 
                             (seq(0, length(config$si_distr) - 1))^2) - 
                         final_mean_si^2)
  
  time_periods_with_no_incidence <- vector()
  for (i in seq_len(nb_time_periods))
  {
    if (sum(incid[seq(config$t_start[i],config$t_end[i])]) == 0) {
      time_periods_with_no_incidence <- c(time_periods_with_no_incidence, i)
    }
  }
  if (length(time_periods_with_no_incidence) > 0) {
    config$t_start <- config$t_start[-time_periods_with_no_incidence]
    config$t_end <- config$t_end[-time_periods_with_no_incidence]
    nb_time_periods <- length(config$t_start)
  }
  
  Onset <- vector()
  for (t in seq_len(T)) {
    Onset <- c(Onset, rep(t, incid[t]))
  }
  NbCases <- length(Onset)
  
  delay <- outer(seq_len(T), seq_len(T), "-")
  si_delay <- apply(delay, 2, function(x) 
    config$si_distr[pmin(pmax(x + 1, 1), length(config$si_distr))])
  sum_on_col_si_delay_tmp <- vnapply(seq_len(nrow(si_delay)), function(i) 
    sum(si_delay [i, ] * incid, na.rm = TRUE))
  sum_on_col_si_delay <- vector()
  for (t in seq_len(T)) {
    sum_on_col_si_delay <- c(sum_on_col_si_delay, 
                             rep(sum_on_col_si_delay_tmp[t], incid[t]))
  }
  mat_sum_on_col_si_delay <- matrix(rep(sum_on_col_si_delay_tmp, T),
                                    nrow = T, ncol = T)
  p <- si_delay / (mat_sum_on_col_si_delay)
  p[which(is.na(p))] <- 0
  p[which(is.infinite(p))] <- 0
  
  mean_r_per_index_case_date <- vnapply(seq_len(ncol(p)), function(j) 
    sum(p[, j] * incid, na.rm = TRUE))
  mean_r_per_date_wt <- vnapply(seq_len(nb_time_periods), function(i) 
    mean(rep(mean_r_per_index_case_date[which((seq_len(T) >= 
                                                 config$t_start[i]) * 
                                                (seq_len(T) <= 
                                                   config$t_end[i]) == 1)],
             incid[which((seq_len(T) >= config$t_start[i]) * 
                           (seq_len(T) <= config$t_end[i]) == 1)])))
  
  if(config$n_sim>0)
  {
    possible_ances_time <- lapply(seq_len(T), function(t) 
      (t - (which(config$si_distr != 0)) + 
         1)[which(t - (which(config$si_distr != 0)) + 1 > 0)])
    
    
    ancestries_time <- t(vapply(seq_len(config$n_sim), function(i) 
      draw_one_set_of_ancestries(), numeric(sum(incid))))
    
    r_sim <- vapply(seq_len(nb_time_periods), function(i) 
      rowSums((ancestries_time[, ] >= config$t_start[i]) * 
                (ancestries_time[, ] <= config$t_end[i]), na.rm = TRUE) / 
        sum(incid[seq(config$t_start[i], config$t_end[i])]), 
      numeric(config$n_sim))
    
    r025_wt <- apply(r_sim, 2, quantile, 0.025, na.rm = TRUE)
    r025_wt <- r025_wt[which(!is.na(r025_wt))]
    r975_wt <- apply(r_sim, 2, quantile, 0.975, na.rm = TRUE)
    r975_wt <- r975_wt[which(!is.na(r975_wt))]
    std_wt <- apply(r_sim, 2, sd, na.rm = TRUE)
    std_wt <- std_wt[which(!is.na(std_wt))]
  }else
  {
    r025_wt <- rep(NA, length(mean_r_per_date_wt))
    r975_wt <- rep(NA, length(mean_r_per_date_wt))
    std_wt <- rep(NA, length(mean_r_per_date_wt))
  }
  
  results <- list(R = as.data.frame(cbind(config$t_start,
                                          config$t_end, mean_r_per_date_wt,
                                          std_wt, r025_wt, r975_wt)))
  
  names(results$R) <- c("t_start", "t_end", "Mean(R)", "Std(R)",
                        "Quantile.0.025(R)", "Quantile.0.975(R)")
  
  results$method <- method
  results$si_distr <- config$si_distr
  
  results$SI.Moments <- as.data.frame(cbind(final_mean_si, final_std_si))
  names(results$SI.Moments) <- c("Mean", "Std")
  
  if (!is.null(dates)) {
    results$dates <- dates
  }
  results$I <- incid
  results$I_local <- incid
  results$I_local[1] <- 0
  results$I_imported <- c(incid[1], rep(0, length(incid) - 1))
  
  class(results) <- "estimate_R"
  return(results)
}
