#' Estimated Instantaneous Reproduction Number
#'
#' \code{estimate_R} estimates the reproduction number of an epidemic, given the
#' incidence time series and the serial interval distribution.
#'
#' @param incid One of the following
#' \itemize{
#'
#' \item{A vector (or a dataframe with a single column) of non-negative integers
#' containing the incidence time series; these can be aggregated at any time
#' unit as specified by argument \code{dt}}
#'
#' \item{A dataframe of non-negative integers with either i) \code{incid$I}
#' containing the total incidence, or ii) two columns, so that
#' \code{incid$local} contains the incidence of cases due to local transmission
#' and \code{incid$imported} contains the incidence of imported cases (with
#' \code{incid$local + incid$imported} the total incidence). If the dataframe
#' contains a column \code{incid$dates}, this is used for plotting.
#' \code{incid$dates} must contains only dates in a row.}
#'
#' \item{An object of class \code{\link{incidence}}}
#'
#' }
#'
#' Note that the cases from the first time step are always all assumed to be
#' imported cases.
#'
#' @param method One of "non_parametric_si", "parametric_si", "uncertain_si",
#'   "si_from_data" or "si_from_sample" (see details).
#'
#' @param si_sample For method "si_from_sample" ; a matrix where each column
#'   gives one distribution of the serial interval to be explored (see details).
#'
#' @param si_data For method "si_from_data" ; the data on dates of symptoms of
#'   pairs of infector/infected individuals to be used to estimate the serial
#'   interval distribution should be a dataframe with 5 columns:
#'   \itemize{
#'   \item{EL: the lower bound
#'   of the symptom onset date of the infector (given as an integer)}
#'   \item{ER:
#'   the upper bound of the symptom onset date of the infector (given as an
#'   integer). Should be such that ER>=EL. If the dates are known exactly use
#'   ER = EL}
#'   \item{SL: the lower bound of the
#'   symptom onset date of the infected individual (given as an integer)}
#'   \item{SR: the upper bound of the symptom onset date of the infected
#'   individual (given as an integer). Should be such that SR>=SL. If the dates
#'   are known exactly use SR = SL}
#'   \item{type
#'   (optional): can have entries 0, 1, or 2, corresponding to doubly
#'   interval-censored, single interval-censored or exact observations,
#'   respectively, see Reich et al. Statist. Med. 2009. If not specified, this
#'   will be automatically computed from the dates}
#'   }
#'
#' @param config An object of class \code{estimate_R_config}, as returned by
#' function \code{make_config}.
#'
#' @param dt length of temporal aggregations of the incidence data. This should
#' be an integer or vector of integers. If a vector, this can either match the
#' length of the incidence data supplied, or it will be recycled. For
#' example, \code{dt = c(3L, 4L)} would correspond to alternating incidence
#' aggregation windows of 3 and 4 days. The default value is 1 time unit
#' (typically day).
#'
#' @param dt_out length of the sliding windows used for R estimates (integer,
#' 7 time units (typically days) by default).
#' Only used if \code{dt > 1};
#' in this case this will superseed config$t_start and config$t_end,
#' see \code{\link{estimate_R_agg}}.
#'
#' @param recon_opt one of "naive" or "match", see \code{\link{estimate_R_agg}}.
#'
#' @param iter number of iterations of the EM algorithm used to reconstruct
#' incidence at 1-time-unit intervals(integer, 10 by default).
#' Only used if \code{dt > 1}, see \code{\link{estimate_R_agg}}.
#'
#' @param tol tolerance used in the convergence check (numeric, 1e-6 by default),
#' see \code{\link{estimate_R_agg}}.
#'
#' @param grid named list containing "precision", "min", and "max" which are
#' used to define a grid of growth rate parameters that are used inside the EM
#' algorithm used to reconstruct incidence at 1-time-unit intervals.
#' Only used if \code{dt > 1}, see \code{\link{estimate_R_agg}}.
#'
#' @param backimputation_window Length of the window used to impute incidence
#'    prior to the first reported cases. The default value is 0, meaning that no
#'    back-imputation is performed. If a positive integer is provided, the
#'    incidence is imputed for the first \code{backimputation_window} time 
#'    units.
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
#' \item{I_imputed}{: the time series of incidence of imputed cases}

#' \item{dates}{: a vector of dates corresponding to the incidence time series}
#'
#' \item{MCMC_converged}{ (only for method \code{si_from_data}): a boolean
#' showing whether the Gelman-Rubin MCMC convergence diagnostic was successful
#' (\code{TRUE}) or not (\code{FALSE})}
#' }
#' }
#'
#' @details
#' Analytical estimates of the reproduction number for an epidemic over
#' predefined time windows can be obtained within a Bayesian framework,
#' for a given discrete distribution of the serial interval (see references).
#'
#' Several methods are available to specify the serial interval distribution.
#'
#' In short there are five methods to specify the serial interval distribution
#' (see help for function \code{make_config} for more detail on each method).
#' In the first two methods, a unique serial interval distribution is
#' considered, whereas in the last three, a range of serial interval
#' distributions are integrated over:
#' \itemize{
#' \item{In method "non_parametric_si" the user specifies the discrete
#' distribution of the serial interval}
#' \item{In method "parametric_si" the user specifies the mean and sd of the
#' serial interval}
#' \item{In method "uncertain_si" the mean and sd of the serial interval are
#' each drawn from truncated normal distributions, with parameters specified by
#' the user}
#' \item{In method "si_from_data", the serial interval distribution is directly
#' estimated, using MCMC, from interval censored exposure data, with data
#' provided by the user together with a choice of parametric distribution for
#' the serial interval}
#' \item{In method "si_from_sample", the user directly provides the sample of
#' serial interval distribution to use for estimation of R. This can be a useful
#'  alternative to the previous method, where the MCMC estimation of the serial
#'  interval distribution could be run once, and the same estimated SI
#'  distribution then used in estimate_R in different contexts, e.g. with
#'  different time windows, hence avoiding to rerun the MCMC every time
#'  estimate_R is called.}
#' }
#'
#' R is estimated within a Bayesian framework, using a Gamma distributed prior,
#' with mean and standard deviation which can be set using the `mean_prior`
#' and `std_prior` arguments within the `make_config` function, which can then
#' be used to specify `config` in the `estimate_R` function.
#' Default values are a mean prior of 5 and standard deviation of 5.
#' This was set to a high prior value with large uncertainty so that if one
#' estimates R to be below 1, the result is strongly data-driven.
#'
#' R is estimated on time windows specified through the `config` argument.
#' These can be overlapping or not (see `make_config` function and vignette
#' for examples).
#'
#' @seealso \itemize{
#'  \item{\code{\link{make_config}}}{ for general settings of the estimation}
#'  \item{\code{\link{discr_si}}}{ to build serial interval distributions}
#'  \item{\code{\link{sample_posterior_R}}}{ to draw samples of R values from
#' the posterior distribution from the output of \code{estimate_R()}
#' }}
#'
#'
#' @author Anne Cori \email{a.cori@imperial.ac.uk}
#' @references {
#' Cori, A. et al. A new framework and software to estimate time-varying
#' reproduction numbers during epidemics (AJE 2013).
#' Wallinga, J. and P. Teunis. Different epidemic curves for severe acute
#' respiratory syndrome reveal similar impacts of control measures (AJE 2004).
#' Reich, N.G. et al. Estimating incubation period distributions with coarse
#' data (Statis. Med. 2009)
#' }
#' @importFrom coarseDataTools dic.fit.mcmc
#' @importFrom coda as.mcmc.list as.mcmc
#' @importFrom incidence incidence
#' @export
#' @examples
#' ## load data on pandemic flu in a school in 2009
#' data("Flu2009")
#'
#' ## estimate the reproduction number (method "non_parametric_si")
#' ## when not specifying t_start and t_end in config, they are set to estimate
#' ## the reproduction number on sliding weekly windows                      
#' res <- estimate_R(incid = Flu2009$incidence,
#'                   method = "non_parametric_si",
#'                   config = make_config(list(si_distr = Flu2009$si_distr)))
#' plot(res)
#'
#' ## the second plot produced shows, at each each day,
#' ## the estimate of the reproduction number over the 7-day window
#' ## finishing on that day.
#'
#' ## to specify t_start and t_end in config, e.g. to have biweekly sliding
#' ## windows  
#' t_start <- seq(2, nrow(Flu2009$incidence)-13)
#' t_end <- t_start + 13             
#' res <- estimate_R(incid = Flu2009$incidence,
#'                   method = "non_parametric_si",
#'                   config = make_config(list(
#'                       si_distr = Flu2009$si_distr,
#'                       t_start = t_start,
#'                       t_end = t_end)))
#' plot(res)
#'
#' ## the second plot produced shows, at each each day,
#' ## the estimate of the reproduction number over the 14-day window
#' ## finishing on that day.
#'
#' ## example with an incidence object
#'
#' ## create fake data
#' library(incidence)
#' data <- c(0,1,1,2,1,3,4,5,5,5,5,4,4,26,6,7,9)
#' location <- sample(c("local","imported"), length(data), replace=TRUE)
#' location[1] <- "imported" # forcing the first case to be imported
#'
#' ## get incidence per group (location)
#' incid <- incidence(data, groups = location)
#'
#' ## Estimate R with assumptions on serial interval
#' res <- estimate_R(incid, method = "parametric_si",
#'                   config = make_config(list(
#'                   mean_si = 2.6, std_si = 1.5)))
#' plot(res)
#' ## the second plot produced shows, at each each day,
#' ## the estimate of the reproduction number over the 7-day window
#' ## finishing on that day.
#'
#' ## estimate the reproduction number (method "parametric_si")
#' res <- estimate_R(Flu2009$incidence, method = "parametric_si",
#'                   config = make_config(list(mean_si = 2.6, std_si = 1.5)))
#' plot(res)
#' ## the second plot produced shows, at each each day,
#' ## the estimate of the reproduction number over the 7-day window
#' ## finishing on that day.
#'
#' ## estimate the reproduction number (method "uncertain_si")
#' res <- estimate_R(Flu2009$incidence, method = "uncertain_si",
#'                   config = make_config(list(
#'                   mean_si = 2.6, std_mean_si = 1,
#'                   min_mean_si = 1, max_mean_si = 4.2,
#'                   std_si = 1.5, std_std_si = 0.5,
#'                   min_std_si = 0.5, max_std_si = 2.5,
#'                   n1 = 100, n2 = 100)))
#' plot(res)
#' ## the bottom left plot produced shows, at each each day,
#' ## the estimate of the reproduction number over the 7-day window
#' ## finishing on that day.
#'
#' ## Example with back-imputation:
#' ## here we use the first 6 days of incidence to impute cases that preceded
#' ## the first reported cases:
#' 
#' res_bi <- estimate_R(incid = Flu2009$incidence,
#'                  method = "parametric_si",
#'                  backimputation_window = 6,
#'                  config = make_config(list(
#'                       mean_si = 2.6, 
#'                       std_si = 1,
#'                       t_start = t_start,
#'                       t_end = t_end)))
#' plot(res_bi, "R")
#' 
#' ## We can see that early estimates of R are lower when back-imputation is 
#' ## used, even though the difference is marginal in this case.
#'
#' \dontrun{
#' ## Note the following examples use an MCMC routine
#' ## to estimate the serial interval distribution from data,
#' ## so they may take a few minutes to run
#'
#' ## load data on rotavirus
#' data("MockRotavirus")
#'
#'################
#' mcmc_control <- make_mcmc_control(
#'   burnin = 1000, # first 1000 iterations discarded as burn-in
#'   thin = 10, # every 10th iteration will be kept, the rest discarded
#'   seed = 1) # set the seed to make the process reproducible
#'
#' R_si_from_data <- estimate_R(
#'   incid = MockRotavirus$incidence,
#'   method = "si_from_data",
#'   si_data = MockRotavirus$si_data, # symptom onset data
#'   config = make_config(si_parametric_distr = "G", # gamma dist. for SI
#'      mcmc_control = mcmc_control,
#'      n1 = 500, # number of posterior samples of SI dist.
#'      n2 = 50, # number of posterior samples of Rt dist.
#'      seed = 2)) # set seed for reproducibility
#'
#' ## compare with version with no uncertainty
#' R_Parametric <- estimate_R(MockRotavirus$incidence,
#'                           method = "parametric_si",
#'                           config = make_config(list(
#'                           mean_si = mean(R_si_from_data$SI.Moments$Mean),
#'                              std_si = mean(R_si_from_data$SI.Moments$Std))))
#' ## generate plots
#' p_uncertainty <- plot(R_si_from_data, "R", options_R=list(ylim=c(0, 1.5)))
#' p_no_uncertainty <- plot(R_Parametric, "R", options_R=list(ylim=c(0, 1.5)))
#' gridExtra::grid.arrange(p_uncertainty, p_no_uncertainty,ncol=2)
#'
#' ## the left hand side graph is with uncertainty in the SI distribution, the
#' ## right hand side without.
#' ## The credible intervals are wider when accounting for uncertainty in the SI
#' ## distribution.
#'
#' ## estimate the reproduction number (method "si_from_sample")
#' MCMC_seed <- 1
#' overall_seed <- 2
#' SI.fit <- coarseDataTools::dic.fit.mcmc(dat = MockRotavirus$si_data,
#'                  dist = "G",
#'                  init.pars = init_mcmc_params(MockRotavirus$si_data, "G"),
#'                  burnin = 1000,
#'                  n.samples = 5000,
#'                  seed = MCMC_seed)
#' si_sample <- coarse2estim(SI.fit, thin = 10)$si_sample
#' R_si_from_sample <- estimate_R(MockRotavirus$incidence,
#'                                method = "si_from_sample",
#'                                si_sample = si_sample,
#'                                config = make_config(list(n2 = 50,
#'                                seed = overall_seed)))
#' plot(R_si_from_sample)
#'
#' ## check that R_si_from_sample is the same as R_si_from_data
#' ## since they were generated using the same MCMC algorithm to generate the SI
#' ## sample (either internally to EpiEstim or externally)
#' all(R_si_from_sample$R$`Mean(R)` == R_si_from_data$R$`Mean(R)`)
#' }
#'
estimate_R <- function(incid,
                       method = c(
                         "non_parametric_si", "parametric_si",
                         "uncertain_si", "si_from_data",
                         "si_from_sample"
                       ),
                       si_data = NULL,
                       si_sample = NULL,
                       config = make_config(incid = incid, method = method),
                       dt = 1L, # aggregation window of the data
                       dt_out = 7L, # desired sliding window length
                       recon_opt = "naive",
                       iter = 10L,
                       tol = 1e-6,
                       grid = list(precision = 0.001, min = -1, max = 1),
                       backimputation_window = 0
                       ) {

  method <- match.arg(method)

  # switch between the standard estimate_R version and that which disaggregates
  # coarsely aggregated incidence data
  if (any(dt >= 2)) {
    
    msg <- "backimputation_window is currently not supported when dt > 1"
    if(backimputation_window > 0) stop(msg)
    
    out <- estimate_R_agg(incid, dt = dt, dt_out = dt_out, iter = iter,
                          config = config, method = method, grid = grid)
    return(out)
  }

  # Impute missed generations of infections if required
  if (backimputation_window) {
    incid <- backimpute_I(incid, window_b = backimputation_window)
  }

  config <- make_config(incid = incid, method = method, config = config)
  config <- process_config(config)
  check_config(config, method)


  # If the serial interval distribution is not provided, estimate it
  if (method == "si_from_data") {
    ## Warning if the expected set of parameters is not adequate
    si_data <- process_si_data(si_data)
    config <- process_config_si_from_data(config, si_data)

    ## estimate serial interval from serial interval data first
    if (!is.null(config$mcmc_control$seed)) {
      cdt <- dic.fit.mcmc(
        dat = si_data,
        dist = config$si_parametric_distr,
        burnin = config$mcmc_control$burnin,
        n.samples = config$n1 * config$mcmc_control$thin,
        init.pars = config$mcmc_control$init_pars,
        seed = config$mcmc_control$seed
      )
    } else {
      cdt <- dic.fit.mcmc(
        dat = si_data,
        dist = config$si_parametric_distr,
        burnin = config$mcmc_control$burnin,
        n.samples = config$n1 * config$mcmc_control$thin,
        init.pars = config$mcmc_control$init_pars
      )
    }

    ## check convergence of the MCMC and print warning if not converged
    MCMC_conv <- check_cdt_samples_convergence(cdt@samples)

    ## thin the chain, and turn the two parameters of the SI distribution into a
    ## whole discrete distribution
    c2e <- coarse2estim(cdt, thin = config$mcmc_control$thin)

    cat(paste(
      "\n\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@",
      "\nEstimating the reproduction number for these serial interval",
      "estimates...\n",
      "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
    ))
    si_sample <- c2e$si_sample
  }

  ## estimate R whether or not si_sample is simulated
  if (!is.null(config$seed)) {
    set.seed(config$seed)
  }
  out <- estimate_R_func(
    incid = incid,
    method = method, 
    si_sample = si_sample,
    config = config
  )
  
  # Add extra fields based on method
  if(method == "si_from_data"){
    out[["MCMC_converged"]] <- MCMC_conv
  }

  return(out)
}

##########################################################
## estimate_R_func: Doing the heavy work in estimate_R  ##
##########################################################

#'
#' @importFrom stats median qgamma quantile rnorm sd
#'
#' @importFrom incidence as.incidence
#'
estimate_R_func <- function(incid,
                            si_sample,
                            method = c(
                              "non_parametric_si", "parametric_si",
                              "uncertain_si", "si_from_data", "si_from_sample"
                            ),
                            config) {

  #########################################################
  # Calculates the cumulative incidence over time steps   #
  #########################################################

  calc_incidence_per_time_step <- function(incid, t_start, t_end) {
    nb_time_periods <- length(t_start)
    incidence_per_time_step <- vnapply(seq_len(nb_time_periods), function(i) {
      sum(incid[seq(t_start[i], t_end[i]), c("local", "imported")])
    })
    return(incidence_per_time_step)
  }

  #########################################################
  # Calculates the parameters of the Gamma posterior      #
  # distribution from the discrete SI distribution        #
  #########################################################

  posterior_from_si_distr <- function(incid, si_distr, a_prior, b_prior,
                                      t_start, t_end) {
    nb_time_periods <- length(t_start)
    lambda <- overall_infectivity(incid, si_distr)
    final_mean_si <- sum(si_distr * (seq(0, length(si_distr) -
                                           1)))
    a_posterior <- vector()
    b_posterior <- vector()
    a_posterior <- vnapply(seq_len(nb_time_periods), function(t) if (t_end[t] >
                                                                     final_mean_si) {
      a_prior + sum(incid[seq(t_start[t], t_end[t]), "local"])
      ## only counting local cases on the "numerator"
    } else {
      NA
    })
    b_posterior <- vnapply(seq_len(nb_time_periods), function(t) if (t_end[t] >
                                                                     final_mean_si) {
      1 / (1 / b_prior + sum(lambda[seq(t_start[t], t_end[t])]))
    } else {
      NA
    })
    return(list(a_posterior, b_posterior))
  }

  #########################################################
  # Samples from the Gamma posterior distribution for a   #
  # given mean SI and std SI                              #
  #########################################################

  sample_from_posterior <- function(sample_size, incid, mean_si, std_si,
                                    si_distr = NULL,
                                    a_prior, b_prior, t_start, t_end) {
    nb_time_periods <- length(t_start)

    if (is.null(si_distr)) {
      si_distr <- discr_si(seq(0, T - 1), mean_si, std_si)
    }

    final_mean_si <- sum(si_distr * (seq(0, length(si_distr) -
                                           1)))
    lambda <- overall_infectivity(incid, si_distr)
    a_posterior <- vector()
    b_posterior <- vector()
    a_posterior <- vnapply(seq_len(nb_time_periods), function(t) if (t_end[t] >
                                                                     final_mean_si) {
      a_prior + sum(incid[seq(t_start[t], t_end[t]), "local"])
      ## only counting local cases on the "numerator"
    } else {
      NA
    })
    b_posterior <- vnapply(seq_len(nb_time_periods), function(t) if (t_end[t] >
                                                                     final_mean_si) {
      1 / (1 / b_prior + sum(lambda[seq(t_start[t], t_end[t])], na.rm = TRUE))
    } else {
      NA
    })
    sample_r_posterior <- vapply(seq_len(nb_time_periods), function(t)
      if (!is.na(a_posterior[t])) {
        rgamma(sample_size,
               shape = unlist(a_posterior[t]),
               scale = unlist(b_posterior[t])
        )
      } else {
        rep(NA, sample_size)
      }, numeric(sample_size))
    if (sample_size == 1L) {
      sample_r_posterior <- matrix(sample_r_posterior, nrow = 1)
    }
    return(list(sample_r_posterior, si_distr))
  }

  method <- match.arg(method)
  incid <- process_I(incid)
  idx_raw_incid <- as.integer(rownames(incid)) > 0
  T <- sum(idx_raw_incid)
  T_imputed <- nrow(incid) - T
  
  check_times(config$t_start, config$t_end, T)
  nb_time_periods <- length(config$t_start)
  t_start_imputed <- config$t_start + T_imputed
  t_end_imputed <- config$t_end + T_imputed

  a_prior <- (config$mean_prior / config$std_prior)^2
  b_prior <- config$std_prior^2 / config$mean_prior

  if (method == "si_from_sample") {
    if (is.null(config$n2)) {
      stop("method si_from_sample requires to specify the config$n2 argument.")
    }
    si_sample <- process_si_sample(si_sample)
  }

  min_nb_cases_per_time_period <- ceiling(1 / config$cv_posterior^2 - a_prior)
  incidence_per_time_step <- calc_incidence_per_time_step(
    incid,
    t_start_imputed, 
    t_end_imputed
  )
  if (incidence_per_time_step[1] < min_nb_cases_per_time_period) {
    warning("You're estimating R too early in the epidemic to get the desired
            posterior CV.")
  }

  if (method == "non_parametric_si") {
    si_uncertainty <- "N"
    parametric_si <- "N"
  }
  if (method == "parametric_si") {
    si_uncertainty <- "N"
    parametric_si <- "Y"
  }
  if (method == "uncertain_si") {
    si_uncertainty <- "Y"
    parametric_si <- "Y"
  }
  if (method %in% c("si_from_data", "si_from_sample")) {
    si_uncertainty <- "Y"
    parametric_si <- "N"
  }
  if (si_uncertainty == "Y") {
    if (parametric_si == "Y") {
      mean_si_sample <- rep(-1, config$n1)
      std_si_sample <- rep(-1, config$n1)
      for (k in seq_len(config$n1)) {
        while (mean_si_sample[k] < config$min_mean_si || mean_si_sample[k] >
               config$max_mean_si) {
          mean_si_sample[k] <- rnorm(1,
                                     mean = config$mean_si,
                                     sd = config$std_mean_si
          )
        }
        while (std_si_sample[k] < config$min_std_si || std_si_sample[k] >
               config$max_std_si) {
          std_si_sample[k] <- rnorm(1, mean = config$std_si,
                                    sd = config$std_std_si)
        }
      }
      temp <- lapply(seq_len(config$n1), function(k) { sample_from_posterior(config$n2,
                                                                           incid, mean_si_sample[k], std_si_sample[k],
                                                                           si_distr = NULL, a_prior,
                                                                           b_prior, t_start_imputed, t_end_imputed
      )})
      config$si_distr <- cbind(
        t(vapply(seq_len(config$n1), function(k) (temp[[k]])[[2]], numeric(T))),
        rep(0, config$n1)
      )
      r_sample <- matrix(NA, config$n2 * config$n1, nb_time_periods)
      for (k in seq_len(config$n1)) {
        r_sample[seq((k - 1) * config$n2 + 1, k * config$n2), which(T_imputed + config$t_end >
                                                                      mean_si_sample[k])] <- (temp[[k]])[[1]][, which(T_imputed + config$t_end >
                                                                                                                        mean_si_sample[k])]
      }
      mean_posterior <- colMeans(r_sample, na.rm = TRUE)
      std_posterior <- apply(r_sample, 2, sd, na.rm = TRUE)
      quantile_0.025_posterior <- apply(r_sample, 2, quantile,
                                        0.025,
                                        na.rm = TRUE
      )
      quantile_0.05_posterior <- apply(r_sample, 2, quantile,
                                       0.05,
                                       na.rm = TRUE
      )
      quantile_0.25_posterior <- apply(r_sample, 2, quantile,
                                       0.25,
                                       na.rm = TRUE
      )
      median_posterior <- apply(r_sample, 2, median, na.rm = TRUE)
      quantile_0.75_posterior <- apply(r_sample, 2, quantile,
                                       0.75,
                                       na.rm = TRUE
      )
      quantile_0.95_posterior <- apply(r_sample, 2, quantile,
                                       0.95,
                                       na.rm = TRUE
      )
      quantile_0.975_posterior <- apply(r_sample, 2, quantile,
                                        0.975,
                                        na.rm = TRUE
      )
    } else {
      config$n1 <- dim(si_sample)[2]
      mean_si_sample <- rep(-1, config$n1)
      std_si_sample <- rep(-1, config$n1)
      for (k in seq_len(config$n1)) {
        mean_si_sample[k] <- sum((seq_len(dim(si_sample)[1]) - 1) *
                                   si_sample[, k])
        std_si_sample[k] <- sqrt(sum(si_sample[, k] *
                                       ((seq_len(dim(si_sample)[1]) - 1) -
                                          mean_si_sample[k])^2))
      }
      temp <- lapply(seq_len(config$n1), function(k) sample_from_posterior(config$n2,
                                                                           incid,
                                                                           mean_si = NULL, std_si = NULL, si_sample[, k], a_prior,
                                                                           b_prior, t_start_imputed, t_end_imputed
      ))
      config$si_distr <- cbind(
        t(vapply(seq_len(config$n1), function(k) (temp[[k]])[[2]],
                 numeric(nrow(si_sample)))),
        rep(0, config$n1)
      )
      r_sample <- matrix(NA, config$n2 * config$n1, nb_time_periods)
      for (k in seq_len(config$n1)) {
        r_sample[seq((k - 1) * config$n2 + 1,k * config$n2), which(T_imputed + config$t_end >
                                                                     mean_si_sample[k])] <- (temp[[k]])[[1]][, which(T_imputed + config$t_end >
                                                                                                                       mean_si_sample[k])]
      }
      mean_posterior <- colMeans(r_sample, na.rm = TRUE)
      std_posterior <- apply(r_sample, 2, sd, na.rm = TRUE)
      quantile_0.025_posterior <- apply(r_sample, 2, quantile,
                                        0.025,
                                        na.rm = TRUE
      )
      quantile_0.05_posterior <- apply(r_sample, 2, quantile,
                                       0.05,
                                       na.rm = TRUE
      )
      quantile_0.25_posterior <- apply(r_sample, 2, quantile,
                                       0.25,
                                       na.rm = TRUE
      )
      median_posterior <- apply(r_sample, 2, median, na.rm = TRUE)
      quantile_0.75_posterior <- apply(r_sample, 2, quantile,
                                       0.75,
                                       na.rm = TRUE
      )
      quantile_0.95_posterior <- apply(r_sample, 2, quantile,
                                       0.95,
                                       na.rm = TRUE
      )
      quantile_0.975_posterior <- apply(r_sample, 2, quantile,
                                        0.975,
                                        na.rm = TRUE
      )
    }
  } else {
    # CertainSI
    if (parametric_si == "Y") {
      config$si_distr <- discr_si(seq(0,T - 1), config$mean_si, config$std_si)
    }
    if (length(config$si_distr) < T + 1) {
      config$si_distr[seq(length(config$si_distr) + 1,T + 1)] <- 0
    }
    final_mean_si <- sum(config$si_distr * (seq(0,length(config$si_distr) -
                                                  1)))
    Finalstd_si <- sqrt(sum(config$si_distr * (seq(0,length(config$si_distr) -
                                                     1))^2) - final_mean_si^2)
    post <- posterior_from_si_distr(
      incid, config$si_distr, a_prior, b_prior,
      t_start_imputed, t_end_imputed
    )
    a_posterior <- unlist(post[[1]])
    b_posterior <- unlist(post[[2]])
    mean_posterior <- a_posterior * b_posterior
    std_posterior <- sqrt(a_posterior) * b_posterior
    quantile_0.025_posterior <- qgamma(0.025,
                                       shape = a_posterior,
                                       scale = b_posterior, lower.tail = TRUE, log.p = FALSE
    )
    quantile_0.05_posterior <- qgamma(0.05,
                                      shape = a_posterior,
                                      scale = b_posterior, lower.tail = TRUE, log.p = FALSE
    )
    quantile_0.25_posterior <- qgamma(0.25,
                                      shape = a_posterior,
                                      scale = b_posterior, lower.tail = TRUE, log.p = FALSE
    )
    median_posterior <- qgamma(0.5,
                               shape = a_posterior,
                               scale = b_posterior, lower.tail = TRUE, log.p = FALSE
    )
    quantile_0.75_posterior <- qgamma(0.75,
                                      shape = a_posterior,
                                      scale = b_posterior, lower.tail = TRUE, log.p = FALSE
    )
    quantile_0.95_posterior <- qgamma(0.95,
                                      shape = a_posterior,
                                      scale = b_posterior, lower.tail = TRUE, log.p = FALSE
    )
    quantile_0.975_posterior <- qgamma(0.975,
                                       shape = a_posterior,
                                       scale = b_posterior, lower.tail = TRUE, log.p = FALSE
    )
  }

  results <- list(R = as.data.frame(cbind(
    config$t_start, config$t_end, mean_posterior,
    std_posterior, quantile_0.025_posterior, quantile_0.05_posterior,
    quantile_0.25_posterior, median_posterior, quantile_0.75_posterior,
    quantile_0.95_posterior, quantile_0.975_posterior
  )))

  non_na_rows <- !is.na(results$R$mean_posterior)
  results$R <- results$R[non_na_rows, ]

  names(results$R) <- c(
    "t_start", "t_end", "Mean(R)", "Std(R)",
    "Quantile.0.025(R)", "Quantile.0.05(R)", "Quantile.0.25(R)",
    "Median(R)", "Quantile.0.75(R)", "Quantile.0.95(R)",
    "Quantile.0.975(R)"
  )
  results$method <- method
  results$si_distr <- config$si_distr
  if (is.matrix(results$si_distr)) {
    colnames(results$si_distr) <- paste0("t", seq(0,ncol(results$si_distr) - 1))
  } else {
    names(results$si_distr) <- paste0("t", seq(0,length(results$si_distr) - 1))
  }
  if (si_uncertainty == "Y") {
    results$SI.Moments <- as.data.frame(cbind(
      mean_si_sample,
      std_si_sample
    ))
  } else {
    results$SI.Moments <- as.data.frame(cbind(
      final_mean_si,
      Finalstd_si
    ))
  }
  names(results$SI.Moments) <- c("Mean", "Std")

  if (!is.null(incid$dates)) {
    results$dates <- check_dates(incid)
  } else {
    results$dates <- seq_len(T)
  }

  tot_incid <- rowSums(incid[, c("local", "imported")])
  results$I <- tot_incid[idx_raw_incid]
  results$I_local <- incid$local[idx_raw_incid]
  results$I_imported <- incid$imported[idx_raw_incid]

  if (T_imputed > 0) {
    results$I_imputed <- tot_incid[!idx_raw_incid]
  }

  class(results) <- "estimate_R"
  return(results)
}
