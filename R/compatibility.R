#' Function to ensure compatibility with EpiEstim versions <2.0
#'
#' Please only use for compatibility;
#' Prefer the new estimate_r function instead
#'
#' @param I see \code{I} in \code{estimate_r}
#' @param T.Start see \code{config$t_start} in \code{estimate_r}
#' @param T.End see \code{config$t_end} in \code{estimate_r}
#' @param method see method in \code{estimate_r} (but EstimateR uses CamelCase where estimate_r uses snake_case for the method names)
#' @param n1 see \code{n1} in \code{estimate_r}
#' @param n2 see \code{n2} in \code{estimate_r}
#' @param Mean.SI see \code{config$mean_si} in \code{estimate_r}
#' @param Std.SI see \code{config$std_si} in \code{estimate_r}
#' @param Std.Mean.SI see \code{config$std_mean_si} in \code{estimate_r}
#' @param Min.Mean.SI see \code{config$min_mean_si} in \code{estimate_r}
#' @param Max.Mean.SI see \code{config$max_mean_si} in \code{estimate_r}
#' @param Std.Std.SI see \code{config$std_std_si} in \code{estimate_r}
#' @param Min.Std.SI see \code{config$min_std_si} in \code{estimate_r}
#' @param Max.Std.SI see \code{config$max_std_si} in \code{estimate_r}
#' @param SI.Distr see \code{config$si_distr} in \code{estimate_r}
#' @param SI.Data see \code{si_data} in \code{estimate_r}
#' @param SI.parametricDistr see \code{config$si_parametric_distr} in \code{estimate_r}
#' @param MCMC.control see \code{config$mcmc_control} in \code{estimate_r}
#' @param SI.Sample see \code{si_sample} in \code{estimate_r}
#' @param seed see \code{config$seed} in \code{estimate_r}
#' @param Mean.Prior see \code{config$mean_prior} in \code{estimate_r}
#' @param Std.Prior see \code{config$std_prior} in \code{estimate_r}
#' @param CV.Posterior see \code{config$cv_posterior} in \code{estimate_r}
#' @param plot see \code{config$plot} in \code{estimate_r}
#' @param legend see \code{config$legend} in \code{estimate_r}
#'
#' @export
#'
EstimateR <- function(I, T.Start, T.End, method = c("NonParametricSI", "ParametricSI",
                                                    "UncertainSI", "SIFromData", "SIFromSample"), 
                      n1 = NULL, n2 = NULL, Mean.SI = NULL, Std.SI = NULL,
                      Std.Mean.SI = NULL, Min.Mean.SI = NULL, Max.Mean.SI = NULL,
                      Std.Std.SI = NULL, Min.Std.SI = NULL, Max.Std.SI = NULL,
                      SI.Distr = NULL, 
                      SI.Data = NULL, SI.parametricDistr = c("G", "W", "L", "off1G", "off1W", "off1L"),  
                      MCMC.control = list(init_pars = NULL, burnin = 3000, thin=10, seed = as.integer(Sys.time())), 
                      SI.Sample = NULL, 
                      seed = NULL,
                      Mean.Prior = 5, Std.Prior = 5, CV.Posterior = 0.3,
                      plot = FALSE, legend = FALSE) {
  .Deprecated("estimate_r")
  
  ### TO DO change I to incidence
  
  method_tr <- c("NonParametricSI" = "non_parametric_si", 
                 "ParametricSI" = "parametric_si",
                 "UncertainSI" = "uncertain_si", 
                 "SIFromData" = "si_from_data", 
                 "SIFromSample" = "si_from_sample")
  method <- method_tr[[method]]

  config <- list(
    t_start = T.Start,
    t_end = T.End,
    n1 = n1,
    n2 = n2,
    mean_si = Mean.SI,
    std_si = Std.SI,
    std_mean_si = Std.Mean.SI,
    min_mean_si = Min.Mean.SI,
    max_mean_si = Max.Mean.SI,
    std_std_si = Std.Std.SI,
    min_std_si = Min.Std.SI,
    max_std_si = Max.Std.SI,
    si_distr = SI.Distr,
    si_parametric_distr = SI.parametricDistr,
    mcmc_control = MCMC.control,
    seed = seed,
    mean_prior = Mean.Prior,
    std_prior = Std.Prior, ### TO DO: change to sd_prior
    cv_posterior = CV.Posterior,
    plot = plot, ### TO DO: plot outside the function (alsways!)
    legend = legend
  )
    
  estimate_r(I, method, si_data = SI.Data, si_sample = SI.Sample, 
             config = config)
}




#' Function to ensure compatibility with EpiEstim versions <2.0
#'
#' Please only use for compatibility;
#' Prefer the new discr_si function instead
#'
#' @param k see \code{k} in \code{discr_si}
#' @param mu see \code{mu} in \code{discr_si}
#' @param sigma see \code{sigma} in \code{discr_si}
#'
#' @export
#'
DiscrSI <- function(k,mu,sigma)
{
  .Deprecated("discr_si")
  discr_si(k,mu,sigma)
}



#' Function to ensure compatibility with EpiEstim versions <2.0
#'
#' Please only use for compatibility;
#' Prefer the new overall_infectivity function instead
#'
#' @param I see \code{I} in \code{overall_infectivity}
#' @param SI.Distr see \code{si_distr} in \code{overall_infectivity}

#' @export
OverallInfectivity <- function(I, SI.Distr)
{
  .Deprecated("overall_infectivity")
  overall_infectivity(I = I, si_distr = SI.Distr)
}




#' Function to ensure compatibility with EpiEstim versions <2.0
#'
#' Please only use for compatibility;
#' Prefer the new wallinga_teunis function instead
#'
#' @param I see \code{I} in \code{wallinga_teunis}
#' @param T.Start see \code{config$t_start} in \code{wallinga_teunis}
#' @param T.End see \code{config$t_end} in \code{wallinga_teunis}
#' @param method see method in \code{wallinga_teunis} (but WT uses CamelCase where wallinga_teunis uses snake_case for the method names)
#' @param Mean.SI see \code{config$mean_si} in \code{wallinga_teunis}
#' @param Std.SI see \code{config$std_si} in \code{wallinga_teunis}
#' @param SI.Distr see \code{config$si_distr} in \code{wallinga_teunis}
#' @param nSim see \code{config$n_sim} in \code{wallinga_teunis}
#' @param plot see \code{config$plot} in \code{wallinga_teunis}
#'
#' @export
WT <- function(I, T.Start, T.End,
               method=c("NonParametricSI", "ParametricSI"),
               Mean.SI=NULL, Std.SI=NULL, 
               SI.Distr=NULL, nSim=10, 
               plot=FALSE)
{
  .Deprecated("wallinga_teunis")
  method_tr <- c("NonParametricSI" = "non_parametric_si", 
                 "ParametricSI" = "parametric_si")
  method <- method_tr[[method]]
  
  config = list(t_start = T.Start, 
                t_end = T.End,
                mean_si = Mean.SI, 
                std_si = Std.SI, 
                si_distr = SI.Distr, 
                n_sim = nSim, 
                plot = plot)
  
  wallinga_teunis(I, method = method, config)
}




  
