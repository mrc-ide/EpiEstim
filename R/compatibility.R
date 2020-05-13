#' Function to ensure compatibility with EpiEstim versions <2.0
#'
#' Please only use for compatibility;
#' Prefer the new estimate_R function instead
#'
#' @param I see \code{incid} in \code{estimate_R}
#' @param T.Start see \code{config$t_start} in \code{estimate_R}
#' @param T.End see \code{config$t_end} in \code{estimate_R}
#' @param method see method in \code{estimate_R} (but EstimateR uses CamelCase 
#' where estimate_R uses snake_case for the method names)
#' @param n1 see \code{n1} in \code{estimate_R}
#' @param n2 see \code{n2} in \code{estimate_R}
#' @param Mean.SI see \code{config$mean_si} in \code{estimate_R}
#' @param Std.SI see \code{config$std_si} in \code{estimate_R}
#' @param Std.Mean.SI see \code{config$std_mean_si} in \code{estimate_R}
#' @param Min.Mean.SI see \code{config$min_mean_si} in \code{estimate_R}
#' @param Max.Mean.SI see \code{config$max_mean_si} in \code{estimate_R}
#' @param Std.Std.SI see \code{config$std_std_si} in \code{estimate_R}
#' @param Min.Std.SI see \code{config$min_std_si} in \code{estimate_R}
#' @param Max.Std.SI see \code{config$max_std_si} in \code{estimate_R}
#' @param SI.Distr see \code{config$si_distr} in \code{estimate_R}
#' @param Mean.Prior see \code{config$mean_prior} in \code{estimate_R}
#' @param Std.Prior see \code{config$std_prior} in \code{estimate_R}
#' @param CV.Posterior see \code{config$cv_posterior} in \code{estimate_R}
#' @param plot Not used anymore, only there for compatibility
#' @param leg.pos Not used anymore, only there for compatibility
#'
#' @export
#'
EstimateR <- function(I, T.Start, T.End,
                      method = c("NonParametricSI", "ParametricSI", 
                                 "UncertainSI"),
                      n1 = NULL, n2 = NULL, Mean.SI = NULL, Std.SI = NULL,
                      Std.Mean.SI = NULL, Min.Mean.SI = NULL, 
                      Max.Mean.SI = NULL,
                      Std.Std.SI = NULL, Min.Std.SI = NULL, Max.Std.SI = NULL,
                      SI.Distr = NULL,
                      Mean.Prior = 5, Std.Prior = 5, CV.Posterior = 0.3,
                      plot = FALSE, leg.pos = "topright") {
  .Deprecated("estimate_R")

  method_tr <- c(
    "NonParametricSI" = "non_parametric_si",
    "ParametricSI" = "parametric_si",
    "UncertainSI" = "uncertain_si"
  )
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
    mean_prior = Mean.Prior,
    std_prior = Std.Prior, ### TO DO: change to sd_prior
    cv_posterior = CV.Posterior
  )

  estimate_R(
    incid = I, method = method,
    config = config
  )
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
DiscrSI <- function(k, mu, sigma) {
  .Deprecated("discr_si")
  discr_si(k, mu, sigma)
}



#' Function to ensure compatibility with EpiEstim versions <2.0
#'
#' Please only use for compatibility;
#' Prefer the new overall_infectivity function instead
#'
#' @param I see \code{incid} in \code{overall_infectivity}
#' @param SI.Distr see \code{si_distr} in \code{overall_infectivity}

#' @export
OverallInfectivity <- function(I, SI.Distr) {
  # We will eventually properly deprecate this, but it is being used in other 
  # packages at the moment, so soft deprecate for now.
  if (interactive()) {
    .Deprecated("overall_infectivity")
  }
  overall_infectivity(incid = I, si_distr = SI.Distr)
}




#' Function to ensure compatibility with EpiEstim versions <2.0
#'
#' Please only use for compatibility;
#' Prefer the new wallinga_teunis function instead
#'
#' @param I see \code{incid} in \code{wallinga_teunis}
#' @param T.Start see \code{config$t_start} in \code{wallinga_teunis}
#' @param T.End see \code{config$t_end} in \code{wallinga_teunis}
#' @param method see method in \code{wallinga_teunis} (but WT uses CamelCase 
#' where wallinga_teunis uses snake_case for the method names)
#' @param Mean.SI see \code{config$mean_si} in \code{wallinga_teunis}
#' @param Std.SI see \code{config$std_si} in \code{wallinga_teunis}
#' @param SI.Distr see \code{config$si_distr} in \code{wallinga_teunis}
#' @param nSim see \code{config$n_sim} in \code{wallinga_teunis}
#' @param plot Not used anymore, only there for compatibility
#' @param leg.pos Not used anymore, only there for compatibility
#'
#' @export
WT <- function(I, T.Start, T.End,
               method = c("NonParametricSI", "ParametricSI"),
               Mean.SI = NULL, Std.SI = NULL,
               SI.Distr = NULL, nSim = 10,
               plot = FALSE, leg.pos = "topright") {
  .Deprecated("wallinga_teunis")
  method_tr <- c(
    "NonParametricSI" = "non_parametric_si",
    "ParametricSI" = "parametric_si"
  )
  method <- method_tr[[method]]

  config <- list(
    t_start = T.Start,
    t_end = T.End,
    mean_si = Mean.SI,
    std_si = Std.SI,
    si_distr = SI.Distr,
    n_sim = nSim
  )

  wallinga_teunis(incid = I, method = method, config = config)
}
