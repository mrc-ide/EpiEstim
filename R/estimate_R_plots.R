#' Wrapper for plot.estimate_R
#'
#' This wrapper has been created so that several \code{estimate_R} objects can 
#' be plotted at the same time. 
#'
#' @param ... Arguments of 
#'   \code{\link{plot.estimate_R}}, but in addition,
#'   parameter \code{x} can be a objects of class \code{estimate_R} (obtained as 
#'   outputs of functions \code{\link{estimate_R}} or 
#'   \code{\link{wallinga_teunis}}.  
#'   If \code{x} is a list, and \code{what='R'} or \code{what='all'}, 
#'   all estimates of R are plotted on a
#'   single graph. This will only work if all the \code{estimate_R} objects in 
#'   the list were computed using the same \code{config$t_start} and 
#'   \code{config$t_end}
#'
#' @param legend A boolean (TRUE by default) governing the presence / absence of
#'   legends on the plots
#'   
#' @return a plot (if \code{what = "incid"}, \code{"R"}, or \code{"SI"}) or a
#'   \code{\link[grid]{grob}} object (if \code{what = "all"}).
#'
#' @seealso \code{\link{plot.estimate_R}}
#'
#' @author Anne Cori, Zhian Kamvar
#'
#' @export
#'
#' @examples
#' ## load data on pandemic flu in a school in 2009
#' data("Flu2009")
#'
#' #### COMPARE THE INSTANTANEOUS AND CASE REPRODUCTION NUMBERS ####
#'
#' ## estimate the instantaneous reproduction number
#' ## (method "non_parametric_si")
#' R_instantaneous <- estimate_R(Flu2009$incidence,
#'                   method = "non_parametric_si",
#'                   config = list(t_start = seq(2, 26), 
#'                                 t_end = seq(8, 32), 
#'                                 si_distr = Flu2009$si_distr
#'                                )
#'                  )
#'
#' ## estimate the case reproduction number
#' R_case <- wallinga_teunis(Flu2009$incidence,
#'                   method = "non_parametric_si",
#'                   config = list(t_start = seq(2, 26), 
#'                                 t_end = seq(8, 32), 
#'                                 si_distr = Flu2009$si_distr
#'                   )
#'                  )
#'
#' ## visualise R estimates on the same plot
#' estimate_R_plots(list(R_instantaneous, R_case), what = "R",
#'                  options_R = list(col = c("blue", "red")), legend = TRUE)
#'                  
#' #### COMPARE THE INSTANTANEOUS R ON SLIDING WEEKLY OR BIWEEKLY WINDOWS ####
#'
#' R_weekly <- estimate_R(Flu2009$incidence,
#'                   method = "non_parametric_si",
#'                   config = list(t_start = seq(9, 26), 
#'                                 t_end = seq(15, 32), 
#'                                 si_distr = Flu2009$si_distr
#'                                )
#'                  )
#'
#' R_biweekly <- estimate_R(Flu2009$incidence,
#'                   method = "non_parametric_si",
#'                   config = list(t_start = seq(2, 19), 
#'                                 t_end = seq(15, 32),  
#'                                 si_distr = Flu2009$si_distr
#'                   )
#'                  )
#'
#' ## visualise R estimates on the same plot
#' estimate_R_plots(list(R_weekly, R_biweekly), what = "R",
#'                  options_R = list(col = c("blue", "red")), legend = TRUE)
estimate_R_plots <- function(..., legend = FALSE) {
  plot.estimate_R(..., legend = legend)
}
