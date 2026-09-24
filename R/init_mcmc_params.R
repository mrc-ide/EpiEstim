#' Find starting points for the estimation of the serial interval
#'
#' This function is deprecated. Use [si_start_values()] instead.
#'
#' @inheritParams si_start_values
#' @return A vector containing the starting values for the parameters of the
#'   distribution of the serial interval, as given by [si_start_values()].
#'
#' @seealso [si_start_values()]
#'
#' @author Anne Cori
#'
#' @export
init_mcmc_params <- function(si_data, dist) {
  .Deprecated("si_start_values")
  unname(unlist(si_start_values(si_data, dist)))
}
