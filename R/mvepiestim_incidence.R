#' Process incidence input for multivariant analyses
#' 
#' Process incidence input for multivariant analyses with [estimate_advantage()]
#'
#' @param incid a multidimensional array containing values of the incidence
#'   for each time step (1st dimension), location (2nd dimension) and
#'   pathogen/strain/variant (3rd dimension)
#'
#' @param incid_imported an optional multidimensional array containing values
#'   of the incidence of imported cases
#'   for each time step (1st dimension), location (2nd dimension) and
#'   pathogen/strain/variant (3rd dimension). `incid - incid_imported` is
#'   therefore the incidence of locally infected cases. If `incid_imported` is
#'   `NULL` this means there are no
#'   known imported cases and all cases other than on those from the first
#'   time step will be considered locally infected.
#'
#' @return a list with two multidimensional elements each with three dimensions:
#'  timestep, location and pathogen/strain/variant:
#' - `local`: an array of the incidence of locally infected cases
#' - `imported`: an array of the incidence of imported cases
#'
#' @export
#'
#' @examples
#' n_v <- 3 # 3 variants
#' n_loc <- 1 # 1 location
#' T <- 100 # 100 time steps
#' # constant incidence 10 per day everywhere
#' incid <- array(10, dim = c(T, n_loc, n_v))
#' process_I_multivariant(incid)

process_I_multivariant <- function(incid, incid_imported = NULL) {
  if (is.null(incid_imported)) {
    incid_imported <- incid
    ## only cases at first time step are imported
    incid_imported[-1, , ] <- 0
  }
  dim1 <- dim(incid)
  dim2 <- dim(incid_imported)
  if (length(dim1) != length(dim2) || !all(dim1 == dim2)) {
    stop("'incid' and 'incid_imported' have incompatible dimensions")
  }
  incid_local <- incid - incid_imported
  res <- list(local = incid_local, imported = incid_imported)
  class(res) <- "incid_multivariant"
  res
}
