#' Description for the `incid` argument of [estimate_R()]
#' @noRd
param_incid_doc <- function(prefix = NULL, suffix = NULL) {

  incid_desc <- "One of the following
  
   - A vector (or a dataframe with a single column) of non-negative integers
     containing the incidence time series; these can be aggregated at any time
     unit as specified by argument `dt`
  
   - A dataframe of non-negative integers with either i) `incid$I`
     containing the total incidence, or ii) two columns, so that
     `incid$local` contains the incidence of cases due to local transmission
     and `incid$imported` contains the incidence of imported cases (with
     `incid$local + incid$imported` the total incidence). If the dataframe
     contains a column `incid$dates`, this is used for plotting.
     `incid$dates` must contains only dates in a row.
  
   - An object of class [incidence::incidence()]
  
   Note that the cases from the first time step are always all assumed to be
   imported cases."
  
  if (!is.null(prefix)) incid_desc <- paste(prefix, incid_desc)
  if (!is.null(suffix)) incid_desc <- paste(incid_desc, suffix)

  incid_desc
}
