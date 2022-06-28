#' Aggregating daily incidence to longer time windows
#'
#' @param incid a vector of daily incidence values
#' @param dt a positive integer indicating the length of the desired aggregation window
#'
#' @return a vector of incidence values, aggregated to dt
#' @export
#'
#' @examples
#' data("SARS2003")
#' incid <- SARS2003$incidence
#' dt <- 7
#' aggregate_inc(incid, dt)
aggregate_inc <- function(incid, dt = 7L)
{
  if(dt < 2) {stop("dt should be an integer >=2")}
  if(!is.integer(dt)) {stop("dt should be an integer >=2")}
  if(!is.vector(incid)) {stop("incid should be a vector of integer values")}
  
  ndays <- length(incid)
  start <- seq(from = 1, to = ndays - (dt - 1), by=dt)
  end <- start + dt - 1
  
  weekly_inc <- numeric(length = length(start))
  
  for (i in 1:length(start)){
    weekly_inc[i] <- sum(incid[start[i]:end[i]])
  }
  
  weekly_inc <- as.integer(weekly_inc)
  weekly_inc
}