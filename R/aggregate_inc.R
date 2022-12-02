#' Aggregating daily incidence to longer time windows
#'
#' @param incid a vector of daily incidence values
#' @param dt a positive integer, or vector thereof, indicating the length(s) of 
#' the desired aggregation window(s). If a vector, this will be recycled. For 
#' example, \code{dt = c(3L, 4L)} would correspond to alternating aggregation 
#' windows of 3 and 4 days
#'
#' @return a vector of incidence values, aggregated to dt
#' @export
#'
#' @examples
#' ## Constant aggregation e.g. weekly reporting
#' data("SARS2003")
#' incid <- SARS2003$incidence
#' dt <- 7L
#' aggregate_inc(incid, dt)
#' 
#' ## Non-constant aggregation e.g. reporting 3x week
#' #' data("SARS2003")
#' incid <- SARS2003$incidence
#' dt <- c(2L,2L,3L)
#' aggregate_inc(incid, dt)
#' 
aggregate_inc <- function(incid, dt = 7L)
{
  if(all(dt < 2)) {stop("at least one value of dt should be an integer >=2")}
  if(any(!is.integer(dt))) {stop("dt should be an integer or vector of integers e.g. 2L or c(2L,2L,3L)")}
  if(!is.vector(incid)) {stop("incid should be a vector of integer values")}
  
  ndays <- length(incid)
  
  if (length(dt) == 1){
    start <- seq(from = 1, to = ndays - (dt - 1), by=dt)
    end <- start + dt - 1
  } else {
    start <- cumsum(c(1, rep(dt, length.out = ndays/sum(dt)*length(dt))))
    end <- cumsum(rep(dt, length.out = length(start)))
    if(end[length(end)]>ndays){
      start <- start[-length(start)]
      end <- end[-length(end)]
    }
  }
  
  message("Incidence aggregated up to day ", end[length(end)], " of ", ndays)
  
  agg_inc <- numeric(length = length(start))
  
  for (i in 1:length(start)){
    agg_inc[i] <- sum(incid[start[i]:end[i]])
  }
  
  agg_inc <- as.integer(paste(agg_inc))
  agg_inc
}