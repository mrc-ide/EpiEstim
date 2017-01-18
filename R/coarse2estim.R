#######################################################################################################################
# coarse2estim integates CoarseDataTools with EpiEstim using the amended version of EstimateR called EstimateR_func #
#######################################################################################################################

#' coarse2estim TITLE
#' 
#' \code{coarse2estim} DESCRIPTION TO COME. 
#' 
#' @param object XXXXXXX.
#' @param n_samples XXXXXXX.
#' @return XXXXXXX.
#' @details{
#' XXXXXXX. 
#' }
#' @seealso XXXXXXX.
#' @author XXXXXXX.
#' @references XXXXXXX.
# #' @import 
#' @export
#' @examples
#' ## XXXXXXX.
coarse2estim <- function(object, n_samples=1000){
  
  samples0 <- as.matrix(object@samples)
  if(n_samples<nrow(samples0))
  {
    index <- sample(1:nrow(samples0), size= n_samples)
    samples <- samples0[index, ]
  }else
  {
    samples <- samples0
  }
  dist <- object@dist
  
  ##  Probability matrix that will be used in EpiEstim based on which distribution is specified by the user
  if (dist == "G" | dist == "E"){
  	# For each input parameter set, find the 99th percentile, and take the maximum of these as the maximum
    # serial interval that we need to consider
    maxValue <- max( sapply(1:n_samples, function(i) ceiling(qgamma(0.999, shape = samples[i,1], scale = samples[i,2])) ) )
    max_interval <- 1:maxValue
    prob_matrix <- apply(samples, 1, function(x) pgamma(max_interval+0.5, shape=x[1], scale=x[2]) - pgamma(max_interval-0.5, shape=x[1], scale=x[2]))
  } else if (dist == "off1G"){
    # offset gamma distribution with shifted min and max value of max serial interval
    maxValue <- max( sapply(1:n_samples, function(i) ceiling(qgamma(0.999, shape = samples[i,1], scale = samples[i,2])) ) )
    max_interval <- 0:maxValue
    prob_matrix <- apply(samples, 1, function(x) pgamma(max_interval+0.5, shape=x[1], scale=x[2]) - pgamma(max_interval-0.5, shape=x[1], scale=x[2]))
  } else if (dist == "W"){
  	maxValue <- max( sapply(1:n_samples, function(i) ceiling(qweibull(0.999, shape = samples[i,1], scale = samples[i,2])) ) )
    max_interval <- 1:maxValue
    prob_matrix <- apply(samples, 1, function(x) pweibull(max_interval+0.5, shape=x[1], scale=x[2]) - pweibull(max_interval-0.5, shape=x[1], scale=x[2]))
  } else if (dist == "L"){
  	maxValue <- max( sapply(1:n_samples, function(i) ceiling(qlnorm(0.999, meanlog = samples[i,1], sdlog = samples[i,2])) ) )
    max_interval <- 1:maxValue
    prob_matrix <- apply(samples, 1, function(x) plnorm(max_interval+0.5, meanlog=x[1], sdlog=x[2]) - plnorm(max_interval-0.5, meanlog=x[1], sdlog=x[2]))
  } else {
    stop(sprintf("Distribtion (%s) not supported",dist))
  }
  # adding initial 0 for P(SI=0)
  prob_matrix <- rbind(rep(0, n_samples), prob_matrix)
  # renormalising
  prob_matrix <- apply(prob_matrix, 2, function(x) x/sum(x))
  
  out <- list(prob_matrix = prob_matrix, dist = dist)
  
  return(out)
}