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
#' @import stats
#' @export
#' @examples
#' ## XXXXXXX.
coarse2estim <- function(object, n_samples=1000){
  
  samples0 <- as.matrix(object@samples)
  index <- sample(1:nrow(samples0), size= n_samples)
  samples <- samples0[index, ]
  dist <- object@dist
  
  ##  Probability matrix that will be used in EpiEstim based on which distribution is specified by the user
  if (dist == "G"){
    max_interval <- 1:ceiling(qgamma(0.999, shape=object@ests[1,1], scale=object@ests[2,1]))
    prob_matrix <- apply(samples, 1, function(x) pgamma(max_interval+0.5, shape=x[1], scale=x[2]) - pgamma(max_interval-0.5, shape=x[1], scale=x[2]))
  } else if (dist == "off1G"){
    # offset gamma distribution with shifted min and max value of max serial interval
    max_interval <- 0:ceiling(qgamma(0.999, shape=object@ests[1,1], scale=object@ests[2,1]))
    prob_matrix <- apply(samples, 1, function(x) pgamma(max_interval+0.5, shape=x[1], scale=x[2]) - pgamma(max_interval-0.5, shape=x[1], scale=x[2]))
  }
  else if (dist == "W"){
    max_interval <- 1:ceiling(qweibull(0.999, shape=object@ests[1,1], scale=object@ests[2,1]))
    prob_matrix <- apply(samples, 1, function(x) pweibull(max_interval+0.5, shape=x[1], scale=x[2]) - pweibull(max_interval-0.5, shape=x[1], scale=x[2]))
  } else if (dist == "L"){
    max_interval <- 1:ceiling(qlnorm(0.999, meanlog=object@ests[1,1], sdlog=object@ests[2,1]))
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