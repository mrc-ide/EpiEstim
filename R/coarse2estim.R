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
    max_interval <- c(10^(-4), 1:ceiling(qgamma(0.999, shape=object@ests[1,1], scale=object@ests[2,1])))
    prob_matrix <- apply(samples, 1, function(x) dgamma(max_interval, shape=x[1], scale=x[2]))
    
  } else if (dist == "off1G"){
    # offset gamma distribution with shifted min and max value of max serial interval
    max_interval <- c(10^(-4), 1:ceiling(qgamma(0.999, shape=object@ests[1,1], scale=object@ests[2,1])))
    prob_matrix <- apply(samples, 1, function(x) dgamma(max_interval, shape=x[1], scale=x[2]))
  }
  else if (dist == "W"){
    max_interval <- c(10^(-4), 1:ceiling(qgamma(0.999, shape=object@ests[1,1], scale=object@ests[2,1])))
    prob_matrix <- apply(samples, 1, function(x) dweibull(max_interval, shape=x[1], scale=x[2]))
    
  } else if (dist == "L"){
    max_interval <- c(10^(-4), 1:ceiling(qgamma(0.999, shape=object@ests[1,1], scale=object@ests[2,1])))
    prob_matrix <- apply(samples, 1, function(x) dlnorm(max_interval, meanlog=x[1], sdlog=x[2]))
  } else {
    stop(sprintf("Distribtion (%s) not supported",dist))
  }
  prob_matrix <- apply(prob_matrix, 2, function(x) x/sum(x))
  
  prob_matrix <- rbind(rep(0, ncol(prob_matrix)), prob_matrix)
  out <- list(prob_matrix = prob_matrix, dist = dist)
  
  return(out)
}