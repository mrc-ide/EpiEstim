## Function integrating CoarseDataTools with EpiEstim using the amended version of EstimateR called EsimateRAmended2

# n_samples is the number of samples that are drawn from the MCMC iterations from dic.fit.mcmc

# coarse2estim is a function reformating an object from a dic.fit.mcmc function in order for it to fit the input in the EpiEstim function
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


EpiCDT <- function(..., CDT = NULL) {
  if (!is.null(CDT)) {
    c2e <- coarse2estim(CDT)
    EstimateRAmended2(..., method = c("NonParametricUncertainSI"),
                      SI.Dist.Matrix = c2e$prob_matrix)
  } else {
    EstimateRAmended2(...)
  }
}







