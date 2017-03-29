#######################################################################################################################
# init_MCMC_params finds clever starting points for the MCMC to be used to estimate the serial interval, when using option SIFromData in EstimateR #
#######################################################################################################################

#' init_MCMC_params TITLE
#' 
#' \code{init_MCMC_params} DESCRIPTION TO COME. 
#' 
#' @param SI.Data XXXXXXX.
#' @param dist XXXXXXX.
#' @return XXXXXXX.
#' @details{
#' XXXXXXX. 
#' }
#' @seealso XXXXXXX.
#' @author XXXXXXX.
#' @references XXXXXXX.
#' @importFrom fitdistrplus fitdist
#' @export
#' @examples
#' ## XXXXXXX.
init_MCMC_params <- function(SI.Data, dist)
{
  naive_SI_obs <- (SI.Data$SR + SI.Data$SL ) / 2 - (SI.Data$ER + SI.Data$EL ) / 2
  mu <- mean(naive_SI_obs)
  sigma <- sd(naive_SI_obs)
  if (dist == "G" | dist == "E"){
    shape <- (mu/sigma)^2
    scale <- sigma^2/mu
    # check this is what we want
    # tmp <- rgamma(10000, shape=shape, scale = scale)
    # mean(tmp)
    # sd(tmp)
    param <- c(shape, scale)
  } else if (dist == "off1G"){
    shape <- ((mu-1)/sigma)^2
    if(shape<=0) shape <- 0.001 # this is to avoid issues when the mean SI is <1
    scale <- sigma^2/(mu-1)
    # check this is what we want
    # tmp <- 1+rgamma(10000, shape=shape, scale = scale)
    # mean(tmp)
    # sd(tmp)
    param <- c(shape, scale)
  } else if (dist == "W"){
    fit.w <- fitdist(naive_SI_obs+0.1, "weibull") # using +0.1 to avoid issues with zero
    shape <- fit.w$estimate["shape"]
    scale <- fit.w$estimate["scale"]
    # check this is what we want
    # tmp <- rweibull(10000, shape=shape, scale = scale)
    # mean(tmp)
    # sd(tmp)
    param <- c(shape, scale)
  } else if (dist == "L"){
    sdlog <- sqrt(log(sigma^2/(mu^2)+1))
    meanlog <- log(mu) - sdlog^2/2
    # check this is what we want
    # tmp <- rlnorm(10000, meanlog=meanlog, sdlog = sdlog)
    # mean(tmp)
    # sd(tmp)
    param <- c(meanlog, param)
  } else {
    stop(sprintf("Distribtion (%s) not supported",dist))
  }
  if(any(is.na(param)))
  {
    error("NA result. Check that SI.Data is in the right format. ")
  }  
  return(param)
}