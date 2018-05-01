#######################################################################################################################
# init_MCMC_params finds clever starting points for the MCMC to be used to estimate the serial interval, when using option si_from_data in estimate_r #
#######################################################################################################################

#' init_MCMC_params TITLE
#' 
#' \code{init_MCMC_params} DESCRIPTION TO COME. 
#' 
#' @param si_data XXXXXXX.
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
init_MCMC_params <- function(si_data, dist = c("G", "W", "L", "off1G", "off1W", "off1L"))
{
  dist <- match.arg(dist)
  naive_SI_obs <- (si_data$SR + si_data$SL ) / 2 - (si_data$ER + si_data$EL ) / 2
  mu <- mean(naive_SI_obs)
  sigma <- sd(naive_SI_obs)
  if (dist == "G"){
    shape <- (mu/sigma)^2
    scale <- sigma^2/mu
    # check this is what we want
    # tmp <- rgamma(10000, shape=shape, scale = scale)
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
    param <- c(meanlog, sdlog)
  } else if (dist == "off1G"){
    shape <- ((mu-1)/sigma)^2
    if(shape<=0) shape <- 0.001 # this is to avoid issues when the mean SI is <1
    scale <- sigma^2/(mu-1)
    # check this is what we want
    # tmp <- 1+rgamma(10000, shape=shape, scale = scale)
    # mean(tmp)
    # sd(tmp)
    param <- c(shape, scale)
  } else if (dist == "off1W"){
    fit.w <- fitdist(naive_SI_obs-1+0.1, "weibull") # using +0.1 to avoid issues with zero
    shape <- fit.w$estimate["shape"]
    scale <- fit.w$estimate["scale"]
    # check this is what we want
    # tmp <- 1+rweibull(10000, shape=shape, scale = scale)
    # mean(tmp)
    # sd(tmp)
    param <- c(shape, scale)
  } else if (dist == "off1L"){
    sdlog <- sqrt(log(sigma^2/((mu-1)^2)+1))
    meanlog <- log(mu-1) - sdlog^2/2
    # check this is what we want
    # tmp <- 1+rlnorm(10000, meanlog=meanlog, sdlog = sdlog)
    # mean(tmp)
    # sd(tmp)
    param <- c(meanlog, sdlog)
  } else {
    stop(sprintf("Distribtion (%s) not supported",dist))
  }
  if(any(is.na(param)))
  {
    stop("NA result. Check that si_data is in the right format. ")
  }  
  return(param)
}