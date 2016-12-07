#########################################################
# Discretized serial interval (assuming a shifted gamma #
# distribution (with shift 1)                           #
#########################################################

#' Discretized Generation Time Distribution Assuming A Shifted Gamma Distribution
#' 
#' code{DiscrSI} computes the discrete distribution of the serial interval, assuming that the serial interval is shifted Gamma distributed, with shift 1. 
#' 
#' @param k Positive integer for which the discrete distribution is desired.
#' @param mu A positive real giving the mean of the Gamma distribution.
#' @param sigma A non-negative real giving the standard deviation of the Gamma distribution.
#' @return Gives the discrete probability \eqn{w_k} that the serial interval is equal to \eqn{k}.
#' @details{
#' Assuming that the serial interval is shifted Gamma distributed with mean \eqn{\mu}, standard deviation \eqn{\sigma} and shift \eqn{1}, 
#' the discrete probability \eqn{w_k} that the serial interval is equal to \eqn{k} is: 
#' \cr 
#' \eqn{w_k = kF_{\{\mu-1,\sigma\}}(k)+(k-2)F_{\{\mu-1,\sigma\}}(k-2)-2(k-1)F_{\{\mu-1,\sigma\}}(k-1)\\+(\mu-1)(2F_{\{\mu-1+\frac{\sigma^2}{\mu-1},\sigma\sqrt{1+\frac{\sigma^2}{\mu-1}}\}}(k-1)-F_{\{\mu-1+\frac{\sigma^2}{\mu-1},\sigma\sqrt{1+\frac{\sigma^2}{\mu-1}}\}}(k-2)-F_{\{\mu-1+\frac{\sigma^2}{\mu-1},\sigma\sqrt{1+\frac{\sigma^2}{\mu-1}}\}}(k))}
#' \cr 
#' where \eqn{F_{\{\mu,\sigma\}}} is the cumulative density function of a Gamma distribution with mean \eqn{\mu} and standard deviation \eqn{\sigma}.
#' }
#' @seealso \code{\link{OverallInfectivity}}, \code{\link{EstimateR}
#' @author Anne Cori \email{a.cori@imperial.ac.uk} 
#' @references Cori, A. et al. A new framework and software to estimate time-varying reproduction numbers during epidemics (AJE 2013).
#' @import stats
#' @export
#' @examples
#' ## Computing the discrete serial interval of influenza
#' MeanFluSI <- 2.6
#' SdFluSI <- 1.5
#' DicreteSIDistr <- vector()
#' for(i in 0:20)
#' {
#' DicreteSIDistr[i+1] <- DiscrSI(i, MeanFluSI, SdFluSI)
#' }
#' plot(0:20, DicreteSIDistr, type="h", lwd=10, lend=1, xlab="time (days)", ylab="frequency")
#' title(main="Discrete distribution of the serial interval of influenza")
DiscrSI<-function(k,mu,sigma)
{
  if(sigma<0)
  {
    stop("sigma must be >=0.")
  }
  a=((mu-1)/sigma)^2
  b=sigma^2/(mu-1)
  CDFGamma<-function(k,a,b)
  {
    return(pgamma(k,shape=a,scale=b))
  }
  res<-k*CDFGamma(k,a,b)+(k-2)*CDFGamma(k-2,a,b)-2*(k-1)*CDFGamma(k-1,a,b)
  res<-res+a*b*(2*CDFGamma(k-1,a+1,b)-CDFGamma(k-2,a+1,b)-CDFGamma(k,a+1,b))
  res<-max(0,res)
  return(res)
}
