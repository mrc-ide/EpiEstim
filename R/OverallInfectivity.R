#########################################################
# Calculates Lambda_t = Sum_1^t I_{t-s} * w_s           #
# with I incidence and w discrete SI distribution       #
#########################################################

#' Overall Infectivity Due To Previously Infected Individuals
#' 
#' \code{OverallInfectivity} computes the overall infectivity due to previously infected individuals. 
#' 
#' @param I Vector of non-negative integers containing an incidence time series.
#' @param SI.Distr Vector of probabilities giving the discrete distribution of the serial interval.
#' @return A vector which contains the overall infectivity \eqn{\lambda_t} at each time step
#' @details{
#' The overall infectivity \eqn{\lambda_t} at time step \eqn{t} is equal to the sum of the previously infected individuals 
#' (given by the incidence vector \eqn{I}), 
#' weigthed by their infectivity at time \eqn{t} (given by the discrete serial interval distribution \eqn{w_k}). 
#' In mathematical terms:   
#' \cr
#' \eqn{\lambda_t = \sum_{k=1}^{t-1}I_{t-k}w_k}
#' \cr
#' }
#' @seealso \code{\link{DiscrSI}}, \code{\link{EstimateR}
#' @author Anne Cori \email{a.cori@imperial.ac.uk} 
#' @references Cori, A. et al. A new framework and software to estimate time-varying reproduction numbers during epidemics (AJE 2013).
#' @import stats
#' @export
#' @examples
#' ## load data on pandemic flu in a school in 2009
#' data("Flu2009")
#' 
#' ## compute overall infectivity
#' lambda <- OverallInfectivity(Flu2009$Incidence, Flu2009$SI.Distr)
#' par(mfrow=c(2,1))
#' plot(Flu2009$Incidence, type="s", xlab="time (days)", ylab="Incidence")
#' title(main="Epidemic curve")
#' plot(lambda, type="s", xlab="time (days)", ylab="Infectivity")
#' title(main="Overall infectivity")
OverallInfectivity <-function (I,SI.Distr)
{
  if(is.vector(I)==FALSE)
  {
    stop("Incidence must be a vector.")
  }
  T<-length(I)
  for(i in 1:T)
  {
    if(I[i]<0)
    {
      stop("Incidence must be a positive vector.")
    }
  }
  if(is.vector(SI.Distr)==FALSE)
  {
    stop("SI.Distr must be a vector.")
  }
  if(SI.Distr[1]!=0)
  {
    stop("SI.Distr[1] needs to be 0.")
  }
  if(length(SI.Distr)>1)
  {
    for(i in 2:length(SI.Distr))
    {
      if(SI.Distr[i]<0)
      {
        stop("SI.Distr must be a positive vector.")
      }
    }
  }
  if(abs(sum(SI.Distr)-1)>0.01)
  {
    stop("SI.Distr must sum to 1.")
  }
  lambda <- vector()
  lambda[1]<-NA
  for (t in 2:length(I))
  {
    lambda[t] <- sum(SI.Distr[1:t]*I[t:1],na.rm=TRUE)
  }
  return(lambda)
}
