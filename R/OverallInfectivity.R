#########################################################
# Calculates Lambda_t = Sum_1^t I_{t-s} * w_s           #
# with I incidence and w discrete SI distribution       #
#########################################################

#' Overall Infectivity Due To Previously Infected Individuals
#' 
#' \code{OverallInfectivity} computes the overall infectivity due to previously infected individuals. 
#' 
#' @param I One of the following
#' \itemize{
#' \item{A vector of non-negative integers containing an incidence time series}
#' \item{A dataframe of non-negative integers with two columns, so that \code{I$local} contains the incidence of cases due to local transmission and \code{I$imported} contains the incidence of imported cases (with \code{I$local + I$imported} the total incidence).}
#' } 
#' Note that the cases from the first time step are always all assumed to be imported cases. 
#' @param SI.Distr Vector of probabilities giving the discrete distribution of the serial interval.
#' @return A vector which contains the overall infectivity \eqn{\lambda_t} at each time step
#' @details{
#' The overall infectivity \eqn{\lambda_t} at time step \eqn{t} is equal to the sum of the previously infected individuals 
#' (given by the incidence vector \eqn{I}, with \code{I=I$local + I$imported} if \eqn{I} is a matrix), 
#' weigthed by their infectivity at time \eqn{t} (given by the discrete serial interval distribution \eqn{w_k}). 
#' In mathematical terms:   
#' \cr
#' \eqn{\lambda_t = \sum_{k=1}^{t-1}I_{t-k}w_k}
#' \cr
#' }
#' @seealso \code{\link{DiscrSI}}, \code{\link{EstimateR}}
#' @author Anne Cori \email{a.cori@@imperial.ac.uk} 
#' @references Cori, A. et al. A new framework and software to estimate time-varying reproduction numbers during epidemics (AJE 2013).
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
  if(is.vector(I))
  {
    I_tmp <- I
    I <- data.frame(local=I_tmp, imported=rep(0, length(I_tmp)))
    I_init <- sum(I[1,])
    I[1,] <- c(0, I_init)
  }else
  {
    if(!is.data.frame(I) | !all(c("local","imported") %in% names(I)) ) 
    {
      stop("I must be a vector or a dataframe with 2 columns called 'local' and 'imported'.")
    }
    if(I$local[1]>0)
    {
      warning("I$local[1] is >0 but must be 0, as all cases on the first time step are assumed imported. This is corrected automatically by cases being transferred to I$imported.")
      I_init <- sum(I[1,])
      I[1,] <- c(0, I_init)
    }
  }
  T<-nrow(I)
  if(any(I<0))
  {
    stop("I must contain only non negative integer values.")
  }
  if(is.vector(SI.Distr)==FALSE)
  {
    stop("SI.Distr must be a vector.")
  }
  if(SI.Distr[1]!=0)
  {
    stop("SI.Distr[1] needs to be 0.")
  }
  if(any(SI.Distr<0))
  {
    stop("SI.Distr must be a positive vector.")
  }
  if(abs(sum(SI.Distr)-1)>0.01)
  {
    stop("SI.Distr must sum to 1.")
  }
  lambda <- vector()
  lambda[1]<-NA
  for (t in 2:T)
  {
    lambda[t] <- sum(SI.Distr[1:t]*rowSums(I[t:1, c("local","imported")]),na.rm=TRUE)
  }
  return(lambda)
}
