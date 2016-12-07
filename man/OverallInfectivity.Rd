\encoding{UTF-8}
\name{OverallInfectivity}
\alias{OverallInfectivity}
\title{Overall Infectivity Due To Previously Infected Individuals}
\description{
  \code{OverallInfectivity} computes the overall infectivity due to previously infected individuals. 
}
\usage{
OverallInfectivity(I, SI.Distr)
}
\arguments{
  \item{I}{vector of non-negative integers containing an incidence time series.}
  \item{SI.Distr}{vector of probabilities giving the discrete distribution of the serial interval.}
}
\value{
  \code{OverallInfectivity(I, SI.Distr)} returns a vector which contains the overall infectivity \eqn{\lambda_t} at each time step \eqn{t}.
}
\seealso{\code{\link{DiscrSI}}, \code{\link{EstimateR}}
}
\references{
  Cori, A. et al. A new framework and software to estimate time-varying reproduction numbers during epidemics (AJE 2013).
}
\details{
  The overall infectivity \eqn{\lambda_t} at time step \eqn{t} is equal to the sum of the previously infected individuals (given by the incidence vector \eqn{I}), 
  weigthed by their infectivity at time \eqn{t} (given by the discrete serial interval distribution \eqn{w_k}). 
  
  In mathematical terms:   
  \cr
  \eqn{\lambda_t = \sum_{k=1}^{t-1}I_{t-k}w_k}
 \cr
}
\author{ Anne Cori \email{a.cori@imperial.ac.uk} }
\examples{
## load data on pandemic flu in a school in 2009
data("Flu2009")

## compute overall infectivity
lambda <- OverallInfectivity(Flu2009$Incidence, Flu2009$SI.Distr)
par(mfrow=c(2,1))
plot(Flu2009$Incidence, type="s", xlab="time (days)", ylab="Incidence")
title(main="Epidemic curve")
plot(lambda, type="s", xlab="time (days)", ylab="Infectivity")
title(main="Overall infectivity")
}