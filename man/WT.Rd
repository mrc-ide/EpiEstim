\encoding{UTF-8}

\name{WT}
\alias{WT}
\title{Estimation of the case reproduction number using the Wallinga and Teunis method}
\description{
  \code{WT} estimates the case reproduction number of an epidemic, given the incidence time series and the serial interval distribution. 
}
\usage{
WT(I, T.Start, T.End, 
	method=c("NonParametricSI","ParametricSI"), 
    Mean.SI=NULL, Std.SI=NULL,
    SI.Distr=NULL, nSim=10, 
    plot=FALSE, leg.pos="topright")
}

\arguments{
  \item{I}{vector of non-negative integers containing the incidence time series.}
  \item{T.Start,T.End}{vectors of positive integers giving the starting end ending times of each window over which the reproduction number will be estimated. These must be in ascending order, and so that for all \code{i}, \code{T.Start[i]<=T.End[i]}. T.Start[1] should be strictly after the first day with non null incidence.}
  \item{method}{one of "NonParametricSI" or "ParametricSI" (see details).}
  \item{Mean.SI}{for method "ParametricSI" ; positive real giving the mean serial interval.}
  \item{Std.SI}{for method "ParametricSI" ; non negative real giving the stadard deviation of the serial interval.}
  \item{SI.Distr}{for method "NonParametricSI" ; vector of probabilities giving the discrete distribution of the serial interval, starting with \code{SI.Distr[1]} (probability that the serial interval is zero), which should be zero.}
  \item{nSim}{a positive integer giving the number of simulated epidemic trees used for computation of the confidence intervals of the case reproduction number (see details).}
  \item{plot}{logical. If \code{TRUE} (default is \code{FALSE}), output is plotted (see value).}
  \item{leg.pos}{one of "\code{bottomright}", "\code{bottom}", "\code{bottomleft}", "\code{left}", "\code{topleft}", "\code{top}", "\code{topright}", "\code{right}", "\code{center}" or \code{\link{xy.coords}(x, y)}, with \code{x} and \code{y} real numbers. 
  This specifies the position of the legend in the plot. Alternatively, \code{locator(1)} can be used ; the user will then need to click where the legend needs to be written.}
}
\value{
	a list with components: 
	\item{R}{a dataframe containing: 
		the times of start and end of each time window considered ; 
		the estimated mean, std, and 0.025 and 0.975 quantiles of the reproduction number for each time window.}
	\item{SIDistr}{a dataframe containing: 
		for method "NonParametricSI", the mean and standard deviation of the discrete serial interval distribution;
		for method "ParametricSI", the discrete serial interval distribution.}
}
\seealso{\code{\link{DiscrSI}}, \code{\link{EstimateR}}}
\references{
  Cori, A. et al. A new framework and software to estimate time-varying reproduction numbers during epidemics. (submitted) 
  Wallinga, J. and P. Teunis (2004). Different epidemic curves for severe acute respiratory syndrome reveal similar impacts of control measures. Am J Epidemiol 160(6): 509-516.
}
\details{
  Estimates of the case reproduction number for an epidemic over predefined time windows can be obtained, 
  for a given discrete distribution of the serial interval, as proposed by Wallinga and Teunis (AJE, 2004). 
  Confidence intervals are obtained by simulating a number (nSim) of possible transmission trees. 
  
  The methods vary in the way the serial interval distribution is specified.
  
  ----------------------- \code{method "NonParametricSI"} -----------------------
  
  The discrete distribution of the serial interval is directly specified in the argument \code{SI.Distr}.
  
  If \code{plot} is \code{TRUE}, 3 plots are produced. 
  The first one shows the epidemic curve. 
  The second one shows the posterior median and 95\% credible interval of the reproduction number. The estimate for a time window is plotted at the end of the time window. 
  The position of the legend on that graph can be monitored by the argument \code{leg.pos} (default is "\code{topright}").
  The third plot shows the discrete distribution of the serial interval. 
  
  ----------------------- \code{method "ParametricSI"} -----------------------
  
  The mean and standard deviation of the continuous distribution of the serial interval are given in the arguments \code{Mean.SI} and \code{Std.SI}.
  The discrete distribution of the serial interval is derived automatically using \code{\link{DiscrSI}}.
  
  If \code{plot} is \code{TRUE}, 3 plots are produced, which are identical to the ones for \code{method "NonParametricSI"} .
  
}
\author{ Anne Cori \email{a.cori@imperial.ac.uk} }
\examples{
## load data on pandemic flu in a school in 2009
data("Flu2009")

## estimate the case reproduction number (method "NonParametricSI")
WT(Flu2009$Incidence, T.Start=2:26, T.End=8:32, method="NonParametricSI", 
    SI.Distr=Flu2009$SI.Distr, plot=TRUE, leg.pos=xy.coords(1,1.75), nSim=100)
# the second plot produced shows, at each each day, 
# the estimate of the cqse reproduction number over the 7-day window finishing on that day.

## estimate the case reproduction number (method "ParametricSI")
WT(Flu2009$Incidence, T.Start=2:26, T.End=8:32, method="ParametricSI", 
    Mean.SI=2.6, Std.SI=1.5, plot=TRUE, nSim=100)
# the second plot produced shows, at each each day, 
# the estimate of the cqse reproduction number over the 7-day window finishing on that day.

}
\keyword{distribution}
