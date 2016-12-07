\encoding{UTF-8}

\name{EstimateR}
\alias{EstimateR}
\title{Estimated Reproduction Number}
\description{
  \code{EstimateR} estimates the reproduction number of an epidemic, given the incidence time series and the serial interval distribution. 
}
\usage{
EstimateR(I, T.Start, T.End, 
    method=c("NonParametricSI","ParametricSI","UncertainSI","NonParametricUncertainSI"), 
    n1=NULL, n2=NULL, Mean.SI=NULL, Std.SI=NULL,
    Std.Mean.SI=NULL, Min.Mean.SI=NULL, Max.Mean.SI=NULL, 
    Std.Std.SI=NULL, Min.Std.SI=NULL, Max.Std.SI=NULL,
    SI.Distr=NULL, 
    Mean.Prior=5, Std.Prior=5, CV.Posterior=0.3, 
    plot=FALSE, leg.pos="topright",CDT = NULL)
}
\arguments{
  \item{I}{vector of non-negative integers containing the incidence time series.}
  \item{T.Start,T.End}{vectors of positive integers giving the starting end ending times of each window over which the reproduction number will be estimated. These must be in ascending order, and so that for all \code{i}, \code{T.Start[i]<=T.End[i]}. T.Start[1] should be strictly after the first day with non null incidence.}
  \item{method}{one of "NonParametricSI", "ParametricSI" or "UncertainSI" (see details).}
  \item{n1}{for method "UncertainSI" ; positive integer giving the size of the sample of pairs (Mean SI (serial interval), Std SI) to be drawn (see details).}
  \item{n2}{for method "UncertainSI" ; positive integer giving the size of the sample drawn from each posterior distribution conditional to a pair (Mean SI, Std SI) (see details).}
  \item{Mean.SI}{for method "ParametricSI" and "UncertainSI" ; positive real giving the mean serial interval (method "ParametricSI") or the average mean serial interval (method "UncertainSI", see details).}
  \item{Std.SI}{for method "ParametricSI" and "UncertainSI" ; non negative real giving the stadard deviation of the serial interval (method "ParametricSI") or the average standard deviation of the serial interval (method "UncertainSI", see details).}
  \item{Std.Mean.SI}{for method "UncertainSI" ; standard deviation of the distribution from which mean serial intervals are drawn (see details).}
  \item{Min.Mean.SI}{for method "UncertainSI" ; lower bound of the distribution from which mean serial intervals are drawn (see details).}
  \item{Max.Mean.SI}{for method "UncertainSI" ; upper bound of the distribution from which mean serial intervals are drawn (see details).}
  \item{Std.Std.SI}{for method "UncertainSI" ; standard deviation of the distribution from which standard deviations of the serial interval are drawn (see details).}
  \item{Min.Std.SI}{for method "UncertainSI" ; lower bound of the distribution from which standard deviations of the serial interval are drawn (see details).}
  \item{Max.Std.SI}{for method "UncertainSI" ; upper bound of the distribution from which standard deviations of the serial interval are drawn (see details).}
  \item{SI.Distr}{for method "NonParametricSI" ; vector of probabilities giving the discrete distribution of the serial interval, starting with \code{SI.Distr[1]} (probability that the serial interval is zero), which should be zero.}
  \item{Mean.Prior}{a positive number giving the mean of the common prior distribution for all reproduction numbers (see details).}
  \item{Std.Prior}{a positive number giving the standard deviation of the common prior distribution for all reproduction numbers (see details).}
  \item{CV.Posterior}{a positive number giving the aimed posterior coefficient of variation (see details).}
  \item{plot}{logical. If \code{TRUE} (default is \code{FALSE}), output is plotted (see value).}
  \item{leg.pos}{one of "\code{bottomright}", "\code{bottom}", "\code{bottomleft}", "\code{left}", "\code{topleft}", "\code{top}", "\code{topright}", "\code{right}", "\code{center}" or \code{\link{xy.coords}(x, y)}, with \code{x} and \code{y} real numbers. 
  This specifies the position of the legend in the plot. Alternatively, \code{locator(1)} can be used ; the user will then need to click where the legend needs to be written.}
  \item{CDT}{for method "NonParametricUncertainSI" ; an object of the S4 class \code{coarseDataTools::cd.fit.mcmc} which will be used to describe the model used to estimate the SI distribution.}
}
\value{
	a list with components: 
	\item{R}{a dataframe containing: 
		the times of start and end of each time window considered ; 
		the posterior mean, std, and 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975 quantiles of the reproduction number for each time window.}
	\item{SIDistr}{a dataframe containing: 
		for method "NonParametricSI", the mean and standard deviation of the discrete serial interval distribution;
		for method "ParametricSI", the discrete serial interval distribution;
		for method "UncertainSI", the means and standard deviations of the serial interval are sampled from truncated normal distributions to account for uncertainty on the serial interval distribution (see details);
		for method "NonParametricUncertainSI", the distribution of the serial interval is sampled directly from the posterior distribution obtained by Bayesian estimation of the serial interval distribution from doubly censored data (relying on the R package coarseDataTools) to account for uncertainty on the serial interval distribution (see details).}
}
\seealso{\code{\link{OverallInfectivity}}, \code{\link{DiscrSI}}
}
\references{
  Cori, A. et al. A new framework and software to estimate time-varying reproduction numbers during epidemics (AJE 2013) .
}
\details{
  Analytical estimates of the reproduction number for an epidemic over predefined time windows can be obtained within a Bayesian framework, 
  for a given discrete distribution of the serial interval (see references). 
  
  The more incident cases are observed over a time window, the smallest the posterior coefficient of variation (CV, ratio of standard deviation over mean) of the reproduction number. 
  An aimed CV can be specified in the argument \code{CV.Posterior} (default is \code{0.3}), and a warning will be produced if the incidence within one of the time windows considered is too low to get this CV. 
  
  The methods vary in the way the serial interval distribution is specified. The plots are also different according to the method used.
  
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
  
  ----------------------- \code{method "UncertainSI"} -----------------------
  
  \code{Method "UncertainSI"} allows accounting for uncertainty on the serial interval distribution (see references). 
  We allow the mean \eqn{\mu} and standard deviation \eqn{\sigma} of the serial interval to vary according to truncated normal distributions. 
  We sample \code{n1} pairs of mean and standard deviations, \eqn{(\mu^{(1)},\sigma^{(1)}),...,(\mu^{(n_2)},\sigma^{(n_2)})}, by first sampling the mean \eqn{\mu^{(k)}} 
  from its truncated normal distribution (with mean \code{Mean.SI}, standard deviation \code{Std.Mean.SI}, minimum \code{Min.Mean.SI} and maximum \code{Max.Mean.SI}), 
  and then sampling the standard deviation \eqn{\sigma^{(k)}} from its truncated normal distribution 
  (with mean \code{Std.SI}, standard deviation \code{Std.Std.SI}, minimum \code{Min.Std.SI} and maximum \code{Max.Std.SI}), but imposing that \eqn{\sigma^{(k)}<\mu^{(k)}}. 
  This constraint ensures that the Gamma probability density function of the serial interval is null at \eqn{t=0}. 
  Warnings are produced when the truncated normal distributions are not symmetric around the mean. 
  For each pair \eqn{(\mu^{(k)},\sigma^{(k)})}, we then draw a sample of size \code{n2} in the posterior distribution of the reproduction number over each time window, conditionnally on this serial interval distribution. 
  After pooling, a sample of size \eqn{\code{n1}\times\code{n2}} of the joint posterior distribution of the reproduction number over each time window is obtained.
  The posterior mean, standard deviation, and 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975 quantiles of the reproduction number for each time window are obtained from this sample.
   
  If \code{plot} is \code{TRUE}, 4 plots are produced.
  The first one shows the epidemic curve. 
  The second one shows the posterior median and 95\% credible interval of the reproduction number. The estimate for a time window is plotted at the end of the time window. 
  The position of the legend on that graph can be monitored by the argument \code{leg.pos} (default is "\code{topright}").
  The third and fourth plots show histograms of the sampled means and standard deviations of the serial interval. 
  
  ----------------------- \code{method "NonParametricUncertainSI"} -----------------------
  
  Details to come...
  
}
\author{ Anne Cori \email{a.cori@imperial.ac.uk} }
\examples{
## load data on pandemic flu in a school in 2009
data("Flu2009")

## estimate the reproduction number (method "NonParametricSI")
EstimateR(Flu2009$Incidence, T.Start=2:26, T.End=8:32, method="NonParametricSI", 
    SI.Distr=Flu2009$SI.Distr, plot=TRUE, leg.pos=xy.coords(1,3))
# the second plot produced shows, at each each day, 
# the estimate of the reproduction number over the 7-day window finishing on that day.

## estimate the reproduction number (method "ParametricSI")
EstimateR(Flu2009$Incidence, T.Start=2:26, T.End=8:32, method="ParametricSI", 
    Mean.SI=2.6, Std.SI=1.5, plot=TRUE)
# the second plot produced shows, at each each day, 
# the estimate of the reproduction number over the 7-day window finishing on that day.

## estimate the reproduction number (method "UncertainSI")
EstimateR(Flu2009$Incidence, T.Start=2:26, T.End=8:32, method="UncertainSI", 
    Mean.SI=2.6, Std.Mean.SI=1, Min.Mean.SI=1, Max.Mean.SI=4.2, 
    Std.SI=1.5, Std.Std.SI=0.5, Min.Std.SI=0.5, Max.Std.SI=2.5, 
    n1=100, n2=100, plot=TRUE)
# the bottom left plot produced shows, at each each day, 
# the estimate of the reproduction number over the 7-day window finishing on that day.

## estimate the reproduction number (method "NonParametricUncertainSI")
## Example to come

}
\keyword{distribution}
