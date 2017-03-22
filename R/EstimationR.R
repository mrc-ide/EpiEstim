#########################################################################################################################
# EstimateR is a wraper which replace the old EstimateR with EstimateR_func and accepts an object of class "cd.fit.mcmc"#
#########################################################################################################################

#' Estimated Instantaneous Reproduction Number
#' 
#' \code{EstimateR} estimates the reproduction number of an epidemic, given the incidence time series and the serial interval distribution. 
#' 
#' @param I One of the following
#' \itemize{
#' \item{A vector (or a dataframe with a single column) of non-negative integers containing the incidence time series}
#' \item{A dataframe of non-negative integers with two columns, so that \code{I$local} contains the incidence of cases due to local transmission and \code{I$imported} contains the incidence of imported cases (with \code{I$local + I$imported} the total incidence).}
#' } 
#' Note that the cases from the first time step are always all assumed to be imported cases. 
#' @param T.Start Vector of positive integers giving the starting times of each window over which the reproduction number will be estimated. These must be in ascending order, and so that for all \code{i}, \code{T.Start[i]<=T.End[i]}. T.Start[1] should be strictly after the first day with non null incidence.
#' @param T.End Vector of positive integers giving the ending times of each window over which the reproduction number will be estimated. These must be in ascending order, and so that for all \code{i}, \code{T.Start[i]<=T.End[i]}. 
#' @param method One of "NonParametricSI", "ParametricSI", "UncertainSI", "SIFromData" or "SIFromSample" (see details).
#' @param n1 For method "UncertainSI" and "SIFromData"; positive integer giving the size of the sample of SI distributions to be drawn (see details).
#' @param n2 For methods "UncertainSI", "SIFromData" and "SIFromSample"; positive integer giving the size of the sample drawn from the posterior distribution of R for each serial interval distribution considered (see details). 
#' @param Mean.SI For method "ParametricSI" and "UncertainSI" ; positive real giving the mean serial interval (method "ParametricSI") or the average mean serial interval (method "UncertainSI", see details).
#' @param Std.SI For method "ParametricSI" and "UncertainSI" ; non negative real giving the stadard deviation of the serial interval (method "ParametricSI") or the average standard deviation of the serial interval (method "UncertainSI", see details).
#' @param Std.Mean.SI For method "UncertainSI" ; standard deviation of the distribution from which mean serial intervals are drawn (see details).
#' @param Min.Mean.SI For method "UncertainSI" ; lower bound of the distribution from which mean serial intervals are drawn (see details).
#' @param Max.Mean.SI For method "UncertainSI" ; upper bound of the distribution from which mean serial intervals are drawn (see details).
#' @param Std.Std.SI For method "UncertainSI" ; standard deviation of the distribution from which standard deviations of the serial interval are drawn (see details).
#' @param Min.Std.SI For method "UncertainSI" ; lower bound of the distribution from which standard deviations of the serial interval are drawn (see details).
#' @param Max.Std.SI For method "UncertainSI" ; upper bound of the distribution from which standard deviations of the serial interval are drawn (see details).
#' @param SI.Distr For method "NonParametricSI" ; vector of probabilities giving the discrete distribution of the serial interval, starting with \code{SI.Distr[1]} (probability that the serial interval is zero), which should be zero.
#' @param SI.Data For method "SIFromData" ; the data on dates of symptoms of pairs of infector/infected individuals to be used to estimate the serial interval distribution (see details).
#' @param SI.parametricDistr For method "SIFromData" ; the parametric distribution to use when estimating the serial interval from data on dates of symptoms of pairs of infector/infected individuals (see details). 
#' Should be one of "G" (Gamma), "E" (Erlang), "off1G" (Gamma shifted by 1), "W" (Weibull), or "L" (Lognormal). 
#' @param MCMC.control For method "SIFromData" ; a list containing the following (see details):
#' \describe{
#'   \item{init.pars}{vector of size 2 corresponding to the initial values of parameters to use for the SI distribution. This is the shape and scale for all but the lognormal distribution, for which it is the meanlog and sdlog. If not specified these are chosen automatically using function \code{\link{init_MCMC_params}}.}
#'   \item{burnin}{a positive integer giving the burnin used in the MCMC when estimating the serial interval distribution.}
#'   \item{thin}{a positive integer corresponding to thinning parameter; the MCMC will be run for \code{burnin+n1*thin iterations}; 1 in \code{thin} iterations will be recorded, after the burnin phase, so the posterior sample size is n1.}
#' }
#' @param SI.Sample For method "SIFromSample" ; a matrix where each column gives one distribution of the serial interval to be explored (see details).
#' @param Mean.Prior A positive number giving the mean of the common prior distribution for all reproduction numbers (see details).
#' @param Std.Prior A positive number giving the standard deviation of the common prior distribution for all reproduction numbers (see details).
#' @param CV.Posterior A positive number giving the aimed posterior coefficient of variation (see details).
#' @param plot Logical. If \code{TRUE} (default is \code{FALSE}), output is plotted (see value).
#' @param leg.pos One of "\code{bottomright}", "\code{bottom}", "\code{bottomleft}", "\code{left}", "\code{topleft}", "\code{top}", "\code{topright}", "\code{right}", "\code{center}" or \code{\link{xy.coords}(x, y)}, with \code{x} and \code{y} real numbers. 
#  This specifies the position of the legend in the plot. Alternatively, \code{locator(1)} can be used ; the user will then need to click where the legend needs to be written.
#' @return {
#' a list with components: 
#' \itemize{
#' \item{R}{: a dataframe containing: 
#' the times of start and end of each time window considered ; 
#' the posterior mean, std, and 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975 quantiles of the reproduction number for each time window.}
#' 	\item{method}{: the method used to estimate R, one of "NonParametricSI", "ParametricSI", "UncertainSI", "SIFromData" or "SIFromSample"}
#' 	\item{SI.Distr}{: a vector or dataframe (depending on the method) containing the discrete serial interval distribution(s) used for estimation}
#' 	\item{SI.Moments}{: a vector or dataframe (depending on the method) containing the mean and std of the discrete serial interval distribution(s) used for estimation}
#' 	\item{I}{: the time series of total incidence}
#' 	\item{I_local}{: the time series of incidence of local cases (so that \code{I_local + I_imported = I})}
#' 	\item{I_imported}{: the time series of incidence of imported cases (so that \code{I_local + I_imported = I})}
#' 	\item{MCMC_converged}{ (only for method \code{SIFromData}): a boolean showing whether the Gelman-Rubin MCMC convergence diagnostic was successful (\code{TRUE}) or not (\code{FALSE})}
#' }
#' }
#' @details{
#' Analytical estimates of the reproduction number for an epidemic over predefined time windows can be obtained within a Bayesian framework, 
#' for a given discrete distribution of the serial interval (see references). 
#' 
#' The more incident cases are observed over a time window, the smallest the posterior coefficient of variation (CV, ratio of standard deviation over mean) of the reproduction number. 
#' An aimed CV can be specified in the argument \code{CV.Posterior} (default is \code{0.3}), and a warning will be produced if the incidence within one of the time windows considered is too low to get this CV. 
#' 
#' The methods vary in the way the serial interval distribution is specified. The plots are also different according to the method used.
#' 
#' In short there are five methods to specify the serial interval distribution (see below for more detail on each method). 
#' In the first two methods, a unique serial interval distribution is considered, whereas in the last three, a range of serial interval distributions are integrated over:
#' \itemize{
#' \item{In method "NonParametricSI" the user specifies the discrete distribution of the serial interval}
#' \item{In method "ParametricSI" the user specifies the mean and sd of the serial interval}
#' \item{In method "UncertainSI" the mean and sd of the serial interval are each drawn from truncated normal distributions, with parameters specified by the user}
#' \item{In method "SIFromData", the serial interval distribution is directly estimated, using MCMC, from interval censored exposure data, with data provided by the user together with a choice of parametric distribution for the serial interval}
#' \item{In method "SIFromSample", the user directly provides the sample of serial interval distribution to use for estimation of R. This can be a useful alternative to the previous method, where the MCMC estimation of the serial interval distribution could be run once, and the same estimated SI distribution then used in EstimateR in different contexts, e.g. with different time windows, hence avoiding to rerun the MCMC everytime EstimateR is called.}
#' }
#' 
#' ----------------------- \code{method "NonParametricSI"} -----------------------
#'   
#' The discrete distribution of the serial interval is directly specified in the argument \code{SI.Distr}.
#' 
#' If \code{plot} is \code{TRUE}, 3 plots are produced. 
#' The first one shows the epidemic curve. 
#' The second one shows the posterior median and 95\% credible interval of the reproduction number. The estimate for a time window is plotted at the end of the time window. 
#' The position of the legend on that graph can be monitored by the argument \code{leg.pos} (default is "\code{topright}").
#' The third plot shows the discrete distribution of the serial interval. 
#' 
#' ----------------------- \code{method "ParametricSI"} -----------------------
#'   
#' The mean and standard deviation of the continuous distribution of the serial interval are given in the arguments \code{Mean.SI} and \code{Std.SI}.
#' The discrete distribution of the serial interval is derived automatically using \code{\link{DiscrSI}}.
#' 
#' If \code{plot} is \code{TRUE}, 3 plots are produced, which are identical to the ones for \code{method "NonParametricSI"} .
#' 
#' ----------------------- \code{method "UncertainSI"} -----------------------
#'    
#' \code{Method "UncertainSI"} allows accounting for uncertainty on the serial interval distribution as described in Cori et al. AJE 2013.
#' We allow the mean \eqn{\mu} and standard deviation \eqn{\sigma} of the serial interval to vary according to truncated normal distributions. 
#' We sample \code{n1} pairs of mean and standard deviations, \eqn{(\mu^{(1)},\sigma^{(1)}),...,(\mu^{(n_2)},\sigma^{(n_2)})}, by first sampling the mean \eqn{\mu^{(k)}} 
#' from its truncated normal distribution (with mean \code{Mean.SI}, standard deviation \code{Std.Mean.SI}, minimum \code{Min.Mean.SI} and maximum \code{Max.Mean.SI}), 
#' and then sampling the standard deviation \eqn{\sigma^{(k)}} from its truncated normal distribution 
#' (with mean \code{Std.SI}, standard deviation \code{Std.Std.SI}, minimum \code{Min.Std.SI} and maximum \code{Max.Std.SI}), but imposing that \eqn{\sigma^{(k)}<\mu^{(k)}}. 
#' This constraint ensures that the Gamma probability density function of the serial interval is null at \eqn{t=0}. 
#' Warnings are produced when the truncated normal distributions are not symmetric around the mean. 
#' For each pair \eqn{(\mu^{(k)},\sigma^{(k)})}, we then draw a sample of size \code{n2} in the posterior distribution of the reproduction number over each time window, conditionnally on this serial interval distribution. 
#' After pooling, a sample of size \eqn{\code{n1}\times\code{n2}} of the joint posterior distribution of the reproduction number over each time window is obtained.
#' The posterior mean, standard deviation, and 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975 quantiles of the reproduction number for each time window are obtained from this sample.
#' 
#' If \code{plot} is \code{TRUE}, 4 plots are produced.
#' The first one shows the epidemic curve. 
#' The second one shows the posterior median and 95\% credible interval of the reproduction number. The estimate for a time window is plotted at the end of the time window. 
#' The position of the legend on that graph can be monitored by the argument \code{leg.pos} (default is "\code{topright}").
#' The third and fourth plots show histograms of the sampled means and standard deviations of the serial interval. 
#' 
#' ----------------------- \code{method "SIFromData"} -----------------------
#'   
#' \code{Method "SIFromData"} allows accounting for uncertainty on the serial interval distribution. 
#' Unlike method "UncertainSI", where we arbitrarily vary the mean and std of the SI in truncated normal distributions, 
#' here, the scope of serial interval distributions considered is directly informed by data
#' on the (potentially censored) dates of symptoms of pairs of infector/infected individuals. 
#' This data, specified in argument \code{SI.Data}, should be a dataframe with 5 columns:
#' \itemize{
#' \item{EL: the lower bound of the symptom onset date of the infector (given as an integer)}
#' \item{ER: the upper bound of the symptom onset date of the infector (given as an integer). Should be such that ER>=EL}
#' \item{SL: the lower bound of the symptom onset date of the infected indivdiual (given as an integer)}
#' \item{SR: the upper bound of the symptom onset date of the infected indivdiual (given as an integer). Should be such that SR>=SL}
#' \item{type (optional): can have entries 0, 1, or 2, corresponding to doubly interval-censored, single interval-censored or exact observations, respsectively, see Reich et al. Statist. Med. 2009. If not specified, this will be automatically computed from the dates}
#' }
#' Assuming a given parametric distribution for the serial interval distribution (specified in SI.parametricDistr), 
#' the posterior distribution of the serial interval is estimated directly fom these data using MCMC methods implemented in the package \code{coarsedatatools}. 
#' The argument \code{MCMC.control} is a list of characteristics which control the MCMC. 
#' The MCMC is run for a total number of iterations of \code{MCMC.control$burnin + n1*MCMC.control$thin};
#' but the output is only recorded after the burnin, and only 1 in every \code{MCMC.control$thin} iterations, 
#' so that the posterior sample size is \code{n1}.
#' For each element in the posterior sample of serial interval distribution, 
#' we then draw a sample of size \code{n2} in the posterior distribution of the reproduction number over each time window, 
#' conditionnally on this serial interval distribution. 
#' After pooling, a sample of size \eqn{\code{n1}\times\code{n2}} of the joint posterior distribution of 
#' the reproduction number over each time window is obtained.
#' The posterior mean, standard deviation, and 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975 quantiles of the reproduction number for each time window are obtained from this sample.
#' 
#' If \code{plot} is \code{TRUE}, 4 plots are produced.
#' The first one shows the epidemic curve. 
#' The second one shows the posterior median and 95\% credible interval of the reproduction number. The estimate for a time window is plotted at the end of the time window. 
#' The position of the legend on that graph can be monitored by the argument \code{leg.pos} (default is "\code{topright}").
#' The third and fourth plots show histograms of the sampled means and standard deviations of the serial interval. 
#' 
#' #' ----------------------- \code{method "SIFromSample"} -----------------------
#' \code{Method "SIFromSample"} also allows accounting for uncertainty on the serial interval distribution. 
#' Unlike methods "UncertainSI" and "SIFromData", the user directly provides (in argument \code{SI.Sample}) a sample of serial interval distribution to be explored. 
#' }
#' @seealso \code{\link{DiscrSI}}
#' @author Anne Cori \email{a.cori@imperial.ac.uk} 
#' @references {
#' Cori, A. et al. A new framework and software to estimate time-varying reproduction numbers during epidemics (AJE 2013).
#' Wallinga, J. and P. Teunis. Different epidemic curves for severe acute respiratory syndrome reveal similar impacts of control measures (AJE 2004).
#' Reich, N.G. et al. Estimating incubation period distributions with coarse data (Statis. Med. 2009)
#' }
#' @importFrom coarseDataTools dic.fit.mcmc
#' @importFrom coda as.mcmc.list as.mcmc
#' @export
#' @examples
#' ## load data on pandemic flu in a school in 2009
#' data("Flu2009")
#' 
#' ## estimate the reproduction number (method "NonParametricSI")
#' EstimateR(Flu2009$Incidence, T.Start=2:26, T.End=8:32, method="NonParametricSI", 
#'           SI.Distr=Flu2009$SI.Distr, plot=TRUE, leg.pos=xy.coords(1,3))
#' # the second plot produced shows, at each each day, 
#' # the estimate of the reproduction number over the 7-day window finishing on that day.
#' 
#' ## estimate the reproduction number (method "ParametricSI")
#' EstimateR(Flu2009$Incidence, T.Start=2:26, T.End=8:32, method="ParametricSI", 
#'           Mean.SI=2.6, Std.SI=1.5, plot=TRUE)
#' # the second plot produced shows, at each each day, 
#' # the estimate of the reproduction number over the 7-day window finishing on that day.
#' 
#' ## estimate the reproduction number (method "UncertainSI")
#' EstimateR(Flu2009$Incidence, T.Start=2:26, T.End=8:32, method="UncertainSI", 
#'           Mean.SI=2.6, Std.Mean.SI=1, Min.Mean.SI=1, Max.Mean.SI=4.2, 
#'           Std.SI=1.5, Std.Std.SI=0.5, Min.Std.SI=0.5, Max.Std.SI=2.5, 
#'           n1=100, n2=100, plot=TRUE)
#' # the bottom left plot produced shows, at each each day, 
#' # the estimate of the reproduction number over the 7-day window finishing on that day.
#' 
#' \dontrun{
#' ## Note the following examples use an MCMC routine 
#' ## to estimate the serial interval distribution from data, 
#' ## so they may take a few minutes to run
#' 
#' ## load data on rotavirus
#' data("MockRotavirus")
#' 
#' ## estimate the reproduction number (method "SIFromData")
#' set.seed(1)
#' R_SIFromData <- EstimateR(MockRotavirus$Incidence, 
#'                                         T.Start=2:47, T.End=8:53, 
#'                                         method="SIFromData", 
#'                                         SI.Data=MockRotavirus$SI.Data, 
#'                                         SI.parametricDistr = "G", 
#'                                         MCMC.control = list(burnin = 1000, thin=10), 
#'                                         n1 = 500, n2 = 50,
#'                                         plot=TRUE, leg.pos=xy.coords(1,3))
#' ## compare with version with no uncertainty
#' R_Parametric <- EstimateR(MockRotavirus$Incidence, 
#'                           T.Start=2:47, T.End=8:53, 
#'                           method="ParametricSI", 
#'                           Mean.SI = mean(R_SIFromData$SI.Moments$Mean), 
#'                           Std.SI = mean(R_SIFromData$SI.Moments$Std), 
#'                           plot=TRUE)
#' ## generate plots
#' p_uncertainty <- plots(R_SIFromData, "R", ylim=c(0, 1.5))
#' p_no_uncertainty <- plots(R_Parametric, "R", ylim=c(0, 1.5))
#' gridExtra::grid.arrange(p_uncertainty, p_no_uncertainty,ncol=2)
#' # the left hand side graph is with uncertainty in the SI distribution, the right hand side without. 
#' # The credible intervals are wider when accounting for uncertainty in the SI distribution. 
#' 
#' #' ## estimate the reproduction number (method "SIFromSample")
#' set.seed(1)
#' SI.fit <- coarseDataTools::dic.fit.mcmc(dat = MockRotavirus$SI.Data, 
#'                              dist="G", 
#'                              init.pars=init_MCMC_params(MockRotavirus$SI.Data, "G"),
#'                              burnin = 1000, 
#'                              n.samples = 5000)
#' SI.Sample <- coarse2estim(SI.fit, thin=10)$SI.Sample
#' R_SIFromSample <- EstimateR(MockRotavirus$Incidence, 
#'                             T.Start=2:47, T.End=8:53, 
#'                             method="SIFromSample", SI.Sample=SI.Sample,
#'                             n2 = 50,
#'                             plot=TRUE, leg.pos=xy.coords(1,3))
#' 
#' # check that R_SIFromSample is the same as R_SIFromData 
#' # since they were generated using the same MCMC algorithm to generate the SI sample
#' # (either internally to EpiEstim or externally)
#' all(R_SIFromSample$R$`Mean(R)` == R_SIFromData$R$`Mean(R)`) 
#' }
#' 
EstimateR <- function(I, T.Start, T.End, method = c("NonParametricSI", "ParametricSI",
                                                    "UncertainSI", "SIFromData", "SIFromSample"), 
                      n1 = NULL, n2 = NULL, Mean.SI = NULL, Std.SI = NULL,
                      Std.Mean.SI = NULL, Min.Mean.SI = NULL, Max.Mean.SI = NULL,
                      Std.Std.SI = NULL, Min.Std.SI = NULL, Max.Std.SI = NULL,
                      SI.Distr = NULL, 
                      SI.Data = NULL, SI.parametricDistr = c("G", "E", "off1G", "W", "L"),  
                      MCMC.control = list(init.pars = NULL, burnin = 3000, thin=10), 
                      SI.Sample = NULL, 
                      Mean.Prior = 5, Std.Prior = 5, CV.Posterior = 0.3,
                      plot = FALSE, leg.pos = "topright") {
  
  ### Need to add warnings if method="SIFromData" and CDT is not null 
  
  if (method=="SIFromData") {
    # Warning if the expected set of parameters is not adequate
    SI.Data <- process_SI.Data(SI.Data)
    
    SI.parametricDistr <- match.arg(SI.parametricDistr)
    if (is.null(n1)) {
      stop("method UncertainSI requires to specify the n1 argument.")
    }
    if (is.null(n2)) {
      stop("method UncertainSI requires to specify the n2 argument.")
    }
    if (n2 <= 0 || n2%%1 != 0) {
      stop("method UncertainSI requires a >0 integer value for n2.")
    }
    if (n1 <= 0 || n1%%1 != 0) {
      stop("method UncertainSI requires a >0 integer value for n1.")
    }
    if(is.null(MCMC.control$init.pars)) MCMC.control$init.pars <- init_MCMC_params(SI.Data, SI.parametricDistr)
    if(SI.parametricDistr=="off1G" & any(SI.Data$SR-SI.Data$EL<=1))
    {
      stop("You cannot fit a Gamma distribution with offset 1 to this SI dataset, because for some data points the maximum serial interval is <=1.\nChoose a different distribution")
    }
    
    CDT <- dic.fit.mcmc(dat = SI.Data, dist=SI.parametricDistr, burnin = MCMC.control$burnin, n.samples = n1*MCMC.control$thin, init.pars=MCMC.control$init.pars)
    
    # check convergence of the MCMC and print warning if not converged
    MCMC_conv <- check_CDTsamples_convergence(CDT@samples)
    
    # thin the chain, and turn the two parameters of the SI distribution into a whole discrete distribution
    c2e <- coarse2estim(CDT, thin=MCMC.control$thin)
    
    cat("\n\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\nEstimating the reproduction number for these serial interval estimates...\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
    
    out <- EstimateR_func(I=I, T.Start=T.Start, T.End=T.End, method = "SIFromData", n1=n1 , n2=n2 , Mean.SI=NULL , Std.SI=NULL ,
                   Std.Mean.SI=NULL , Min.Mean.SI=NULL , Max.Mean.SI=NULL ,
                   Std.Std.SI=NULL , Min.Std.SI=NULL , Max.Std.SI=NULL ,
                   SI.Distr=NULL , SI.Sample= c2e$SI.Sample , Mean.Prior=Mean.Prior , Std.Prior=Std.Prior, CV.Posterior=CV.Posterior ,
                   plot=plot , leg.pos=leg.pos)
    out[["MCMC_converged"]] <- MCMC_conv
  } else {
    out <- EstimateR_func(I=I, T.Start=T.Start, T.End=T.End, method = method, n1=n1 , n2=n2 , Mean.SI=Mean.SI , Std.SI=Std.SI ,
                   Std.Mean.SI=Std.Mean.SI , Min.Mean.SI=Min.Mean.SI , Max.Mean.SI=Max.Mean.SI ,
                   Std.Std.SI=Std.Std.SI , Min.Std.SI=Min.Std.SI , Max.Std.SI=Max.Std.SI ,
                   SI.Distr=SI.Distr , SI.Sample= SI.Sample, Mean.Prior=Mean.Prior , Std.Prior=Std.Prior, CV.Posterior=CV.Posterior ,
                   plot=plot , leg.pos=leg.pos)
  }
  return(out)
}

#########################################################
# EstimateR_func: Doing the heavy work in EstimateR     #
#########################################################

#' @import reshape2 grid gridExtra
#' @importFrom ggplot2 last_plot ggplot aes geom_step ggtitle geom_ribbon geom_line xlab ylab xlim geom_hline ylim geom_histogram
#' @importFrom plotly layout mutate arrange rename summarise filter ggplotly
#' @importFrom stats median pgamma plnorm pweibull qgamma qlnorm quantile qweibull rgamma rmultinom rnorm sd
EstimateR_func <- function (I, T.Start, T.End, method = c("NonParametricSI", "ParametricSI",
                                                          "UncertainSI", "SIFromData", "SIFromSample"), 
                            n1 = NULL, n2 = NULL, Mean.SI = NULL, Std.SI = NULL,
                            Std.Mean.SI = NULL, Min.Mean.SI = NULL, Max.Mean.SI = NULL,
                            Std.Std.SI = NULL, Min.Std.SI = NULL, Max.Std.SI = NULL,
                            SI.Distr = NULL, 
                            SI.Sample = NULL, 
                            Mean.Prior = 5, Std.Prior = 5, CV.Posterior = 0.3,
                            plot = FALSE, leg.pos = "topright")
{
  
  #########################################################
  # Calculates the cumulative incidence over time steps   #
  #########################################################
  
  CalculIncidencePerTimeStep <- function(I, T.Start, T.End) {
    NbTimePeriods <- length(T.Start)
    IncidencePerTimeStep <- sapply(1:NbTimePeriods, function(i) sum(I[T.Start[i]:T.End[i],]))
    return(IncidencePerTimeStep)
  }
  
  #########################################################
  # Calculates the parameters of the Gamma posterior      #
  # distribution from the discrete SI distribution        #
  #########################################################
  
  PosteriorFromSIDistr <- function(I, SI.Distr, a.Prior, b.Prior,
                                   T.Start, T.End) {
    NbTimePeriods <- length(T.Start)
    lambda <- OverallInfectivity(I, SI.Distr)
    FinalMean.SI <- sum(SI.Distr * (0:(length(SI.Distr) -
                                         1)))
    a.Posterior <- vector()
    b.Posterior <- vector()
    a.Posterior <- sapply(1:(NbTimePeriods), function(t) if (T.End[t] >
                                                             FinalMean.SI) {
      a.Prior + sum(I[T.Start[t]:T.End[t],])
    }
    else {
      NA
    })
    b.Posterior <- sapply(1:(NbTimePeriods), function(t) if (T.End[t] >
                                                             FinalMean.SI) {
      1/(1/b.Prior + sum(lambda[T.Start[t]:T.End[t]]))
    }
    else {
      NA
    })
    return(list(a.Posterior, b.Posterior))
  }
  
  #########################################################
  # Samples from the Gamma posterior distribution for a   #
  # given mean SI and std SI                              #
  #########################################################
  
  SampleFromPosterior <- function(SampleSize, I, Mean.SI, Std.SI, SI.Distr=NULL, 
                                  a.Prior, b.Prior, T.Start, T.End) {
    NbTimePeriods <- length(T.Start)
    
    if(is.null(SI.Distr))
      SI.Distr <- sapply(1:T, function(t) DiscrSI(t - 1, Mean.SI, Std.SI))
    
    FinalMean.SI <- sum(SI.Distr * (0:(length(SI.Distr) -
                                         1)))
    lambda <- OverallInfectivity(I, SI.Distr)
    a.Posterior <- vector()
    b.Posterior <- vector()
    a.Posterior <- sapply(1:(NbTimePeriods), function(t) if (T.End[t] >
                                                             FinalMean.SI) {
      a.Prior + sum(I[T.Start[t]:T.End[t],])
    }
    else {
      NA
    })
    b.Posterior <- sapply(1:(NbTimePeriods), function(t) if (T.End[t] >
                                                             FinalMean.SI) {
      1/(1/b.Prior + sum(lambda[T.Start[t]:T.End[t]], na.rm = TRUE))
    }
    else {
      NA
    })
    SampleR.Posterior <- sapply(1:(NbTimePeriods), function(t) if (!is.na(a.Posterior[t])) {
      rgamma(SampleSize, shape = unlist(a.Posterior[t]),
             scale = unlist(b.Posterior[t]))
    }
    else {
      rep(NA, SampleSize)
    })
    return(list(SampleR.Posterior, SI.Distr))
  }
  method <- match.arg(method)
  
  I <- process_I(I)
  T<-nrow(I)
  
  if (Mean.Prior <= 0) {
    stop("Mean.Prior must be >0.")
  }
  if (Std.Prior <= 0) {
    stop("Std.Prior must be >0.")
  }
  a.Prior <- (Mean.Prior/Std.Prior)^2
  b.Prior <- Std.Prior^2/Mean.Prior
  if (!is.vector(T.Start)) {
    stop("T.Start must be a vector.")
  }
  if (!is.vector(T.End)) {
    stop("T.End must be a vector.")
  }
  if (length(T.Start) != length(T.End)) {
    stop("T.Start and T.End must have the same length.")
  }
  NbTimePeriods <- length(T.Start)
  if (any(T.Start > T.End)) {
    stop("T.Start[i] must be <= T.End[i] for all i.")
  }
  if (any(T.Start < 1 || T.Start%%1 != 0)) {
    stop("T.Start must be a vector of >0 integers.")
  }
  if (any(T.End < 1 || T.End%%1 != 0)) {
    stop("T.End must be a vector of >0 integers.")
  }
  if (method == "NonParametricSI") {
    if (is.null(SI.Distr)) {
      stop("method NonParametricSI requires to specify the SI.Distr argument.")
    }
    if (!is.vector(SI.Distr)) {
      stop("method NonParametricSI requires that SI.Distr must be a vector.")
    }
    if (SI.Distr[1] != 0) {
      stop("method NonParametricSI requires that SI.Distr[1] = 0.")
    }
    if (any(SI.Distr < 0)) {
      stop("method NonParametricSI requires that SI.Distr must be a positive vector.")
    }
    if (abs(sum(SI.Distr) - 1) > 0.01) {
      stop("method NonParametricSI requires that SI.Distr must sum to 1.")
    }
  }
  if (method == "ParametricSI") {
    if (is.null(Mean.SI)) {
      stop("method NonParametricSI requires to specify the Mean.SI argument.")
    }
    if (is.null(Std.SI)) {
      stop("method NonParametricSI requires to specify the Std.SI argument.")
    }
    if (Mean.SI <= 1) {
      stop("method ParametricSI requires a value >1 for Mean.SI.")
    }
    if (Std.SI <= 0) {
      stop("method ParametricSI requires a >0 value for Std.SI.")
    }
  }
  if (method == "UncertainSI") {
    if (is.null(Mean.SI)) {
      stop("method UncertainSI requires to specify the Mean.SI argument.")
    }
    if (is.null(Std.SI)) {
      stop("method UncertainSI requires to specify the Std.SI argument.")
    }
    if (is.null(n1)) {
      stop("method UncertainSI requires to specify the n1 argument.")
    }
    if (is.null(n2)) {
      stop("method UncertainSI requires to specify the n2 argument.")
    }
    if (is.null(Std.Mean.SI)) {
      stop("method UncertainSI requires to specify the Std.Mean.SI argument.")
    }
    if (is.null(Min.Mean.SI)) {
      stop("method UncertainSI requires to specify the Min.Mean.SI argument.")
    }
    if (is.null(Max.Mean.SI)) {
      stop("method UncertainSI requires to specify the Max.Mean.SI argument.")
    }
    if (is.null(Std.Std.SI)) {
      stop("method UncertainSI requires to specify the Std.Std.SI argument.")
    }
    if (is.null(Min.Std.SI)) {
      stop("method UncertainSI requires to specify the Min.Std.SI argument.")
    }
    if (is.null(Max.Std.SI)) {
      stop("method UncertainSI requires to specify the Max.Std.SI argument.")
    }
    if (Mean.SI <= 0) {
      stop("method UncertainSI requires a >0 value for Mean.SI.")
    }
    if (Std.SI <= 0) {
      stop("method UncertainSI requires a >0 value for Std.SI.")
    }
    if (n2 <= 0 || n2%%1 != 0) {
      stop("method UncertainSI requires a >0 integer value for n2.")
    }
    if (n1 <= 0 || n1%%1 != 0) {
      stop("method UncertainSI requires a >0 integer value for n1.")
    }
    if (Std.Mean.SI <= 0) {
      stop("method UncertainSI requires a >0 value for Std.Mean.SI.")
    }
    if (Min.Mean.SI < 1) {
      stop("method UncertainSI requires a value >=1 for Min.Mean.SI.")
    }
    if (Max.Mean.SI < Mean.SI) {
      stop("method UncertainSI requires that Max.Mean.SI >= Mean.SI.")
    }
    if (Mean.SI < Min.Mean.SI) {
      stop("method UncertainSI requires that Mean.SI >= Min.Mean.SI.")
    }
    if (signif(Max.Mean.SI - Mean.SI, 3) != signif(Mean.SI -
                                                   Min.Mean.SI, 3)) {
      warning("The distribution you chose for the mean SI is not centered around the mean.")
    }
    if (Std.Std.SI <= 0) {
      stop("method UncertainSI requires a >0 value for Std.Std.SI.")
    }
    if (Min.Std.SI <= 0) {
      stop("method UncertainSI requires a >0 value for Min.Std.SI.")
    }
    if (Max.Std.SI < Std.SI) {
      stop("method UncertainSI requires that Max.Std.SI >= Std.SI.")
    }
    if (Std.SI < Min.Std.SI) {
      stop("method UncertainSI requires that Std.SI >= Min.Std.SI.")
    }
    if (signif(Max.Std.SI - Std.SI, 3) != signif(Std.SI -
                                                 Min.Std.SI, 3)) {
      warning("The distribution you chose for the std of the SI is not centered around the mean.")
    }
  }
  if(method == "SIFromSample")
  {
    if (is.null(SI.Sample)) {
      stop("method SIFromSample requires to specify the SI.Sample argument.")
    }
    if (!is.matrix(SI.Sample)) {
      stop("method SIFromSample requires that SI.Sample must be a is.matrix")
    }
    if (any(SI.Sample[1,] != 0)) {
      stop("method SIFromSample requires that SI.Sample[1,] contains only 0.")
    }
    if (any(SI.Sample < 0)) {
      stop("method SIFromSample requires that SI.Sample must contain only non negtaive values.")
    }
    if (any(abs(colSums(SI.Sample) - 1) > 0.01)) {
      stop("method SIFromSample requires the sum of each column in SI.Sample to be 1.")
    }
  }
  if (CV.Posterior < 0) {
    stop("CV.Posterior must be >0.")
  }
  MinNbCasesPerTimePeriod <- ceiling(1/CV.Posterior^2 - a.Prior)
  IncidencePerTimeStep <- CalculIncidencePerTimeStep(I, T.Start,
                                                     T.End)
  if (IncidencePerTimeStep[1] < MinNbCasesPerTimePeriod) {
    warning("You're estimating R too early in the epidemic to get the desired posterior CV.")
  }
  if (plot != TRUE && plot != FALSE) {
    stop("plot must be TRUE or FALSE.")
  }
  if (method == "NonParametricSI") {
    SIUncertainty <- "N"
    ParametricSI <- "N"
  }
  if (method == "ParametricSI") {
    SIUncertainty <- "N"
    ParametricSI <- "Y"
  }
  if (method == "UncertainSI") {
    SIUncertainty <- "Y"
    ParametricSI <- "Y"
  }
  if (method == "SIFromData" | method == "SIFromSample") {
    SIUncertainty <- "Y"
    ParametricSI <- "N"
  }
  if (SIUncertainty == "Y") {
    if  (ParametricSI == "Y") {
      Mean.SI.sample <- rep(-1, n1)
      Std.SI.sample <- rep(-1, n1)
      for (k in 1:n1) {
        while (Mean.SI.sample[k] < Min.Mean.SI || Mean.SI.sample[k] >
               Max.Mean.SI) {
          Mean.SI.sample[k] <- rnorm(1, mean = Mean.SI,
                                     sd = Std.Mean.SI)
        }
        while (Std.SI.sample[k] < Min.Std.SI || Std.SI.sample[k] >
               Max.Std.SI || Std.SI.sample[k] > Mean.SI.sample[k]) {
          Std.SI.sample[k] <- rnorm(1, mean = Std.SI, sd = Std.Std.SI)
        }
      }
      temp <- lapply(1:n1, function(k) SampleFromPosterior(n2,
                                                           I, Mean.SI.sample[k], Std.SI.sample[k], SI.Distr=NULL, a.Prior,
                                                           b.Prior, T.Start, T.End))
      SI.Distr <- cbind(t(sapply(1:n1, function(k) (temp[[k]])[[2]])),
                        rep(0, n1))
      Rsample <- matrix(NA, n2 * n1, NbTimePeriods)
      for (k in 1:n1) {
        Rsample[((k - 1) * n2 + 1):(k * n2), which(T.End >
                                                     Mean.SI.sample[k])] <- (temp[[k]])[[1]][, which(T.End >
                                                                                                       Mean.SI.sample[k])]
      }
      Mean.Posterior <- apply(Rsample, 2, mean, na.rm = TRUE)
      Std.Posterior <- apply(Rsample, 2, sd, na.rm = TRUE)
      Quantile.0.025.Posterior <- apply(Rsample, 2, quantile,
                                        0.025, na.rm = TRUE)
      Quantile.0.05.Posterior <- apply(Rsample, 2, quantile,
                                       0.05, na.rm = TRUE)
      Quantile.0.25.Posterior <- apply(Rsample, 2, quantile,
                                       0.25, na.rm = TRUE)
      Median.Posterior <- apply(Rsample, 2, median, na.rm = TRUE)
      Quantile.0.75.Posterior <- apply(Rsample, 2, quantile,
                                       0.75, na.rm = TRUE)
      Quantile.0.95.Posterior <- apply(Rsample, 2, quantile,
                                       0.95, na.rm = TRUE)
      Quantile.0.975.Posterior <- apply(Rsample, 2, quantile,
                                        0.975, na.rm = TRUE)
    }
    else {
      n1<-dim(SI.Sample)[2]
      Mean.SI.sample <- rep(-1, n1)
      Std.SI.sample <- rep(-1, n1)
      for (k in 1:n1) {
        Mean.SI.sample[k] <- sum((1:dim(SI.Sample)[1]-1)*SI.Sample[,k])
        Std.SI.sample[k] <- sqrt(sum(SI.Sample[,k]*((1:dim(SI.Sample)[1]-1) - Mean.SI.sample[k])^2))
      }
      temp <- lapply(1:n1, function(k) SampleFromPosterior(n2,
                                                           I, Mean.SI=NULL, Std.SI=NULL, SI.Sample[,k], a.Prior,
                                                           b.Prior, T.Start, T.End))
      SI.Distr <- cbind(t(sapply(1:n1, function(k) (temp[[k]])[[2]])),
                        rep(0, n1))
      Rsample <- matrix(NA, n2 * n1, NbTimePeriods)
      for (k in 1:n1) {
        Rsample[((k - 1) * n2 + 1):(k * n2), which(T.End >
                                                     Mean.SI.sample[k])] <- (temp[[k]])[[1]][, which(T.End >
                                                                                                       Mean.SI.sample[k])]
      }
      Mean.Posterior <- apply(Rsample, 2, mean, na.rm = TRUE)
      Std.Posterior <- apply(Rsample, 2, sd, na.rm = TRUE)
      Quantile.0.025.Posterior <- apply(Rsample, 2, quantile,
                                        0.025, na.rm = TRUE)
      Quantile.0.05.Posterior <- apply(Rsample, 2, quantile,
                                       0.05, na.rm = TRUE)
      Quantile.0.25.Posterior <- apply(Rsample, 2, quantile,
                                       0.25, na.rm = TRUE)
      Median.Posterior <- apply(Rsample, 2, median, na.rm = TRUE)
      Quantile.0.75.Posterior <- apply(Rsample, 2, quantile,
                                       0.75, na.rm = TRUE)
      Quantile.0.95.Posterior <- apply(Rsample, 2, quantile,
                                       0.95, na.rm = TRUE)
      Quantile.0.975.Posterior <- apply(Rsample, 2, quantile,
                                        0.975, na.rm = TRUE)
    }
  }else{
    # CertainSI
    if (ParametricSI == "Y") {
      SI.Distr <- sapply(1:T, function(t) DiscrSI(t - 1,
                                                  Mean.SI, Std.SI))
    }
    if (length(SI.Distr) < T + 1) {
      SI.Distr[(length(SI.Distr) + 1):(T + 1)] <- 0
    }
    FinalMean.SI <- sum(SI.Distr * (0:(length(SI.Distr) -
                                         1)))
    FinalStd.SI <- sqrt(sum(SI.Distr * (0:(length(SI.Distr) -
                                             1))^2) - FinalMean.SI^2)
    post <- PosteriorFromSIDistr(I, SI.Distr, a.Prior, b.Prior,
                                 T.Start, T.End)
    a.Posterior <- unlist(post[[1]])
    b.Posterior <- unlist(post[[2]])
    Mean.Posterior <- a.Posterior * b.Posterior
    Std.Posterior <- sqrt(a.Posterior) * b.Posterior
    Quantile.0.025.Posterior <- qgamma(0.025, shape = a.Posterior,
                                       scale = b.Posterior, lower.tail = TRUE, log.p = FALSE)
    Quantile.0.05.Posterior <- qgamma(0.05, shape = a.Posterior,
                                      scale = b.Posterior, lower.tail = TRUE, log.p = FALSE)
    Quantile.0.25.Posterior <- qgamma(0.25, shape = a.Posterior,
                                      scale = b.Posterior, lower.tail = TRUE, log.p = FALSE)
    Median.Posterior <- qgamma(0.5, shape = a.Posterior,
                               scale = b.Posterior, lower.tail = TRUE, log.p = FALSE)
    Quantile.0.75.Posterior <- qgamma(0.75, shape = a.Posterior,
                                      scale = b.Posterior, lower.tail = TRUE, log.p = FALSE)
    Quantile.0.95.Posterior <- qgamma(0.95, shape = a.Posterior,
                                      scale = b.Posterior, lower.tail = TRUE, log.p = FALSE)
    Quantile.0.975.Posterior <- qgamma(0.975, shape = a.Posterior,
                                       scale = b.Posterior, lower.tail = TRUE, log.p = FALSE)
  }
  results <- list()
  results$R <- as.data.frame(cbind(T.Start, T.End, Mean.Posterior,
                                   Std.Posterior, Quantile.0.025.Posterior, Quantile.0.05.Posterior,
                                   Quantile.0.25.Posterior, Median.Posterior, Quantile.0.75.Posterior,
                                   Quantile.0.95.Posterior, Quantile.0.975.Posterior))
  names(results$R) <- c("T.Start", "T.End", "Mean(R)", "Std(R)",
                        "Quantile.0.025(R)", "Quantile.0.05(R)", "Quantile.0.25(R)",
                        "Median(R)", "Quantile.0.75(R)", "Quantile.0.95(R)",
                        "Quantile.0.975(R)")
  results$method <- method
  results$SI.Distr <- SI.Distr
  if (SIUncertainty == "Y") {
    results$SI.Moments <- as.data.frame(cbind(Mean.SI.sample,
                                              Std.SI.sample))
    names(results$SI.Moments) <- c("Mean","Std")
  }else {
    results$SI.Moments <- as.data.frame(cbind(FinalMean.SI,
                                              FinalStd.SI))
    names(results$SI.Moments) <- c("Mean", "Std")
  }
  
  results$I <- rowSums(I)
  results$I_local <- I$local
  results$I_imported <- I$imported
  
  if (plot) {
    
    ########################################################################
    ### these few lines are to make CRAN checks happy with ggplot2... ###
    Time <- NULL
    Incidence <- NULL
    value <- NULL
    meanR <- NULL
    group <- NULL
    lower <- NULL
    upper <- NULL
    Times <- NULL
    ..density.. <- NULL
    start <- NULL
    end <- NULL
    ########################################################################
    
    p1 <- ggplot(data.frame(Time=1:T, Incidence=rowSums(I)), aes(x=Time, y=Incidence)) +
      geom_step() +
      ggtitle("Epidemic curve")
    #p1ly <- ggplotly(p1)
    
    # test if intervals overlap 
    time.points <- apply(results$R[,c("T.Start","T.End") ], 1, function(x) x[1]:(x[2]-1)) 
    if (length(time.points) == length(unique(matrix(time.points,ncol=1)))) { 
      
      df <- melt(data.frame(start=T.Start, end=T.End, meanR=Mean.Posterior, lower=Quantile.0.025.Posterior,
                            upper=Quantile.0.975.Posterior), id=c("meanR", "lower", "upper")) 
      df$group <- as.factor(rep(1:length(T.Start), dim(df)[1]/length(T.Start)))
      
      p2 <- ggplot(df, aes(x=as.numeric(value), y=as.numeric(meanR), group=as.factor(group))) +
        geom_ribbon(aes(ymin=lower, ymax=upper), colour=NA, fill="black", alpha=0.2) +
        geom_line() +
        xlab("Time") +
        ylab("R") +
        xlim(c(1,max(T.End))) +
        ggtitle("Estimated R")
      #p2ly <- ggplotly(p2)
      
    } else { 
      
      p2 <- ggplot(data.frame(start=T.Start, end=T.End, meanR=Mean.Posterior, lower=Quantile.0.025.Posterior,
                              upper=Quantile.0.975.Posterior), aes(end, meanR)) +
        geom_ribbon(aes(ymin=lower, ymax=upper, fill="95%CrI")) +
        geom_line(aes(colour="Mean")) +
        geom_hline(yintercept=1, linetype="dotted") +
        xlab("Time") +
        ylab("R") +
        xlim(c(1,max(T.End))) +
        ylim(c(0,max(Quantile.0.975.Posterior, na.rm = TRUE))) +
        ggtitle("Estimated R") +
        scale_colour_manual("",values="black")+
        scale_fill_manual("",values="grey")
      #p2ly <- ggplotly(p2) 
    }
    
    if (SIUncertainty == "Y") {
      p3 <- ggplot(data.frame(Mean.SI.sample), aes(Mean.SI.sample)) +
        geom_histogram(bins=30) +
        xlab("Mean serial interval") +
        ylab("Density") + 
        ggtitle("Explored \n mean serial intervals")
      
      p4 <- ggplot(data.frame(Std.SI.sample), aes(Std.SI.sample)) +
        geom_histogram(bins=30) +
        xlab("Std serial interval") +
        ggtitle("Explored \n std serial intervals")
      
      grid.arrange(p1,p2,p3,p4,ncol=2)
      
    } else {
      SI.Distr.times <- unlist(apply(data.frame(0:(length(SI.Distr) - 1), SI.Distr), 1,
                                     function(x) {if (x[2]!=0) unlist(rep(x[1],round(x[2]*1000)),use.names=FALSE)}))
      names(SI.Distr.times) <- NULL
      
      p3 <- ggplot(data.frame(Times=SI.Distr.times), aes(0.5+Times)) +
        geom_histogram(binwidth=1, aes(y=..density..)) +
        xlab("Time") + 
        xlim(c(0,0.5+max(SI.Distr.times))) + 
        ylab("Frequency") + 
        ggtitle("Serial interval distribution") 
      #p3ly <- ggplotly(p3)
      
      grid.arrange(p1,p3,p2,ncol=1)
      
    }
  }
  
  return(results)
}

