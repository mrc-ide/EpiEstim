#########################################################################################################################
# WT function to estimate Rc the case reproduction number #
#########################################################################################################################

#' Estimation of the case reproduction number using the Wallinga and Teunis method
#' 
#' \code{WT} estimates the case reproduction number of an epidemic, given the incidence time series and the serial interval distribution. 
#' 
#' @param I One of the following
#' \itemize{
#' \item{Vector (or a dataframe with a column named 'I') of non-negative integers containing an incidence time series.
#'  If the dataframe contains a column \code{I$dates}, this is used for plotting. 
#'  \code{I$dates} must contains only dates in a row.}
#' \item{An object of class \code{\link[incidence]{incidence}}}
#'  }
#' @param t_start Vector of positive integers giving the starting times of each window over which the reproduction number will be estimated. These must be in ascending order, and so that for all \code{i}, \code{t_start[i]<=t_end[i]}. t_start[1] should be strictly after the first day with non null incidence.
#' @param t_end Vector of positive integers giving the ending times of each window over which the reproduction number will be estimated. These must be in ascending order, and so that for all \code{i}, \code{t_start[i]<=t_end[i]}.
#' @param method One of "non_parametric_si" or "parametric_si" (see details).
#' @param mean_si For method "parametric_si" ; positive real giving the mean serial interval.
#' @param std_si For method "parametric_si" ; non negative real giving the stadard deviation of the serial interval.
#' @param si_distr For method "non_parametric_si" ; vector of probabilities giving the discrete distribution of the serial interval, starting with \code{si_distr[1]} (probability that the serial interval is zero), which should be zero.
#' @param nSim A positive integer giving the number of simulated epidemic trees used for computation of the confidence intervals of the case reproduction number (see details).
#' @param plot Logical. If \code{TRUE} (default is \code{FALSE}), output is plotted (see value).
#' @param legend A boolean (TRUE by default) governing the presence / absence of legends on the plots
#' This specifies the position of the legend in the plot. Alternatively, \code{locator(1)} can be used ; the user will then need to click where the legend needs to be written.
#' @return {
#' 	a list with components: 
#' 	\itemize{
#' 	\item{R}{: a dataframe containing: 
#' 	    the times of start and end of each time window considered ; 
#' 	    the estimated mean, std, and 0.025 and 0.975 quantiles of the reproduction number for each time window.}
#' 	\item{method}{: the method used to estimate R, one of "non_parametric_si", "parametric_si", "uncertain_si", "si_from_data" or "si_from_sample"}
#' 	\item{si_distr}{: a vector containing the discrete serial interval distribution used for estimation}
#' 	\item{SI.Moments}{: a vector containing the mean and std of the discrete serial interval distribution(s) used for estimation}
#' 	\item{I}{: the time series of total incidence}
#' 	\item{I_local}{: the time series of incidence of local cases (so that \code{I_local + I_imported = I})}
#' 	\item{I_imported}{: the time series of incidence of imported cases (so that \code{I_local + I_imported = I})}
#' 	\item{dates}{: a vector of dates corresponding to the incidence time series}
#' 	}
#' 	}
#' @details{
#' Estimates of the case reproduction number for an epidemic over predefined time windows can be obtained, 
#' for a given discrete distribution of the serial interval, as proposed by Wallinga and Teunis (AJE, 2004). 
#' Confidence intervals are obtained by simulating a number (nSim) of possible transmission trees. 
#' 
#' The methods vary in the way the serial interval distribution is specified.
#' 
#' ----------------------- \code{method "non_parametric_si"} -----------------------
#' 
#' The discrete distribution of the serial interval is directly specified in the argument \code{si_distr}.
#' 
#' If \code{plot} is \code{TRUE}, 3 plots are produced. 
#' The first one shows the epidemic curve. 
#' The second one shows the posterior mean and 95\% credible interval of the reproduction number. The estimate for a time window is plotted at the end of the time window. 
#' The third plot shows the discrete distribution of the serial interval. 
#'
#' ----------------------- \code{method "parametric_si"} -----------------------
#' 
#' The mean and standard deviation of the continuous distribution of the serial interval are given in the arguments \code{mean_si} and \code{std_si}.
#' The discrete distribution of the serial interval is derived automatically using \code{\link{discr_si}}.
#'
#' If \code{plot} is \code{TRUE}, 3 plots are produced, which are identical to the ones for \code{method "non_parametric_si"} .
#' }
#' @seealso \code{\link{discr_si}}, \code{\link{estimate_r}}
#' @author Anne Cori \email{a.cori@imperial.ac.uk} 
#' @references {
#' Cori, A. et al. A new framework and software to estimate time-varying reproduction numbers during epidemics (AJE 2013).
#' Wallinga, J. and P. Teunis. Different epidemic curves for severe acute respiratory syndrome reveal similar impacts of control measures (AJE 2004).
#' }
#' @export
#' @import reshape2 grid gridExtra
#' @importFrom ggplot2 last_plot ggplot aes geom_step ggtitle geom_ribbon geom_line xlab ylab xlim geom_hline ylim geom_histogram
#' @importFrom plotly layout mutate arrange rename summarise filter ggplotly
#' @examples
#' ## load data on pandemic flu in a school in 2009
#' data("Flu2009")
#' 
#' ## estimate the case reproduction number (method "non_parametric_si")
#' WT(Flu2009$incidence, t_start=2:26, t_end=8:32, method="non_parametric_si", 
#'    si_distr=Flu2009$si_distr, plot=TRUE, nSim=100)
#' # the second plot produced shows, at each each day, 
#' # the estimate of the case reproduction number over the 7-day window finishing on that day.
#' 
#' ## estimate the case reproduction number (method "parametric_si")
#' WT(Flu2009$incidence, t_start=2:26, t_end=8:32, method="parametric_si", 
#'    mean_si=2.6, std_si=1.5, plot=TRUE, nSim=100)
#' # the second plot produced shows, at each each day, 
#' # the estimate of the case reproduction number over the 7-day window finishing on that day.
WT <- function(I, t_start, t_end,
               method=c("non_parametric_si","parametric_si"),
               mean_si=NULL, std_si=NULL, si_distr=NULL, nSim=10, 
               plot=FALSE, legend=FALSE)
{
  
  ### Functions ###
  
  #########################################################
  # Draws a possile transmission tree                     #
  #########################################################
  
  DrawOneSetOfAncestries <- function()
  {
    res <- vector()
    for(t in t_start[1]:t_end[length(t_end)])
    {
      if(length(which(Onset==t))>0)
      {
        if(length(possibleAncesTime[[t]])>0)
        {
          prob <- si_distr[t-possibleAncesTime[[t]]+1]*I[possibleAncesTime[[t]]]
          if(any(prob>0))
          {
            res[which(Onset==t)] <- possibleAncesTime[[t]][which(rmultinom(length(which(Onset==t)),size=1,prob=prob)==TRUE,arr.ind=TRUE)[,1]]
          }else
          {
            res[which(Onset==t)] <- NA
          }
        }else
          res[which(Onset==t)] <- NA
      }
    }
    return(res)
  }
  
  ### Error messages ###
  
  method <- match.arg(method)
  
  I <- process_I(I)
  if(!is.null(I$dates)) 
  {
    dates <- check_dates(I)
    I <- process_I_vector(rowSums(I[,c("local","imported")]))
    T<-length(I)
  }else
  {
    I <- process_I_vector(rowSums(I[,c("local","imported")]))
    T<-length(I)
    dates <- 1:T
  }
  
  
  ### Adjusting t_start and t_end so that at least an incident case has been observed before t_start[1] ###
  
  i <- 1
  while(sum(I[1:(t_start[i]-1)])==0)
  {
    i <- i+1
  }
  temp <- which(t_start<i)
  if(length(temp>0))
  {
    t_start <- t_start[-temp]
    t_end <- t_end[-temp]
  }
  
  check_times(t_start, t_end, T)
  nb_time_periods <- length(t_start)
  
  if(method=="non_parametric_si")
  {
    check_si_distr(si_distr)
    si_distr <- c(si_distr,0)
  }
  
  if(method=="parametric_si")
  {
    if(is.null(mean_si))
    {
      stop("method non_parametric_si requires to specify the mean_si argument.")
    }
    if(is.null(std_si))
    {
      stop("method non_parametric_si requires to specify the std_si argument.")
    }
    if(mean_si<1)
    {
      stop("method parametric_si requires a value >1 for mean_si.")
    }
    if(std_si<0)
    {
      stop("method parametric_si requires a >0 value for std_si.")
    }
  }
  
  if(plot!=TRUE && plot!=FALSE)
  {
    stop("plot must be TRUE or FALSE.")
  }
  
  if(!is.numeric(nSim))
  {
    stop("nSim must be a positive integer.")
  }
  if(nSim<=0)
  {
    stop("nSim must be a positive integer.")
  }
  
  ### What does each method do ###
  
  if(method=="non_parametric_si")
  {
    parametric_si<-"N"
  }
  if(method=="parametric_si")
  {
    parametric_si<-"Y"
  }
  
  if(parametric_si=="Y")
  {
    si_distr <- discr_si(0:(T-1),mean_si,std_si)
  }
  if(length(si_distr)<T+1){si_distr[(length(si_distr)+1):(T+1)]<-0}
  final_mean_si<-sum(si_distr*(0:(length(si_distr)-1)))
  Finalstd_si<-sqrt(sum(si_distr*(0:(length(si_distr)-1))^2)-final_mean_si^2)
  
  TimePeriodsWithNoincidence <- vector()
  for(i in 1:nb_time_periods)
  {
    if(sum(I[t_start[i]:t_end[i]])==0)
    {
      TimePeriodsWithNoincidence <- c(TimePeriodsWithNoincidence,i)
    }
  }
  if(length(TimePeriodsWithNoincidence)>0)
  {
    t_start <- t_start[-TimePeriodsWithNoincidence]
    t_end <- t_end[-TimePeriodsWithNoincidence]
    nb_time_periods <- length(t_start)
  }
  
  Onset <- vector()
  for (t in 1:T) {
    Onset <- c(Onset, rep(t, I[t]))
  }
  NbCases <- length(Onset)
  
  delay <- outer (1:T, 1:T, "-")
  SIdelay <- apply(delay, 2, function(x) si_distr[pmin(pmax(x + 1, 1), length(si_distr))])
  SumOnColSIDelay_tmp <- sapply(1:nrow(SIdelay), function(i) sum(SIdelay [i,]* I, na.rm=TRUE) )
  SumOnColSIDelay <- vector()
  for (t in 1:T) {
    SumOnColSIDelay <- c(SumOnColSIDelay, rep(SumOnColSIDelay_tmp[t], I[t]))
  }
  MatSumOnColSIDelay <- matrix(rep(SumOnColSIDelay_tmp, T), nrow = T, ncol = T)
  p <- SIdelay/(MatSumOnColSIDelay)
  p[which(is.na(p))] <- 0
  p[which(is.infinite(p))] <- 0
  MeanRperIndexCaseDate <- sapply(1:ncol(p), function(j) sum(p[,j]*I, na.rm=TRUE))
  MeanRperDate.WT <- sapply(1:nb_time_periods, function(i) mean(rep(MeanRperIndexCaseDate[which((1:T >= t_start[i]) * (1:T <= t_end[i]) == 1)], I[which((1:T >= t_start[i]) * (1:T <= t_end[i]) == 1)]) ) )
  
  possibleAncesTime <- sapply(1:T,function(t) (t-(which(si_distr!=0))+1)[which(t-(which(si_distr!=0))+1>0)])
  ancestriesTime <- t(sapply(1:nSim , function(i) DrawOneSetOfAncestries()))
  
  Rsim <- sapply(1:nb_time_periods,function(i) rowSums((ancestriesTime[,]>=t_start[i]) * (ancestriesTime[,]<=t_end[i]),na.rm=TRUE)/sum(I[t_start[i]:t_end[i]]))
  
  R025.WT <- apply(Rsim, 2, quantile,0.025,na.rm=TRUE)
  R025.WT <- R025.WT[which(!is.na(R025.WT))]
  R975.WT <- apply(Rsim, 2, quantile,0.975,na.rm=TRUE)
  R975.WT <- R975.WT[which(!is.na(R975.WT))]
  std.WT <- apply(Rsim, 2, sd,na.rm=TRUE)
  std.WT <- std.WT[which(!is.na(std.WT))]
  
  results<-list()
  
  results$R<-as.data.frame(cbind(t_start,t_end,MeanRperDate.WT,std.WT,R025.WT,R975.WT))
  names(results$R)<-c("t_start","t_end","Mean(R)","Std(R)","Quantile.0.025(R)","Quantile.0.975(R)")
  
  results$method <- method
  results$si_distr <- si_distr
  
  results$SI.Moments<-as.data.frame(cbind(final_mean_si,Finalstd_si))
  names(results$SI.Moments)<-c("Mean","Std")
  
  if(!is.null(dates)) 
    results$dates <- dates
  results$I <- I
  results$I_local <- I
  results$I_local[1] <- 0
  results$I_imported <- c(I[1], rep(0, length(I)-1))
  
  if(plot)
  {
    plots(results, what="all", legend = legend)
  }
  
  return(results)
  
}
