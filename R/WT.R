#########################################################################################################################
# WT function to estimate Rc the case reproduction number #
#########################################################################################################################

#' Estimation of the case reproduction number using the Wallinga and Teunis method
#' 
#' \code{WT} estimates the case reproduction number of an epidemic, given the incidence time series and the serial interval distribution. 
#' 
#' @param I Vector (or a dataframe with a column named 'I') of non-negative integers containing an incidence time series.
#'  If the dataframe contains a column \code{I$dates}, this is used for plotting. 
#' \code{I$dates} must be in date format, and the dates must be in a row.
#' @param T.Start Vector of positive integers giving the starting times of each window over which the reproduction number will be estimated. These must be in ascending order, and so that for all \code{i}, \code{T.Start[i]<=T.End[i]}. T.Start[1] should be strictly after the first day with non null incidence.
#' @param T.End Vector of positive integers giving the ending times of each window over which the reproduction number will be estimated. These must be in ascending order, and so that for all \code{i}, \code{T.Start[i]<=T.End[i]}.
#' @param method One of "NonParametricSI" or "ParametricSI" (see details).
#' @param Mean.SI For method "ParametricSI" ; positive real giving the mean serial interval.
#' @param Std.SI For method "ParametricSI" ; non negative real giving the stadard deviation of the serial interval.
#' @param SI.Distr For method "NonParametricSI" ; vector of probabilities giving the discrete distribution of the serial interval, starting with \code{SI.Distr[1]} (probability that the serial interval is zero), which should be zero.
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
#' 	\item{method}{: the method used to estimate R, one of "NonParametricSI", "ParametricSI", "UncertainSI", "SIFromData" or "SIFromSample"}
#' 	\item{SI.Distr}{: a vector containing the discrete serial interval distribution used for estimation}
#' 	\item{SI.Moments}{: a vector containing the mean and std of the discrete serial interval distribution(s) used for estimation}
#' 	\item{I}{: the time series of total incidence}
#' 	\item{I_local}{: the time series of incidence of local cases (so that \code{I_local + I_imported = I})}
#' 	\item{I_imported}{: the time series of incidence of imported cases (so that \code{I_local + I_imported = I})}
#' 	\item{dates}{: if dates were specified in I, a vector of dates corresponding to the incidence time series}
#' 	}
#' 	}
#' @details{
#' Estimates of the case reproduction number for an epidemic over predefined time windows can be obtained, 
#' for a given discrete distribution of the serial interval, as proposed by Wallinga and Teunis (AJE, 2004). 
#' Confidence intervals are obtained by simulating a number (nSim) of possible transmission trees. 
#' 
#' The methods vary in the way the serial interval distribution is specified.
#' 
#' ----------------------- \code{method "NonParametricSI"} -----------------------
#' 
#' The discrete distribution of the serial interval is directly specified in the argument \code{SI.Distr}.
#' 
#' If \code{plot} is \code{TRUE}, 3 plots are produced. 
#' The first one shows the epidemic curve. 
#' The second one shows the posterior mean and 95\% credible interval of the reproduction number. The estimate for a time window is plotted at the end of the time window. 
#' The third plot shows the discrete distribution of the serial interval. 
#'
#' ----------------------- \code{method "ParametricSI"} -----------------------
#' 
#' The mean and standard deviation of the continuous distribution of the serial interval are given in the arguments \code{Mean.SI} and \code{Std.SI}.
#' The discrete distribution of the serial interval is derived automatically using \code{\link{DiscrSI}}.
#'
#' If \code{plot} is \code{TRUE}, 3 plots are produced, which are identical to the ones for \code{method "NonParametricSI"} .
#' }
#' @seealso \code{\link{DiscrSI}}, \code{\link{EstimateR}}
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
#' ## estimate the case reproduction number (method "NonParametricSI")
#' WT(Flu2009$Incidence, T.Start=2:26, T.End=8:32, method="NonParametricSI", 
#'    SI.Distr=Flu2009$SI.Distr, plot=TRUE, nSim=100)
#' # the second plot produced shows, at each each day, 
#' # the estimate of the cqse reproduction number over the 7-day window finishing on that day.
#' 
#' ## estimate the case reproduction number (method "ParametricSI")
#' WT(Flu2009$Incidence, T.Start=2:26, T.End=8:32, method="ParametricSI", 
#'    Mean.SI=2.6, Std.SI=1.5, plot=TRUE, nSim=100)
#' # the second plot produced shows, at each each day, 
#' # the estimate of the case reproduction number over the 7-day window finishing on that day.
WT <- function(I, T.Start, T.End,
               method=c("NonParametricSI","ParametricSI"),
               Mean.SI=NULL, Std.SI=NULL, SI.Distr=NULL, nSim=10, 
               plot=FALSE, legend=FALSE)
{
  
  ### Functions ###
  
  #########################################################
  # Draws a possile transmission tree                     #
  #########################################################
  
  DrawOneSetOfAncestries <- function()
  {
    res <- vector()
    for(t in T.Start[1]:T.End[length(T.End)])
    {
      if(length(which(Onset==t))>0)
      {
        if(length(possibleAncesTime[[t]])>0)
        {
          prob <- SI.Distr[t-possibleAncesTime[[t]]+1]*I[possibleAncesTime[[t]]]
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
  
  if(!is.null(I$dates)) 
  {
    dates <- check_dates(I)
  }else
  {
    dates <- NULL
  }
  I <- process_I_vector(I)
  T<-length(I)
  
  ### Adjusting T.Start and T.End so that at least an incident case has been observed before T.Start[1] ###
  
  i <- 1
  while(sum(I[1:(T.Start[i]-1)])==0)
  {
    i <- i+1
  }
  temp <- which(T.Start<i)
  if(length(temp>0))
  {
    T.Start <- T.Start[-temp]
    T.End <- T.End[-temp]
  }
  
  check_times(T.Start, T.End, T)
  NbTimePeriods <- length(T.Start)
  
  if(method=="NonParametricSI")
  {
    check_SI.Distr(SI.Distr)
    SI.Distr <- c(SI.Distr,0)
  }
  
  if(method=="ParametricSI")
  {
    if(is.null(Mean.SI))
    {
      stop("method NonParametricSI requires to specify the Mean.SI argument.")
    }
    if(is.null(Std.SI))
    {
      stop("method NonParametricSI requires to specify the Std.SI argument.")
    }
    if(Mean.SI<1)
    {
      stop("method ParametricSI requires a value >1 for Mean.SI.")
    }
    if(Std.SI<0)
    {
      stop("method ParametricSI requires a >0 value for Std.SI.")
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
  
  if(method=="NonParametricSI")
  {
    ParametricSI<-"N"
  }
  if(method=="ParametricSI")
  {
    ParametricSI<-"Y"
  }
  
  if(ParametricSI=="Y")
  {
    SI.Distr <- sapply(1:T, function(t) DiscrSI(t-1,Mean.SI,Std.SI))
  }
  if(length(SI.Distr)<T+1){SI.Distr[(length(SI.Distr)+1):(T+1)]<-0}
  FinalMean.SI<-sum(SI.Distr*(0:(length(SI.Distr)-1)))
  FinalStd.SI<-sqrt(sum(SI.Distr*(0:(length(SI.Distr)-1))^2)-FinalMean.SI^2)
  
  TimePeriodsWithNoIncidence <- vector()
  for(i in 1:NbTimePeriods)
  {
    if(sum(I[T.Start[i]:T.End[i]])==0)
    {
      TimePeriodsWithNoIncidence <- c(TimePeriodsWithNoIncidence,i)
    }
  }
  if(length(TimePeriodsWithNoIncidence)>0)
  {
    T.Start <- T.Start[-TimePeriodsWithNoIncidence]
    T.End <- T.End[-TimePeriodsWithNoIncidence]
    NbTimePeriods <- length(T.Start)
  }
  
  Onset <- vector()
  for (t in 1:T) {
    Onset <- c(Onset, rep(t, I[t]))
  }
  NbCases <- length(Onset)
  
  delay <- outer (1:T, 1:T, "-")
  SIdelay <- apply(delay, 2, function(x) SI.Distr[pmin(pmax(x + 1, 1), length(SI.Distr))])
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
  MeanRperDate.WT <- sapply(1:NbTimePeriods, function(i) mean(rep(MeanRperIndexCaseDate[which((1:T >= T.Start[i]) * (1:T <= T.End[i]) == 1)], I[which((1:T >= T.Start[i]) * (1:T <= T.End[i]) == 1)]) ) )
  
  possibleAncesTime <- sapply(1:T,function(t) (t-(which(SI.Distr!=0))+1)[which(t-(which(SI.Distr!=0))+1>0)])
  ancestriesTime <- t(sapply(1:nSim , function(i) DrawOneSetOfAncestries()))
  
  Rsim <- sapply(1:NbTimePeriods,function(i) rowSums((ancestriesTime[,]>=T.Start[i]) * (ancestriesTime[,]<=T.End[i]),na.rm=TRUE)/sum(I[T.Start[i]:T.End[i]]))
  
  R025.WT <- apply(Rsim, 2, quantile,0.025,na.rm=TRUE)
  R025.WT <- R025.WT[which(!is.na(R025.WT))]
  R975.WT <- apply(Rsim, 2, quantile,0.975,na.rm=TRUE)
  R975.WT <- R975.WT[which(!is.na(R975.WT))]
  std.WT <- apply(Rsim, 2, sd,na.rm=TRUE)
  std.WT <- std.WT[which(!is.na(std.WT))]
  
  results<-list()
  
  results$R<-as.data.frame(cbind(T.Start,T.End,MeanRperDate.WT,std.WT,R025.WT,R975.WT))
  names(results$R)<-c("T.Start","T.End","Mean(R)","Std(R)","Quantile.0.025(R)","Quantile.0.975(R)")
  
  results$method <- method
  results$SI.Distr <- SI.Distr
  
  results$SI.Moments<-as.data.frame(cbind(FinalMean.SI,FinalStd.SI))
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
