#########################################################################################################################
# WT function to estimate Rc the case reproduction number #
#########################################################################################################################

#' Estimation of the case reproduction number using the Wallinga and Teunis method
#' 
#' \code{WT} estimates the case reproduction number of an epidemic, given the incidence time series and the serial interval distribution. 
#' 
#' @param I Vector of non-negative integers containing an incidence time series.
#' @param T.Start Vector of positive integers giving the starting times of each window over which the reproduction number will be estimated. These must be in ascending order, and so that for all \code{i}, \code{T.Start[i]<=T.End[i]}. T.Start[1] should be strictly after the first day with non null incidence.
#' @param T.End Vector of positive integers giving the ending times of each window over which the reproduction number will be estimated. These must be in ascending order, and so that for all \code{i}, \code{T.Start[i]<=T.End[i]}.
#' @param method One of "NonParametricSI" or "ParametricSI" (see details).
#' @param Mean.SI For method "ParametricSI" ; positive real giving the mean serial interval.
#' @param Std.SI For method "ParametricSI" ; non negative real giving the stadard deviation of the serial interval.
#' @param SI.Distr For method "NonParametricSI" ; vector of probabilities giving the discrete distribution of the serial interval, starting with \code{SI.Distr[1]} (probability that the serial interval is zero), which should be zero.
#' @param nSim A positive integer giving the number of simulated epidemic trees used for computation of the confidence intervals of the case reproduction number (see details).
#' @param plot Logical. If \code{TRUE} (default is \code{FALSE}), output is plotted (see value).
#' @param leg.pos One of "\code{bottomright}", "\code{bottom}", "\code{bottomleft}", "\code{left}", "\code{topleft}", "\code{top}", "\code{topright}", "\code{right}", "\code{center}" or \code{\link{xy.coords}(x, y)}, with \code{x} and \code{y} real numbers. 
#' This specifies the position of the legend in the plot. Alternatively, \code{locator(1)} can be used ; the user will then need to click where the legend needs to be written.
#' @return {
#' 	a list with components: 
#' 	\itemize{
#' 	\item{R}{: a dataframe containing: 
#' 	    the times of start and end of each time window considered ; 
#' 	  the estimated mean, std, and 0.025 and 0.975 quantiles of the reproduction number for each time window.}
#' 	\item{SIDistr}{: a dataframe containing: 
#' 	    for method "NonParametricSI", the mean and standard deviation of the discrete serial interval distribution;
#' 	for method "ParametricSI", the discrete serial interval distribution.}
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
#'    SI.Distr=Flu2009$SI.Distr, plot=TRUE, leg.pos=xy.coords(1,1.75), nSim=100)
#' # the second plot produced shows, at each each day, 
#' # the estimate of the cqse reproduction number over the 7-day window finishing on that day.
#' 
#' ## estimate the case reproduction number (method "ParametricSI")
#' WT(Flu2009$Incidence, T.Start=2:26, T.End=8:32, method="ParametricSI", 
#'    Mean.SI=2.6, Std.SI=1.5, plot=TRUE, nSim=100)
#' # the second plot produced shows, at each each day, 
#' # the estimate of the cqse reproduction number over the 7-day window finishing on that day.
WT <- function(I,T.Start,T.End,method=c("NonParametricSI","ParametricSI"),Mean.SI=NULL,Std.SI=NULL,SI.Distr=NULL,nSim=10,plot=FALSE,leg.pos="topright")
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
  
  ### Error messages ###
  
  method <- match.arg(method)
  
  if(!is.vector(I))
  {
    stop("I must be a vector.")
  }
  T<-length(I)
  for(i in 1:T)
  {
    if(I[i]<0)
    {
      stop("I must be a positive vector.")
    }
  }
  I[which(is.na(I))] <- 0
  
  if(!is.vector(T.Start))
  {
    stop("T.Start must be a vector.")
  }
  if(!is.vector(T.End))
  {
    stop("T.End must be a vector.")
  }
  if(length(T.Start)!=length(T.End))
  {
    stop("T.Start and T.End must have the same length.")
  }
  NbTimePeriods<-length(T.Start)
  for(i in 1:NbTimePeriods)
  {
    if(T.Start[i]>T.End[i])
    {
      stop("T.Start[i] must be <= T.End[i] for all i.")
    }
    if(T.Start[i]<1 || T.Start[i]%%1!=0)
    {
      stop("T.Start must be a vector of >0 integers.")
    }
    if(T.End[i]<1 || T.End[i]%%1!=0)
    {
      stop("T.End must be a vector of >0 integers.")
    }
  }
  
  if(method=="NonParametricSI")
  {
    if(is.null(SI.Distr))
    {
      stop("method NonParametricSI requires to specify the SI.Distr argument.")
    }
    if(!is.vector(SI.Distr))
    {
      stop("method NonParametricSI requires that SI.Distr must be a vector.")
    }
    if(SI.Distr[1]!=0)
    {
      stop("method NonParametricSI requires that SI.Distr[1] = 0.")
    }
    if(length(SI.Distr)>1)
    {
      for(i in 2:length(SI.Distr))
      {
        if(SI.Distr[i]<0)
        {
          stop("method NonParametricSI requires that SI.Distr must be a positive vector.")
        }
      }
    }
    if(abs(sum(SI.Distr)-1)>0.01)
    {
      stop("method NonParametricSI requires that SI.Distr must sum to 1.")
    }
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
  
  if(ParametricSI=="Y") # method "ParametricSI"
  {
    if(length(which(abs(cumsum(SI.Distr)-1)<0.01))==0)
    {
      warning("The serial interval distribution you have chosen is very wide compared to the duration of the epidemic.\nEstimation will be performed anyway but restults should be interpreted with care.")
      MaxT <- length(cumsum(SI.Distr))
    }else
    {
      MaxT<-min(which(abs(cumsum(SI.Distr)-1)<0.01))
    }
    results$SIDistr<-as.data.frame(cbind(0:(MaxT-1),SI.Distr[1:MaxT]))
    names(results$SIDistr)<-c("k","w[k]")
  }else # method "NonParametricSI"
  {
    results$SIDistr<-as.data.frame(cbind(FinalMean.SI,FinalStd.SI))
    names(results$SIDistr)<-c("Mean Discrete SI","Std Discrete SI")
  }
  
  if(plot)
  {
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
    
    p1 <- ggplot(data.frame(Time=1:T, Incidence=I), aes(x=Time, y=Incidence)) +
      geom_step() +
      ggtitle("Epidemic curve")
    p1ly <- ggplotly(p1)
    
    # test if intervals overlap 
    time.points <- apply(results$R[,c("T.Start","T.End") ], 1, function(x) x[1]:(x[2]-1)) 
    if (length(time.points) == length(unique(matrix(time.points,ncol=1)))) { 
      
      df <- melt(data.frame(start=T.Start, end=T.End, meanR=MeanRperDate.WT, lower=R025.WT,
                            upper=R975.WT), id=c("meanR", "lower", "upper")) 
      df$group <- as.factor(rep(1:length(T.Start), dim(df)[1]/length(T.Start)))
      
      p2 <- ggplot(df, aes(x=as.numeric(value), y=as.numeric(meanR), group=as.factor(group))) +
        geom_ribbon(aes(ymin=lower, ymax=upper), colour=NA, fill="black", alpha=0.2) +
        geom_line() +
        xlab("Time") +
        ylab("R") +
        xlim(c(1,max(T.End))) +
        ggtitle("Estimated Rc")
      p2ly <- ggplotly(p2)
      
    } else { 
      
      p2 <- ggplot(data.frame(start=T.Start, end=T.End, meanR=MeanRperDate.WT, lower=R025.WT,
                              upper=R975.WT), aes(end, meanR)) +
        geom_ribbon(aes(ymin=lower, ymax=upper), fill="grey") +
        geom_line() +
        geom_hline(yintercept=1, linetype="dotted") +
        xlab("Time") +
        ylab("R") +
        xlim(c(1,max(T.End))) +
        ylim(c(0,max(R975.WT, na.rm = TRUE))) +
        ggtitle("Estimated Rc") 
      p2ly <- ggplotly(p2)
      #+
      #legend(leg.pos, c("Median", "95%CrI"), col = c("Black", 
      #               grey), lwd = c(1, 10), bty = "n", cex = 1.2)
    }
    
    SI.Distr.times <- unlist(apply(data.frame(0:(length(SI.Distr) - 1), SI.Distr), 1,
                                   function(x) {if (x[2]!=0) unlist(rep(x[1],round(x[2]*1000)),use.names=FALSE)}))
    names(SI.Distr.times) <- NULL
    
    p3 <- ggplot(data.frame(Times=SI.Distr.times), aes(0.5+Times)) +
      geom_histogram(binwidth=1, aes(y=..density..)) +
      xlab("Time") + 
      xlim(c(0,0.5+max(SI.Distr.times))) + 
      ylab("Frequency") + 
      ggtitle("Serial interval distribution") 
    p3ly <- ggplotly(p3)
    
    grid.arrange(p1,p3,p2,ncol=1)
  }
  
  results$method <- method
  results$SI.Distr <- SI.Distr
  results$I <- I
  return(results)
  
}
