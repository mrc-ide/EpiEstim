#' Plotting the outputs of functions estimating the reproduction number from incidence time series and assumptions regarding the serial interval distribution
#' 
#' \code{plots} allows plotting the outputs of functions \code{\link{EstimateR}} and \code{\link{WT}}
#' 
#' @param x The output of function \code{\link{EstimateR}} or function \code{\link{WT}}
#' @param what A string specifying what to plot, namely the incidence time series (\code{what='I'}), the estimated reproduction number (\code{what='R'}), or the serial interval distribution (\code{what='SI'}). 
#' @return a plot or a list of plots (if \code{what == "SI"} and \code{x$method == "UncertainSI"})
# #' @details
#' @seealso \code{\link{EstimateR}} and \code{\link{WT}}
#' @author Rolina van Gaalen \email{rolina.van.gaalen@rivm.nl} and Anne Cori \email{a.cori@imperial.ac.uk} 
# #' @references 
#' @export
#' @examples 
#' ## load data on pandemic flu in a school in 2009
#' data("Flu2009")
#' 
#' ## estimate the instantaneous reproduction number (method "NonParametricSI")
#' R_i <- EstimateR(Flu2009$Incidence, T.Start=2:26, T.End=8:32, method="NonParametricSI", 
#'           SI.Distr=Flu2009$SI.Distr, plot=FALSE)
#'
#' ## estimate the instantaneous reproduction number (method "NonParametricSI")
#' R_c <- WT(Flu2009$Incidence, T.Start=2:26, T.End=8:32, method="NonParametricSI", 
#'           SI.Distr=Flu2009$SI.Distr, plot=FALSE)
#'
#' ## produce plot of the incidence, 
#'        ## estimated instantaneous and case reproduction numbers 
#'        ## and serial interval distribution used
#' p_I <- plots(R_i, "I") # plots the incidence 
#' p_SI <- plots(R_i, "SI") # plots the serial interval distribution
#' p_Ri <- plots(R_i, "R") # plots the estimated instantaneous reproduction number
#' p_Rc <- plots(R_c, "R") # plots the estimated case reproduction number
#' gridExtra::grid.arrange(p_I,p_SI,p_Ri,p_Rc,ncol=2)
#' 
#' @import reshape2 grid gridExtra
#' @importFrom ggplot2 last_plot ggplot aes geom_step ggtitle geom_ribbon geom_line xlab ylab xlim geom_hline ylim geom_histogram
#' @importFrom plotly layout mutate arrange rename summarise filter ggplotly
plots <- function(x=NULL, what=c("I", "R", "SI")) {
  
  if (is.null(x)) {
    stop("plots requires non NULL x input.")
  }
  
  T.Start <- x$R$T.Start 
  T.End <- x$R$T.End
  Mean.Posterior <- x$R[, "Mean(R)"]
  Quantile.0.025.Posterior <- x$R[, "Quantile.0.025(R)"]
  Quantile.0.975.Posterior <- x$R[, "Quantile.0.975(R)"]
  method <- x$method
  SI.Distr <- x$SI.Distr
  I <- x$I
  
  if(is.vector(I))
  {
    I_tmp <- I
    I <- data.frame(local=I_tmp, imported=rep(0, length(I_tmp)))
  }else
  {
    if(!is.data.frame(I) | !all(c("local","imported") %in% names(I)) ) 
    {
      stop("I must be a vector or a dataframe with 2 columns called 'local' and 'imported'.")
    }
  }
  T<-nrow(I)
  
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
  
  if (method == "UncertainSI") {
    Mean.SI.sample <- x$SIDistr["Mean.SI.sample"]
    Std.SI.sample <- x$SIDistr["Std.SI.sample"]
  }
  what <- match.arg(what)
  if (what == "I") {
    p1 <- ggplot(data.frame(Time=1:T, Incidence=rowSums(I)), aes(x=Time, y=Incidence)) +
      geom_step() +
      ggtitle("Epidemic curve")
    p1ly <- ggplotly(p1)
    print(p1ly)
    return(p1)
  }else if (what == "R") {
    
    time.points <- apply(x$R[,c("T.Start","T.End") ], 1, function(x) x[1]:(x[2]-1)) 
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
      p2ly <- ggplotly(p2)
      
      print(p2ly)
      return(p2)
    } else { 
      
      p2 <- ggplot(data.frame(start=T.Start, end=T.End, meanR=Mean.Posterior, lower=Quantile.0.025.Posterior,
                              upper=Quantile.0.975.Posterior), aes(end, meanR)) +
        geom_ribbon(aes(ymin=lower, ymax=upper), fill="grey") +
        geom_line() +
        geom_hline(yintercept=1, linetype="dotted") +
        xlab("Time") +
        ylab("R") +
        xlim(c(1,max(T.End))) +
        ylim(c(0,max(Quantile.0.975.Posterior, na.rm = TRUE))) +
        ggtitle("Estimated R") 
      p2ly <- ggplotly(p2)
      
      print(p2ly)
      return(p2)
    }
    
  } else if (what == "SI") {
    
    if (method == "UncertainSI") {
      
      p3 <- ggplot(data.frame(Mean.SI.sample), aes(Mean.SI.sample)) +
        geom_histogram(bins=30) +
        xlab("Mean serial interval") +
        ylab("Density") + 
        ggtitle("Explored \n mean serial intervals")
      
      p4 <- ggplot(data.frame(Std.SI.sample), aes(Std.SI.sample)) +
        geom_histogram(bins=30) +
        xlab("Std serial interval") +
        ggtitle("Explored \n std serial intervals")
      
      grid.arrange(p3,p4,ncol=2)
      return(list(p3, p4))
      
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
      p3ly <- ggplotly(p3)
      
      print(p3ly)
      return(p3)
      
    }
  }
  
}
