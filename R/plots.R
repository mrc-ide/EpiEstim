#' Plotting the outputs of functions estimating the reproduction number from incidence time series and assumptions regarding the serial interval distribution
#' 
#' \code{plots} allows plotting the outputs of functions \code{\link{EstimateR}} and \code{\link{WT}}
#' 
#' @param x The output of function \code{\link{EstimateR}} or function \code{\link{WT}}, or a list of such outputs. If a list, and \code{what='R'} or \code{what='all'}, all estimates of R are plotted on a single graph. 
#' @param what A string specifying what to plot, namely the incidence time series (\code{what='I'}), the estimated reproduction number (\code{what='R'}), the serial interval distribution (\code{what='SI'}, or all three (\code{what='all'})). 
#' @param add_imported_cases A boolean to specify whether, on the incidence time series plot, to add the incidence of imported cases. 
#' @param ylim For what = "I" or "R"; a parameter similar to that in \code{par}, to monitor the limits of the vertical axis
#' @param options_I For what = "I" or "all". A list of graphical options: 
#'  \describe{
#' \item{col}{A colour or vector of colours used for plotting I. By default uses the default R colours.}
#' \item{transp}{A numeric value between 0 and 1 used to monitor transparency of the bars plotted. Defaults to 0.7.}
#' } 
#' @param options_R For what = "R" or "all". A list of graphical options: 
#'  \describe{
#' \item{col}{A colour or vector of colours used for plotting R. By default uses the default R colours.}
#' \item{transp}{A numeric value between 0 and 1 used to monitor transparency of the 95\%CrI. Defaults to 0.2.}
#' } 
#' @param options_SI For what = "SI" or "all". A list of graphical options: 
#'  \describe{
#' \item{prob_min}{A numeric value between 0 and 1. The SI distributions explored are only shown from time 0 up to the time t so that each distribution explored has probability < \code{prob_min} to be on any time step after t. Defaults to 0.001.}
#' \item{transp}{A numeric value between 0 and 1 used to monitor transparency of the lines. Defaults to 0.25}
#' } 
#' @param legend A boolean (TRUE by default) governing the presence / absence of legends on the plots
#' @return a plot (if \code{what = "I"}, \code{"R"}, or \code{"SI"}) or a \code{\link{grob}} object (if \code{what = "all"}).
# #' @details
#' @seealso \code{\link{EstimateR}} and \code{\link{WT}}
#' @author Rolina van Gaalen \email{rolina.van.gaalen@rivm.nl} and Anne Cori \email{a.cori@imperial.ac.uk} 
# #' @references 
#' @importFrom ggplot2 aes aes_string theme
#' @importFrom scales alpha
#' @importFrom grDevices palette
#' @export
#' @examples 
#' ## load data on pandemic flu in a school in 2009
#' data("Flu2009")
#' 
#' ## estimate the instantaneous reproduction number (method "NonParametricSI")
#' R_i <- EstimateR(Flu2009$Incidence, T.Start=2:26, T.End=8:32, method="NonParametricSI", 
#'           SI.Distr=Flu2009$SI.Distr, plot=FALSE)
#'
#' ## visualise results
#' plots(R_i)
#'
#' ## estimate the instantaneous reproduction number (method "NonParametricSI")
#' R_c <- WT(Flu2009$Incidence, T.Start=2:26, T.End=8:32, method="NonParametricSI", 
#'           SI.Distr=Flu2009$SI.Distr, plot=FALSE)
#'
#' ## produce plot of the incidence 
#'        ## (with, on top of total incidence, the incidence of imported cases), 
#'        ## estimated instantaneous and case reproduction numbers 
#'        ## and serial interval distribution used
#' p_I <- plots(R_i, "I", add_imported_cases=TRUE) # plots the incidence 
#' p_SI <- plots(R_i, "SI") # plots the serial interval distribution
#' p_Ri <- plots(R_i, "R", ylim=c(0,4)) # plots the estimated instantaneous reproduction number
#' p_Rc <- plots(R_c, "R", ylim=c(0,4)) # plots the estimated case reproduction number
#' gridExtra::grid.arrange(p_I,p_SI,p_Ri,p_Rc,ncol=2)
#' 
#' @import reshape2 grid gridExtra
#' @importFrom ggplot2 last_plot ggplot aes aes_string geom_step ggtitle geom_ribbon geom_line xlab ylab xlim geom_hline ylim geom_histogram scale_colour_manual scale_fill_manual scale_linetype_manual lims
#' @importFrom plotly layout mutate arrange rename summarise filter ggplotly
#' @importFrom graphics plot
#' @importFrom incidence as.incidence
plots <- function(x = NULL, what=c("all", "I", "R", "SI"), add_imported_cases=FALSE, ylim=NULL, 
                  options_I = list(col = palette(), transp = 0.7),
                  options_R = list(col = palette(), transp = 0.2),
                  options_SI = list(prob_min = 0.001, transp = 0.25), 
                  legend = TRUE) {
  
  if (is.null(x)) {
    stop("plots requires non NULL x input.")
  }
  
  # check if x is a single output of EpiEstim or a list of such outputs
  if(is.data.frame(x[[1]])) # x is a single output of EpiEstim
  {
    multiple_input <- FALSE
    options_R$col <- options_R$col[1]
  }else
  {
    multiple_input <- TRUE
    x_list <- x
    x <- x_list[[1]]
    if(length(x_list)>length(col))
    {
      warnings("color vector too short, recycling colors.")
      options_R$col <- rep(options_R$col, ceiling(length(x_list) / length(options_R$col)))
      options_R$col <- options_R$col[1:length(x_list)]
    }else
    {
      options_R$col <- options_R$col[1:length(x_list)]
    }
  }
  
  T.Start <- x$R$T.Start 
  T.End <- x$R$T.End
  Mean.Posterior <- x$R[, "Mean(R)"]
  Quantile.0.025.Posterior <- x$R[, "Quantile.0.025(R)"]
  Quantile.0.975.Posterior <- x$R[, "Quantile.0.975(R)"]
  method <- x$method
  SI.Distr <- x$SI.Distr
  I <- data.frame(local=x$I_local, imported=x$I_imported)
  T<-nrow(I)
  
  ########################################################################
  ### these few lines are to make CRAN checks happy with ggplot2... ###
  Time <- NULL
  Incidence <- NULL
  Incidence_imported <- NULL
  value <- NULL
  meanR <- NULL
  meanR2 <- NULL
  group <- NULL
  lower2 <- NULL
  lower <- NULL
  upper <- NULL
  upper2 <- NULL
  Times <- NULL
  ..density.. <- NULL
  start <- NULL
  end <- NULL
  SI.Distr.1 <- NULL
  ########################################################################
  
  if (method == "UncertainSI" | method == "SIFromData" | method == "SIFromSample") {
    Mean.SI.sample <- x$SI.Moments["Mean"]
    Std.SI.sample <- x$SI.Moments["Std"]
  }
  what <- match.arg(what)
  if (what == "I" | what =="all") {
    if(add_imported_cases)
    {
      p1 <- plot(as.incidence(I), ylab="Incidence", xlab = "Time", color = options_I$col, alpha = options_I$transp) +
        ggtitle("Epidemic curve")
    }else
    {
      p1 <- plot(as.incidence(rowSums(I)), ylab="Incidence", xlab = "Time", color = options_I$col, alpha = options_I$transp) +
        ggtitle("Epidemic curve")
    }
    
    if(!is.null(ylim))
      p1 <- p1 + lims(y=ylim)
  }
  if (what == "R" | what =="all") {
    
    time.points <- apply(x$R[,c("T.Start","T.End") ], 1, function(x) x[1]:(x[2]-1)) 
    if (length(time.points) == length(unique(matrix(time.points,ncol=1)))) { 
      
      if(!multiple_input)
      {
        if(is.null(ylim))
          ylim <- c(0,max(Quantile.0.975.Posterior, na.rm = TRUE))
        
        df <- melt(data.frame(start=T.Start, end=T.End, meanR=Mean.Posterior, lower=Quantile.0.025.Posterior,
                              upper=Quantile.0.975.Posterior), id=c("meanR", "lower", "upper")) 
        df$group <- as.factor(rep(1:length(T.Start), dim(df)[1]/length(T.Start)))
        
        p2 <- ggplot(df, aes(x=as.numeric(value), y=as.numeric(meanR), group=as.factor(group))) +
          geom_ribbon(aes(ymin=lower, ymax=upper, fill="95%CrI")) +
          geom_line(aes(y = meanR, colour="Mean")) +
          xlab("Time") +
          ylab("R") +
          xlim(c(1,max(T.End))) +
          ylim(ylim) +
          scale_colour_manual("",values=options_R$col) +
          scale_fill_manual("",values=alpha(options_R$col, options_R$transp)) +
          ggtitle("Estimated R")
        
      }else
      {
        df_tmp <- data.frame(start=T.Start, end=T.End, meanR=Mean.Posterior, lower=Quantile.0.025.Posterior,
                             upper=Quantile.0.975.Posterior)
        df <- df_tmp
        id_tmp <- c("meanR", "lower", "upper")
        id <- id_tmp
        
        for(i in 2:length(x_list))
        {
          x2 <- x_list[[i]]
          T.Start2 <- x2$R$T.Start 
          T.End2 <- x2$R$T.End
          Mean.Posterior2 <- x2$R[, "Mean(R)"]
          Quantile.0.025.Posterior2 <- x2$R[, "Quantile.0.025(R)"]
          Quantile.0.975.Posterior2 <- x2$R[, "Quantile.0.975(R)"]  
          df_tmp2 <- data.frame(start2=T.Start2, end2=T.End2, meanR2=Mean.Posterior2, lower2=Quantile.0.025.Posterior2,
                                upper2=Quantile.0.975.Posterior2)
          names(df_tmp2) <- paste0(names(df_tmp), i)
          df <- cbind(df, df_tmp2)
          id_tmp2 <- paste0(id, i)
          id <- c(id, id_tmp2)
        }
        
        if(is.null(ylim))
          ylim <- c(0,max(df[,grep("upper", names(df))], na.rm = TRUE))
        
        df <- melt(df, id=id) 
        df$group <- as.factor(rep(1:length(T.Start), dim(df)[1]/length(T.Start)))
        
        p2 <- ggplot(df, aes(x=as.numeric(value), y=as.numeric(meanR), group=as.factor(group))) +
          geom_ribbon(aes(ymin=lower, ymax=upper, fill="95%CrI")) +
          geom_line(aes(y = meanR, colour="Mean")) 
        
        for(i in 2:length(x_list))
        {
          p2 <- p2 + geom_ribbon(aes_string(ymin=paste0("lower",i), ymax=paste0("upper",i), fill=shQuote(paste0("95%CrI",i)))) +
            geom_line(aes_string(y = paste0("meanR",i), colour=shQuote(paste0("Mean",i))))
        }
        
        p2 <- p2 +
          xlab("Time") +
          ylab("R") +
          xlim(c(1,max(T.End))) +
          ylim(ylim) +
          scale_colour_manual("",values=options_R$col) +
          scale_fill_manual("",values=alpha(options_R$col, options_R$transp)) +
          ggtitle("Estimated R")
        
      }
      
    } else { 
      
      if(!multiple_input)
      {
        if(is.null(ylim))
          ylim <- c(0,max(Quantile.0.975.Posterior, na.rm = TRUE))
        
        
        p2 <- ggplot(data.frame(start=T.Start, end=T.End, meanR=Mean.Posterior, lower=Quantile.0.025.Posterior,
                                upper=Quantile.0.975.Posterior), aes(end, meanR)) +
          geom_ribbon(aes(ymin=lower, ymax=upper, fill="95%CrI")) +
          geom_line(aes(colour="Mean")) +
          geom_hline(yintercept=1, linetype="dotted") +
          xlab("Time") +
          ylab("R") +
          xlim(c(1,max(T.End))) +
          ylim(ylim) +
          ggtitle("Estimated R") +
          scale_colour_manual("",values=options_R$col)+
          scale_fill_manual("",values=alpha(options_R$col, options_R$transp))
      }else
      {
        ####
        
        df_tmp <- data.frame(start=T.Start, end=T.End, meanR=Mean.Posterior, lower=Quantile.0.025.Posterior,
                             upper=Quantile.0.975.Posterior)
        df <- df_tmp
        
        for(i in 2:length(x_list))
        {
          x2 <- x_list[[i]]
          T.Start2 <- x2$R$T.Start 
          T.End2 <- x2$R$T.End
          Mean.Posterior2 <- x2$R[, "Mean(R)"]
          Quantile.0.025.Posterior2 <- x2$R[, "Quantile.0.025(R)"]
          Quantile.0.975.Posterior2 <- x2$R[, "Quantile.0.975(R)"]  
          df_tmp2 <- data.frame(start2=T.Start2, end2=T.End2, meanR2=Mean.Posterior2, lower2=Quantile.0.025.Posterior2,
                                upper2=Quantile.0.975.Posterior2)
          names(df_tmp2) <- paste0(names(df_tmp), i)
          df <- cbind(df, df_tmp2)
          
        }
        
        if(is.null(ylim))
          ylim <- c(0,max(df[,grep("upper", names(df))], na.rm = TRUE))
        
        
        p2 <- ggplot(df, aes(end, meanR)) +
          geom_ribbon(aes(ymin=lower, ymax=upper, fill="95%CrI")) +
          geom_line(aes(y = meanR, colour="Mean")) 
        
        for(i in 2:length(x_list))
        {
          p2 <- p2 + 
            geom_ribbon(aes_string(ymin=paste0("lower",i), ymax=paste0("upper",i), fill=shQuote(paste0("95%CrI",i)))) +
            geom_line(aes_string(y = paste0("meanR",i), colour=shQuote(paste0("Mean",i))))
        }
        
        p2 <- p2 +
          geom_hline(yintercept=1, linetype="dotted") +
          xlab("Time") +
          ylab("R") +
          xlim(c(1,max(T.End))) +
          ylim(ylim) +
          ggtitle("Estimated R") +
          scale_colour_manual("",values=options_R$col)+
          scale_fill_manual("",values=alpha(options_R$col, options_R$transp))
      }
      
    }
    
  }
  if (what == "SI" | what == "all") {
    
    if (method == "UncertainSI" | method == "SIFromData" | method == "SIFromSample") {
      
      tmp <- cumsum(apply(SI.Distr,2,max) >= options_SI$prob_min)
      stop_at <- min(which(tmp ==tmp[length(tmp)]))
      
      SI.Distr_for_plot <- SI.Distr[,1:stop_at]
      
      dataL = melt(t(SI.Distr_for_plot))
      p3  <- ggplot(dataL, aes_string(x="Var1", y="value", group="Var2")) + 
        geom_line(alpha=options_SI$transp) +
        ggtitle("Explored SI distributions") + 
        xlab("Time") +
        ylab("Frequency") 
      
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
      
    }
  }
  
  if(what == "I")
  {
    if(!legend) p1 <- p1 + theme(legend.position="none")
    return(p1)
  }
  if(what == "R")
  {
    if(!legend) p2 <- p2 + theme(legend.position="none")
    return(p2)
  }
  if(what == "SI")
  {
    if(!legend) p3 <- p3 + theme(legend.position="none")
    return(p3)
  }
  if(what == "all")
  {
    if(!legend) 
    {
      p1 <- p1 + theme(legend.position="none")
      p2 <- p2 + theme(legend.position="none")
      p3 <- p3 + theme(legend.position="none")
    }
    
    out <- list(I=p1, SI=p3, R=p2)
    out.grid <- arrangeGrob(grobs=out, nrow = 3, ncol=1)
    grid.arrange(out.grid, newpage = FALSE)
    return(out.grid)
  }
  
}
