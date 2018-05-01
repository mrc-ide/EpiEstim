#' Plotting the outputs of functions estimating the reproduction number from incidence time series and assumptions regarding the serial interval distribution
#' 
#' \code{plots} allows plotting the outputs of functions \code{\link{estimate_r}} and \code{\link{WT}}
#' 
#' @param x The output of function \code{\link{estimate_r}} or function \code{\link{WT}}, or a list of such outputs. If a list, and \code{what='R'} or \code{what='all'}, all estimates of R are plotted on a single graph. 
#' @param what A string specifying what to plot, namely the incidence time series (\code{what='I'}), the estimated reproduction number (\code{what='R'}), the serial interval distribution (\code{what='SI'}, or all three (\code{what='all'})). 
#' @param add_imported_cases A boolean to specify whether, on the incidence time series plot, to add the incidence of imported cases. 
#' @param options_I For what = "I" or "all". A list of graphical options: 
#'  \describe{
#' \item{col}{A colour or vector of colours used for plotting I. By default uses the default R colours.}
#' \item{transp}{A numeric value between 0 and 1 used to monitor transparency of the bars plotted. Defaults to 0.7.}
#' \item{xlim}{A parameter similar to that in \code{par}, to monitor the limits of the horizontal axis}
#' \item{ylim}{A parameter similar to that in \code{par}, to monitor the limits of the vertical axis}
#' } 
#' @param options_R For what = "R" or "all". A list of graphical options: 
#'  \describe{
#' \item{col}{A colour or vector of colours used for plotting R. By default uses the default R colours.}
#' \item{transp}{A numeric value between 0 and 1 used to monitor transparency of the 95\%CrI. Defaults to 0.2.}
#' \item{xlim}{A parameter similar to that in \code{par}, to monitor the limits of the horizontal axis}
#' \item{ylim}{A parameter similar to that in \code{par}, to monitor the limits of the vertical axis}
#' } 
#' @param options_SI For what = "SI" or "all". A list of graphical options: 
#'  \describe{
#' \item{prob_min}{A numeric value between 0 and 1. The SI distributions explored are only shown from time 0 up to the time t so that each distribution explored has probability < \code{prob_min} to be on any time step after t. Defaults to 0.001.}
#' \item{col}{A colour or vector of colours used for plotting the SI. Defaults to black.}
#' \item{transp}{A numeric value between 0 and 1 used to monitor transparency of the lines. Defaults to 0.25}
#' \item{xlim}{A parameter similar to that in \code{par}, to monitor the limits of the horizontal axis}
#' \item{ylim}{A parameter similar to that in \code{par}, to monitor the limits of the vertical axis}
#' } 
#' @param legend A boolean (TRUE by default) governing the presence / absence of legends on the plots
#' @return a plot (if \code{what = "I"}, \code{"R"}, or \code{"SI"}) or a \code{\link{grob}} object (if \code{what = "all"}).
# #' @details
#' @seealso \code{\link{estimate_r}} and \code{\link{WT}}
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
#' ## estimate the instantaneous reproduction number (method "non_parametric_si")
#' R_i <- estimate_r(Flu2009$incidence, method="non_parametric_si",
#'                  config=list(t_start=2:26, t_end=8:32, 
#'                              si_distr=Flu2009$si_distr, plot=FALSE)
#'                 )
#'
#' ## visualise results
#' plots(R_i, legend = FALSE)
#'
#' ## estimate the instantaneous reproduction number (method "non_parametric_si")
#' R_c <- WT(Flu2009$incidence, t_start=2:26, t_end=8:32, method="non_parametric_si", 
#'           si_distr=Flu2009$si_distr, plot=FALSE)
#'
#' ## produce plot of the incidence 
#'        ## (with, on top of total incidence, the incidence of imported cases), 
#'        ## estimated instantaneous and case reproduction numbers 
#'        ## and serial interval distribution used
#' p_I <- plots(R_i, "I", add_imported_cases=TRUE) # plots the incidence 
#' p_SI <- plots(R_i, "SI") # plots the serial interval distribution
#' p_Ri <- plots(R_i, "R", 
#'           options_R = list(ylim=c(0,4))) 
#'           # plots the estimated instantaneous reproduction number
#' p_Rc <- plots(R_c, "R", 
#'           list(ylim=c(0,4))) 
#'           # plots the estimated case reproduction number
#' gridExtra::grid.arrange(p_I,p_SI,p_Ri,p_Rc,ncol=2)
#' 
#' @import reshape2 grid gridExtra
#' @importFrom ggplot2 last_plot ggplot aes aes_string geom_step ggtitle geom_ribbon geom_line xlab ylab xlim geom_hline ylim geom_histogram scale_colour_manual scale_fill_manual scale_linetype_manual lims
#' @importFrom plotly layout mutate arrange rename summarise filter ggplotly
#' @importFrom graphics plot
#' @importFrom incidence as.incidence
plots <- function(x = NULL, what=c("all", "I", "R", "SI"), add_imported_cases=FALSE, 
                  options_I = list(col = palette(), transp = 0.7, xlim = NULL, ylim=NULL),
                  options_R = list(col = palette(), transp = 0.2, xlim = NULL, ylim=NULL),
                  options_SI = list(prob_min = 0.001, col = "black", transp = 0.25, xlim = NULL, ylim=NULL), 
                  legend = TRUE) {
  
  if (is.null(x)) {
    stop("plots requires non NULL x input.")
  }
  
  # dealing with the fact that some options may be left to default but others may have been specified by user
  if (is.null(options_I$col)) options_I$col <- palette()
  if (is.null(options_I$transp)) options_I$transp <- 0.7
  
  if (is.null(options_R$col)) options_R$col <- palette()
  if (is.null(options_R$transp)) options_R$transp <- 0.2
  
  if (is.null(options_SI$prob_min)) options_SI$prob_min <- 0.001
  if (is.null(options_SI$col)) options_SI$col <- "black"
  if (is.null(options_SI$transp)) options_SI$transp <- 0.25
  
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
  
  t_start <- x$R$t_start 
  t_end <- x$R$t_end
  mean_posterior <- x$R[, "Mean(R)"]
  quantile_0.025_posterior <- x$R[, "Quantile.0.025(R)"]
  quantile_0.975_posterior <- x$R[, "Quantile.0.975(R)"]
  method <- x$method
  si_distr <- x$si_distr
  I <- data.frame(local=x$I_local, imported=x$I_imported)
  T<-nrow(I)
  if(!is.null(x$dates))
  {
    dates <- x$dates
  }else
  {
    dates <- 1:T
  }
  
  ########################################################################
  ### these few lines are to make CRAN checks happy with ggplot2... ###
  Time <- NULL
  incidence <- NULL
  incidence_imported <- NULL
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
  si_distr.1 <- NULL
  ########################################################################
  
  if (method == "uncertain_si" | method == "si_from_data" | method == "si_from_sample") {
    mean_si.sample <- x$SI.Moments["Mean"]
    std_si.sample <- x$SI.Moments["Std"]
  }
  what <- match.arg(what)
  if (what == "I" | what =="all") {
    if(add_imported_cases)
    {
      p1 <- plot(as.incidence(I, dates = x$dates), ylab="Incidence", xlab = "Time", color = options_I$col, alpha = options_I$transp) +
        ggtitle("Epidemic curve")
    }else
    {
      p1 <- plot(as.incidence(rowSums(I), dates = x$dates), ylab="Incidence", xlab = "Time", color = options_I$col, alpha = options_I$transp) +
        ggtitle("Epidemic curve")
    }
    
    if(!is.null(options_I$xlim))
      p1 <- p1 + lims(x=options_I$xlim)
    
    if(!is.null(options_I$ylim))
      p1 <- p1 + lims(y=options_I$ylim)
  }
  if (what == "R" | what =="all") {
    
    time.points <- apply(x$R[,c("t_start","t_end") ], 1, function(x) x[1]:(x[2]-1)) 
    if (length(time.points) == length(unique(matrix(time.points,ncol=1)))) { 
      
      if(!multiple_input)
      {
        if(is.null(options_R$ylim))
          options_R$ylim <- c(0,max(quantile_0.975_posterior, na.rm = TRUE))
        
        if(is.null(options_R$xlim))
          options_R$xlim <- c(min(dates),max(dates)+1)
        
        df <- melt(data.frame(start=dates[t_start], end=dates[t_end], meanR=mean_posterior, lower=quantile_0.025_posterior,
                              upper=quantile_0.975_posterior), id=c("meanR", "lower", "upper")) 
        df$group <- as.factor(rep(1:length(t_start), dim(df)[1]/length(t_start)))
        
        p2 <- ggplot(df, aes(x=value, y=as.numeric(meanR), group=as.factor(group))) +
          geom_ribbon(aes(ymin=lower, ymax=upper, fill="95%CrI")) +
          geom_line(aes(y = meanR, colour="Mean")) +
          xlab("Time") +
          ylab("R") +
          xlim(options_R$xlim) +
          ylim(options_R$ylim) +
          scale_colour_manual("",values=options_R$col) +
          scale_fill_manual("",values=alpha(options_R$col, options_R$transp)) +
          ggtitle("Estimated R")
        
      }else
      {
        df_tmp <- data.frame(start=dates[t_start], end=dates[t_end], meanR=mean_posterior, lower=quantile_0.025_posterior,
                             upper=quantile_0.975_posterior)
        df <- df_tmp
        id_tmp <- c("meanR", "lower", "upper")
        id <- id_tmp
        
        for(i in 2:length(x_list))
        {
          x2 <- x_list[[i]]
          t_start2 <- x2$R$t_start 
          if(!is.null(x2$dates))
          {
            dates2 <- x2$dates
          }else
          {
            dates2 <- 1:T
          }
          mean_posterior2 <- x2$R[, "Mean(R)"]
          quantile_0.025_posterior2 <- x2$R[, "Quantile.0.025(R)"]
          quantile_0.975_posterior2 <- x2$R[, "Quantile.0.975(R)"]  
          df_tmp2 <- data.frame(start2=dates2[t_start2], end2=dates2[t_end], meanR2=mean_posterior2, lower2=quantile_0.025_posterior2,
                                upper2=quantile_0.975_posterior2)
          names(df_tmp2) <- paste0(names(df_tmp), i)
          df <- cbind(df, df_tmp2)
          id_tmp2 <- paste0(id, i)
          id <- c(id, id_tmp2)
        }
        
        if(is.null(options_R$ylim))
          options_R$ylim <- c(0,max(df[,grep("upper", names(df))], na.rm = TRUE))
        
        if(is.null(options_R$xlim))
          options_R$xlim <- c(min(dates),max(dates)+1)
        
        df <- melt(df, id=id) 
        df$group <- as.factor(rep(1:length(t_start), dim(df)[1]/length(t_start)))
        
        p2 <- ggplot(df, aes(x=value, y=as.numeric(meanR), group=as.factor(group))) +
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
          xlim(options_R$xlim) +
          ylim(options_R$ylim) +
          scale_colour_manual("",values=options_R$col) +
          scale_fill_manual("",values=alpha(options_R$col, options_R$transp)) +
          ggtitle("Estimated R")
        
      }
      
    } else { 
      
      if(!multiple_input)
      {
        if(is.null(options_R$ylim))
          options_R$ylim <- c(0,max(quantile_0.975_posterior, na.rm = TRUE))
        
        if(is.null(options_R$xlim))
          options_R$xlim <- c(min(dates),max(dates)+1)
        
        p2 <- ggplot(data.frame(start=dates[t_start], end=dates[t_end], meanR=mean_posterior, lower=quantile_0.025_posterior,
                                upper=quantile_0.975_posterior), aes(end, meanR)) +
          geom_ribbon(aes(ymin=lower, ymax=upper, fill="95%CrI")) +
          geom_line(aes(colour="Mean")) +
          geom_hline(yintercept=1, linetype="dotted") +
          xlab("Time") +
          ylab("R") +
          xlim(options_R$xlim) +
          ylim(options_R$ylim) +
          ggtitle("Estimated R") +
          scale_colour_manual("",values=options_R$col)+
          scale_fill_manual("",values=alpha(options_R$col, options_R$transp))
      }else
      {
        ####
        
        df_tmp <- data.frame(start=dates[t_start], end=dates[t_end], meanR=mean_posterior, lower=quantile_0.025_posterior,
                             upper=quantile_0.975_posterior)
        df <- df_tmp
        
        for(i in 2:length(x_list))
        {
          x2 <- x_list[[i]]
          t_start2 <- x2$R$t_start 
          if(!is.null(x2$dates))
          {
            dates2 <- x2$dates
          }else
          {
            dates2 <- 1:T
          }
          mean_posterior2 <- x2$R[, "Mean(R)"]
          quantile_0.025_posterior2 <- x2$R[, "Quantile.0.025(R)"]
          quantile_0.975_posterior2 <- x2$R[, "Quantile.0.975(R)"]  
          df_tmp2 <- data.frame(start2=dates2[t_start2], end2=dates2[t_end], meanR2=mean_posterior2, lower2=quantile_0.025_posterior2,
                                upper2=quantile_0.975_posterior2)
          names(df_tmp2) <- paste0(names(df_tmp), i)
          df <- cbind(df, df_tmp2)
          
        }
        
        if(is.null(options_R$ylim))
          options_R$ylim <- c(0,max(df[,grep("upper", names(df))], na.rm = TRUE))
        
        if(is.null(options_R$xlim))
          options_R$xlim <- c(min(dates),max(dates)+1)
        
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
          xlim(options_R$xlim) +
          ylim(options_R$ylim) +
          ggtitle("Estimated R") +
          scale_colour_manual("",values=options_R$col)+
          scale_fill_manual("",values=alpha(options_R$col, options_R$transp))
      }
      
    }
    
  }
  if (what == "SI" | what == "all") {
    
    if (method == "uncertain_si" | method == "si_from_data" | method == "si_from_sample") {
      
      tmp <- cumsum(apply(si_distr,2,max) >= options_SI$prob_min)
      stop_at <- min(which(tmp == tmp[length(tmp)]))
      
      si_distr_for_plot <- si_distr[,1:stop_at] 
      
      dataL <-melt(t(si_distr_for_plot))
      dataL$Var1 <- 0:(ncol(si_distr_for_plot)-1)
      p3  <- ggplot(dataL, aes_string(x="Var1", y="value", group="Var2")) + 
        geom_line(col=options_SI$col, alpha=options_SI$transp) +
        ggtitle("Explored SI distributions") + 
        xlab("Time") +
        ylab("Frequency") 
      
      if(!is.null(options_SI$xlim))
        p3 <- p3 + lims(x=options_SI$xlim)
      
      if(!is.null(options_SI$ylim))
        p3 <- p3 + lims(y=options_SI$ylim)
      
    } else {
      
      tmp <- cumsum(si_distr >= options_SI$prob_min)
      stop_at <- min(which(tmp == tmp[length(tmp)]))
      
      si_distr_for_plot <- si_distr[1:stop_at]
      
      dataL <- data.frame(Times=0:(length(si_distr_for_plot)-1), SIDistr = si_distr_for_plot)
      p3  <- ggplot(dataL, aes_string(x="Times", y="SIDistr")) + 
        geom_line(col=options_SI$col, alpha=options_SI$transp) +
        ggtitle("Explored SI distribution") + 
        xlab("Time") +
        ylab("Frequency") 
      
      if(!is.null(options_SI$xlim))
        p3 <- p3 + lims(x=options_SI$xlim)
      
      if(!is.null(options_SI$ylim))
        p3 <- p3 + lims(y=options_SI$ylim)
      
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
    
    out <- list(I=p1, R=p2, SI=p3)
    out.grid <- arrangeGrob(grobs=out, nrow = 3, ncol=1)
    grid.arrange(out.grid, newpage = FALSE)
    return(out.grid)
  }
  
}
