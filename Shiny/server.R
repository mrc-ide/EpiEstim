# TODO: Add comment
#
# Author: annecori
###############################################################################

library(shiny)
library(EpiEstim)
library(devtools)
## install the hackout version of coarseDataTools
install_github("nickreich/coarseDataTools", ref = "hackout3")
library(coarseDataTools)

## Set your directory to where the shiny apps are
source(file.path(paste0(getwd(), "/rscripts/EstimateRAmended2.R")))
source(file.path(paste0(getwd(), "/rscripts/coarse2estim.R")))

## Since the mcmc simulation from coarseDataTool takes quite some time, it will be good to save
## the result as a global variable, so here is to initiate the object
CDTfit <- NULL

simulateData <- function(n = 200, noise = 2, par1 = 6, par2 = 2) {
  g <- rgamma(n, par1, par2)
  EL <- -floor(runif(n, 0, noise))
  ER <- ceiling(runif(n, 0, noise))
  SL <- floor(g)+1
  SR <- ceiling(g)+1
  data.frame(EL, ER, SL, SR, type=1)
}

data(Measles1861)
data(Flu1918)
data(Smallpox1972)
data(SARS2003)
data(Flu2009)
alldatasets <- list(Measles1861, Flu1918, Smallpox1972, SARS2003, Flu2009)
simulationDataset <- list(simulation = simulateData())

###############################################################################

# Define server logic required to plot various variables against mpg
shinyServer(function(input, output) {

			mydata <- reactive(Data(input))
			mydata1 <- reactive(Data1(input))
			mySI <- reactive(SI(input))
			myR <- reactive(EstimateR(input))
			myPlot <- reactive(plotR(input))

			# view Data
			output$viewData <- renderTable({
						if ( (as.numeric(input$dataset)==0) && (is.null(input$DataFile) )){
							return(NULL)
						}else
						{
							if ( (as.numeric(input$dataset)==0) && (!is.null(input$DataFile) )){
								dataset <- read.csv(input$DataFile$datapath, header=TRUE)
							}else
							{
								dataset <- alldatasets[[as.numeric(input$dataset)]]
							}
							Time <- 1:length(dataset$Incidence)
							Incidence <- dataset$Incidence
							data.frame(Time, Incidence)
						}
					}, include.rownames=FALSE)

			output$viewData1 <- renderTable({
			  if ( (as.numeric(input$dataset1)==0) && (is.null(input$DataFile1) )){
			    return(NULL)
			  }else
			  {
			    if ( (as.numeric(input$dataset1)==0) && (!is.null(input$DataFile1) )){
			      dataset <- read.csv(input$DataFile1$datapath, header=TRUE)
			    } else if (as.numeric(input$dataset1)==2)  {
			      dataset <- simulationDataset$simulation
			    }
			    dataset
			  }
			}, include.rownames=FALSE)

			# Estimate R and plot results
			output$plot <- renderPlot({
						if ( (as.numeric(input$dataset)==0) && (is.null(input$DataFile) )){
							return(NULL)
						}else
						{
							if ( (as.numeric(input$dataset)==0) && (!is.null(input$DataFile) )){
								dataset <- read.csv(input$DataFile$datapath, header=TRUE)
							}else
							{
								dataset <- alldatasets[[as.numeric(input$dataset)]]
							}

						  if (!is.null(CDTfit)) {
						    EpiCDT(I = dataset$Incidence, T.Start=2:(length(dataset$Incidence)-as.numeric(input$width)+1),
						                T.End = (1+as.numeric(input$width)):length(dataset$Incidence), CDT = CDTfit,
						                n2 = input$n2SI, plot = TRUE)
						  } else {
						    if(input$whichR==1)
						    {
						      if(input$method1==1)
						      {
						        EstimateRAmended2(dataset$Incidence, T.Start=2:(length(dataset$Incidence)-as.numeric(input$width)+1),T.End=(1+as.numeric(input$width)):length(dataset$Incidence),method="ParametricSI", Mean.SI=as.numeric(input$meanSI), Std.SI=as.numeric(input$sdSI), Mean.Prior=as.numeric(input$meanPrior), Std.Prior=as.numeric(input$sdPrior), plot=TRUE)
						      }
						      if(input$method1==2)
						      {
						        if(is.null(input$SIFile))
						        {
						          return(NULL)
						        }else
						        {
						          SI.Distr <- read.csv(input$SIFile$datapath, header=TRUE)[,2]
						          EstimateRAmended2(dataset$Incidence, T.Start=2:(length(dataset$Incidence)-as.numeric(input$width)+1),T.End=(1+as.numeric(input$width)):length(dataset$Incidence),method="NonParametricSI", SI.Distr=SI.Distr, Mean.Prior=as.numeric(input$meanPrior), Std.Prior=as.numeric(input$sdPrior), plot=TRUE)

						        }
						      }
						      if(input$method1==3)
						      {
						        EstimateRAmended2(dataset$Incidence, T.Start=2:(length(dataset$Incidence)-as.numeric(input$width)+1),T.End=(1+as.numeric(input$width)):length(dataset$Incidence), method="UncertainSI", Mean.SI=as.numeric(input$meanmeanSI), Std.Mean.SI=as.numeric(input$stdmeanSI), Min.Mean.SI=as.numeric(input$minmeanSI), Max.Mean.SI=as.numeric(input$maxmeanSI), Std.SI=as.numeric(input$meanstdSI), Std.Std.SI=as.numeric(input$stdstdSI), Min.Std.SI=as.numeric(input$minstdSI), Max.Std.SI=as.numeric(input$maxstdSI), n1=input$n1SI, n2=input$n2SI, Mean.Prior=as.numeric(input$meanPrior), Std.Prior=as.numeric(input$sdPrior), plot=TRUE)
						      }
						    }else
						      if(input$whichR==2)
						      {
						        if(input$method2==1)
						        {
						          WT(dataset$Incidence, T.Start=2:(length(dataset$Incidence)-as.numeric(input$width)+1),T.End=(1+as.numeric(input$width)):length(dataset$Incidence),method="ParametricSI", Mean.SI=as.numeric(input$meanSI), Std.SI=as.numeric(input$sdSI), nSim=as.integer(input$nSim), plot=TRUE)
						        }
						        if(input$method2==2)
						        {
						          if(is.null(input$SIFile))
						          {
						            return(NULL)
						          }else
						          {
						            SI.Distr <- read.csv(input$SIFile$datapath, header=TRUE)[,2]
						            WT(dataset$Incidence, T.Start=2:(length(dataset$Incidence)-as.numeric(input$width)+1),T.End=(1+as.numeric(input$width)):length(dataset$Incidence),method="NonParametricSI", SI.Distr=SI.Distr, nSim=as.integer(input$nSim), plot=TRUE)
						          }
						        }
						      }
						  }
						}

					})

			# view table of results
			output$viewR <- renderTable({
						if ( (as.numeric(input$dataset)==0) && (is.null(input$DataFile) )){
							return(NULL)
						}else {
							if ( (as.numeric(input$dataset)==0) && (!is.null(input$DataFile) )){
								dataset <- read.csv(input$DataFile$datapath, header=TRUE)
							}else
							{
								dataset <- alldatasets[[as.numeric(input$dataset)]]
							}

						  if (!is.null(CDTfit)) {
						    R <- EpiCDT(I = dataset$Incidence, T.Start=2:(length(dataset$Incidence)-as.numeric(input$width)+1),
						                           T.End = (1+as.numeric(input$width)):length(dataset$Incidence), CDT = CDTfit,
						                n2 = input$n2SI)$R
						  } else {
						    if(input$whichR==1)
						    {
						      if(input$method1==1)
						      {
						        R <- EpiCDT(dataset$Incidence, T.Start=2:(length(dataset$Incidence)-as.numeric(input$width)+1),T.End=(1+as.numeric(input$width)):length(dataset$Incidence),method="ParametricSI", Mean.SI=as.numeric(input$meanSI), Std.SI=as.numeric(input$sdSI), Mean.Prior=as.numeric(input$meanPrior), Std.Prior=as.numeric(input$sdPrior), plot=FALSE)$R
						      }
						      if(input$method1==2)
						      {
						        if(is.null(input$SIFile))
						        {
						          return(NULL)
						        }else
						        {
						          SI.Distr <- read.csv(input$SIFile$datapath, header=TRUE)[,2]
						          R <- EpiCDT(dataset$Incidence, T.Start=2:(length(dataset$Incidence)-as.numeric(input$width)+1),T.End=(1+as.numeric(input$width)):length(dataset$Incidence),method="NonParametricSI", SI.Distr=SI.Distr, Mean.Prior=as.numeric(input$meanPrior), Std.Prior=as.numeric(input$sdPrior), plot=FALSE)$R

						        }
						      }
						      if(input$method1==3)
						      {
						        R <- EpiCDT(dataset$Incidence, T.Start=2:(length(dataset$Incidence)-as.numeric(input$width)+1),T.End=(1+as.numeric(input$width)):length(dataset$Incidence), method="UncertainSI", Mean.SI=as.numeric(input$meanmeanSI), Std.Mean.SI=as.numeric(input$stdmeanSI), Min.Mean.SI=as.numeric(input$minmeanSI), Max.Mean.SI=as.numeric(input$maxmeanSI), Std.SI=as.numeric(input$meanstdSI), Std.Std.SI=as.numeric(input$stdstdSI), Min.Std.SI=as.numeric(input$minstdSI), Max.Std.SI=as.numeric(input$maxstdSI), n1=input$n1SI, n2=input$n2SI, Mean.Prior=as.numeric(input$meanPrior), Std.Prior=as.numeric(input$sdPrior), plot=FALSE)$R
						      }
						    }else if(input$whichR==2)
						    {
						      if(input$method2==1)
						      {
						        R <- WT(dataset$Incidence, T.Start=2:(length(dataset$Incidence)-as.numeric(input$width)+1),T.End=(1+as.numeric(input$width)):length(dataset$Incidence),method="ParametricSI", Mean.SI=as.numeric(input$meanSI), Std.SI=as.numeric(input$sdSI), nSim=as.integer(input$nSim), plot=FALSE)$R
						      }
						      if(input$method2==2)
						      {
						        if(is.null(input$SIFile))
						        {
						          return(NULL)
						        }else
						        {
						          SI.Distr <- read.csv(input$SIFile$datapath, header=TRUE)[,2]
						          R <- WT(dataset$Incidence, T.Start=2:(length(dataset$Incidence)-as.numeric(input$width)+1),T.End=(1+as.numeric(input$width)):length(dataset$Incidence),method="NonParametricSI", SI.Distr=SI.Distr, nSim=as.integer(input$nSim), plot=FALSE)$R
						        }
						      }
						  }


							}
							R
						}

					}, include.rownames=FALSE)

			output$viewSI <- renderTable({
						if ( (as.numeric(input$dataset)==0) && (is.null(input$DataFile) )){
							return(NULL)
						}else
						{
							if ( (as.numeric(input$dataset)==0) && (!is.null(input$DataFile) )){
								dataset <- read.csv(input$DataFile$datapath, header=TRUE)
							}else
							{
								dataset <- alldatasets[[as.numeric(input$dataset)]]
							}

						  if (!is.null(CDTfit)) {
						    SI <- EpiCDT(I = dataset$Incidence, T.Start=2:(length(dataset$Incidence)-as.numeric(input$width)+1),
						                 T.End = (1+as.numeric(input$width)):length(dataset$Incidence), CDT = CDTfit,
						                 n2 = input$n2SI)$SIDistr
						    return(SI)
						  } else {
						    if(input$whichR==1)
						    {
						      if(input$method1==1)
						      {
						        SI <- EpiCDT(dataset$Incidence, T.Start=2:(length(dataset$Incidence)-as.numeric(input$width)+1),T.End=(1+as.numeric(input$width)):length(dataset$Incidence),method="ParametricSI", Mean.SI=as.numeric(input$meanSI), Std.SI=as.numeric(input$sdSI), Mean.Prior=as.numeric(input$meanPrior), Std.Prior=as.numeric(input$sdPrior), plot=FALSE)$SIDistr
						        return(SI)
						      }
						      if(input$method1==2)
						      {
						        if(is.null(input$SIFile))
						        {
						          return(NULL)
						        }else
						        {
						          SI.Distr <- read.csv(input$SIFile$datapath, header=TRUE)[,2]
						          SI <- EpiCDT(dataset$Incidence, T.Start=2:(length(dataset$Incidence)-as.numeric(input$width)+1),T.End=(1+as.numeric(input$width)):length(dataset$Incidence),method="NonParametricSI", SI.Distr=SI.Distr, Mean.Prior=as.numeric(input$meanPrior), Std.Prior=as.numeric(input$sdPrior), plot=FALSE)$SIDistr
						          return(SI)
						        }
						      }
						      if(input$method1==3)
						      {
						        return(NULL)
						      }
						    }else
						    {
						      if(input$whichR==2)
						      {
						        if(input$method2==1)
						        {
						          SI <- WT(dataset$Incidence, T.Start=2:(length(dataset$Incidence)-as.numeric(input$width)+1),T.End=(1+as.numeric(input$width)):length(dataset$Incidence),method="ParametricSI", Mean.SI=as.numeric(input$meanSI), Std.SI=as.numeric(input$sdSI), nSim=2, plot=FALSE)$SIDistr
						          return(SI)
						        }
						        if(input$method2==2)
						        {
						          if(is.null(input$SIFile))
						          {
						            return(NULL)
						          }else
						          {
						            SI.Distr <- read.csv(input$SIFile$datapath, header=TRUE)[,2]
						            SI <- WT(dataset$Incidence, T.Start=2:(length(dataset$Incidence)-as.numeric(input$width)+1),T.End=(1+as.numeric(input$width)):length(dataset$Incidence),method="NonParametricSI", SI.Distr=SI.Distr, nSim=2, plot=FALSE)$SIDistr
						            return(SI)
						          }
						        }
						      }
						    }
						  }



						}
					}, include.rownames=FALSE)

			output$viewCDT <- renderTable({
			  if ( (as.numeric(input$dataset1)==0) && (is.null(input$DataFile1))){
			    return(NULL)
			  } else if ( (as.numeric(input$dataset1)==0) && (!is.null(input$DataFile1) )){
			      dataset <- read.csv(input$DataFile1$datapath, header=TRUE)
			    } else {
			      dataset <- simulationDataset$simulation
			    }

			  CDTfit <<- dic.fit.mcmc(dat = dataset, dist = "off1G", n.samples = 1000)
			  CDTfit@ests
			}, include.rownames = TRUE)

			output$downloadData <- downloadHandler(
					filename = function() { paste('Data-', input$dataset, '.csv', sep='') },
					content = function(file) {
						if ( !is.null(mydata) ){write.csv(mydata(), file, row.names=F)}
					}
			)

			output$downloadSI <- downloadHandler(
					filename = function() { paste('SI-Distr-Mean-', input$meanSI, '-Std-',input$sdSI,'.csv', sep='') },
					content = function(file) {
						if ( !is.null(mydata) ){write.csv(mySI(), file, row.names=F)}
					}
			)

			output$downloadR <- downloadHandler(
					filename = function() { paste('R-estimates-', input$dataset, '.csv', sep='') },
					content = function(file) {
						if ( !is.null(mydata) ){write.csv(myR(), file, row.names=F)}

					}
			)

			output$savePlot <- downloadHandler(
					filename = function() { paste('R-estimates-', input$dataset, '.png', sep='') },
					content = function(file) {
						png(file)
						print(myPlot())
						dev.off()
					}
			)


		})

###############################################################################

Data <- function(x)
{
	if ( (as.numeric(x$dataset)==0) && (is.null(x$DataFile) )){
		return(NULL)
	}else
	{
		if ( (as.numeric(x$dataset)==0) && (!is.null(x$DataFile) )){
			dataset <- read.csv(x$DataFile$datapath, header=TRUE)
		}else
		{
			dataset <- alldatasets[[as.numeric(x$dataset)]]
		}
		Time <- 1:length(dataset$Incidence)
		Incidence <- dataset$Incidence
		return(data.frame(Time, Incidence))
	}
}

Data1 <- function(x)
{
  if ( (as.numeric(x$dataset1)==0) && (is.null(x$DataFile1) )){
    return(NULL)
  }else
  {
    if ( (as.numeric(x$dataset1)==0) && (!is.null(x$DataFile1) )){
      dataset <- read.csv(x$DataFile1$datapath, header=TRUE)
    }else
    {
      dataset <- alldatasets[[as.numeric(x$dataset1)]]
    }
    Time <- 1:length(dataset$Incidence)
    Incidence <- dataset$Incidence
    return(data.frame(Time, Incidence))
  }
}


SI <- function(x)
{
	if ( (as.numeric(x$dataset)==0) && (is.null(x$DataFile) )){
		return(NULL)
	}else
	{
		if ( (as.numeric(x$dataset)==0) && (!is.null(x$DataFile) )){
			dataset <- read.csv(x$DataFile$datapath, header=TRUE)
		}else
		{
			dataset <- alldatasets[[as.numeric(x$dataset)]]
		}
		T.Start <- 2:(length(dataset$Incidence)-as.numeric(x$width)+1)
		T.End <- (1+as.numeric(x$width)):length(dataset$Incidence)
		return(EpiCDT(dataset$Incidence, T.Start=2:(length(dataset$Incidence)-as.numeric(x$width)+1),T.End=(1+as.numeric(x$width)):length(dataset$Incidence),method="ParametricSI", Mean.SI=as.numeric(x$meanSI), Std.SI=as.numeric(x$sdSI), Mean.Prior=as.numeric(x$meanPrior), Std.Prior=as.numeric(x$sdPrior), plot=TRUE)$SIDistr)
	}
}

EstimateR <- function(x)
{
	if ( (as.numeric(x$dataset)==0) && (is.null(x$DataFile) )){
		return(NULL)
	}else
	{
		if ( (as.numeric(x$dataset)==0) && (!is.null(x$DataFile) )){
			dataset <- read.csv(x$DataFile$datapath, header=TRUE)
		}else
		{
			dataset <- alldatasets[[as.numeric(x$dataset)]]
		}
		if(x$whichR==1)
		{
			T.Start <- 2:(length(dataset$Incidence)-as.numeric(x$width)+1)
			T.End <- (1+as.numeric(x$width)):length(dataset$Incidence)
			R <- EpiCDT(dataset$Incidence, T.Start=2:(length(dataset$Incidence)-as.numeric(x$width)+1),T.End=(1+as.numeric(x$width)):length(dataset$Incidence),method="ParametricSI", Mean.SI=as.numeric(x$meanSI), Std.SI=as.numeric(x$sdSI), Mean.Prior=as.numeric(x$meanPrior), Std.Prior=as.numeric(x$sdPrior), plot=TRUE)$R
			R.Mean <- R$'Mean(R)'
			R.Std <- R$'Std(R)'
			R.Median <- R$'Median(R)'
			R.Q.0.025 <- R$'Quantile.0.025(R)'
			R.Q.0.05 <- R$'Quantile.0.05(R)'
			R.Q.0.25 <- R$'Quantile.0.25(R)'
			R.Q.0.75 <- R$'Quantile.0.75(R)'
			R.Q.0.95 <- R$'Quantile.0.95(R)'
			R.Q.0.975 <- R$'Quantile.0.975(R)'
			R <- data.frame(T.Start, T.End,R.Mean,R.Std,R.Q.0.025,R.Q.0.05,R.Q.0.25,R.Median,R.Q.0.75,R.Q.0.95,R.Q.0.975)
			names(R) <- c("T.Start","T.End","Mean(R)","Std(R)","Quantile.0.025(R)","Quantile.0.05(R)","Quantile.0.25(R)","Median(R)","Quantile.0.75(R)","Quantile.0.95(R)","Quantile.0.975(R)")
		}
		if(x$whichR==2)
		{
			R <- WT(dataset$Incidence, T.Start=2:(length(dataset$Incidence)-as.numeric(x$width)+1),T.End=(1+as.numeric(x$width)):length(dataset$Incidence),method="ParametricSI", Mean.SI=as.numeric(x$meanSI), Std.SI=as.numeric(x$sdSI), nSim=as.integer(x$nSim), plot=TRUE)$R
		}
		return(R)
	}
}

plotR <- function(x)
{
	if ( (as.numeric(x$dataset)==0) && (is.null(x$DataFile) )){
		plot(1,1,col="white",main="",xlab="",ylab="",axes=FALSE)
	}else
	{
		if ( (as.numeric(x$dataset)==0) && (!is.null(x$DataFile) )){
			dataset <- read.csv(x$DataFile$datapath, header=TRUE)
		}else
		{
			dataset <- alldatasets[[as.numeric(x$dataset)]]
		}
		if(x$whichR==1)
		{
			EpiCDT(dataset$Incidence, T.Start=2:(length(dataset$Incidence)-as.numeric(x$width)+1),T.End=(1+as.numeric(x$width)):length(dataset$Incidence),method="ParametricSI", Mean.SI=as.numeric(x$meanSI), Std.SI=as.numeric(x$sdSI), Mean.Prior=as.numeric(x$meanPrior), Std.Prior=as.numeric(x$sdPrior), plot=TRUE)
		}else
		if(x$whichR==2)
		{
			WT(dataset$Incidence, T.Start=2:(length(dataset$Incidence)-as.numeric(x$width)+1),T.End=(1+as.numeric(x$width)):length(dataset$Incidence),method="ParametricSI", Mean.SI=as.numeric(x$meanSI), Std.SI=as.numeric(x$sdSI), nSim=as.integer(x$nSim), plot=TRUE)
		}
	}
}

