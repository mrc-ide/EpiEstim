install.packages("devtools")
library(devtools)
install.packages("quantreg")
install.packages("coda")
install.packages("mcmc")
install.packages("reshape2")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("gridExtra")
install.packages("plotly")
library(MCMCpack)
library(reshape2)
library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)
library(plotly)
#library(EpiEstim)
install_github("nickreich/coarseDataTools", ref = "hackout3", force=TRUE)
library(coarseDataTools)


source("Z:/Hackout3/EpiEstim/R/EstimationR.R")
#source("~/Downloads/dic.fit.mcmc.R")

source("Z:/Hackout3/EpiEstim/R/stochasticSEIRModel3.R")

simulatedData <- simulateTransmissionModel(2)  # contains simulatedData$casesPerDay and simulatedData$dataMatrix
serialIntervalData <- na.omit(simulatedData$dataMatrix)
casesPerDayData <- simulatedData$casesPerDay

names <- c("EL", "ER", "SL", "SR", "type")
colnames(serialIntervalData) <- names
serialIntervalData[,1] = serialIntervalData[,2]
serialIntervalData[,2] = serialIntervalData[,3]
serialIntervalData[,3] = serialIntervalData[,4]
serialIntervalData[,4] = serialIntervalData[,5]
serialIntervalData[,5] = 0

serialIntervalData <- as.data.frame(serialIntervalData)

####  FEED INTO COARSE DATA TOOLS

# Only use 80 host pairs' interval data to estimate the serial interval
fit <- dic.fit.mcmc(dat = serialIntervalData[1:80,], dist="G", init.pars = c(0.5, 0.5))



####  FEED INTO EPIESTIM

est <- EstimateR(casesPerDayData[1:100,2], T.Start=5:29, T.End=16:40, n2 = dim(fit@samples)[2], CDT = fit, plot=TRUE)

#  We are getting an error here from the EstimateR code:  "SI.Distr does not sum to 1".  We haven't yet figured out the source of this error. 






#Flu2009$Incidence, T.Start=2:26, T.End=8:32, method = c("NonParametricUncertainSI"),
#SI.Dist.Matrix = c2e$prob_matrix, n2=small_n



# EpiCDT <- function(..., CDT = NULL) {
#   if (!is.null(CDT)) {
#     c2e <- coarse2estim(CDT)
#     EstimateR(..., method = c("NonParametricUncertainSI"),
#                       SI.Dist.Matrix = c2e$prob_matrix)
#   } else {
#     EstimateR(...)
#   }
# }

