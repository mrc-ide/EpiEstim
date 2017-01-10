rm(list=ls())

library(coarseDataTools)

devtools::document()
devtools::check() 
devtools::build()
devtools::install()

# check simplest version of EpiEstim works well 

data("Flu2009")
R <- EstimateR(Flu2009$Incidence, T.Start=2:26, T.End=8:32, method="NonParametricSI", 
                        SI.Distr=Flu2009$SI.Distr, plot=TRUE, leg.pos=xy.coords(1,3))
R

#source("~/Documents/HackoutCode/dic.fit.R")
#source("~/Documents/HackoutCode/dic.fit.mcmc.R")

serialIntervalData <-  read.csv("~/Dropbox/Hackout code/datasets/RotavirusEcuadorSIData.csv", header = FALSE)

names <- c("EL", "ER", "SL", "SR", "type")
colnames(serialIntervalData) <- names

serialIntervalData <- as.data.frame(serialIntervalData)

####  FEED INTO COARSE DATA TOOLS

fit <- dic.fit.mcmc(dat = serialIntervalData, dist="G") # may take a few minutes

casesPerDayData <- read.csv("~/Dropbox/Hackout code/datasets/GermanyRotavirus1516.csv", header = FALSE)
casesPerDayData <- t(casesPerDayData)
casesPerDayData <- cbind(1:length(casesPerDayData), casesPerDayData)
casesPerDayData <- as.data.frame(casesPerDayData)
names <- c("Time", "Cases")
colnames(casesPerDayData) <- names

####  FEED INTO EPIESTIM

set.seed(1)
res <- EpiEstim::EstimateR(casesPerDayData[1:53,2], T.Start=2:43, T.End=9:50, method="NonParametricUncertainSI", n2 = dim(fit@samples)[2], CDT = fit, plot=TRUE)
# compare with version without uncertainty
res0 <- EpiEstim::EstimateR(casesPerDayData[1:53,2], T.Start=2:43, T.End=9:50, method="ParametricSI", Mean.SI = 2.5, Std.SI = 1.0, n2 = dim(fit@samples)[2], plot=TRUE)

