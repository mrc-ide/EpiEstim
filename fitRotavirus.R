rm(list=ls())

#devtools::install_github("nickreich/coarseDataTools", ref = "hackout3", force=TRUE)
#library(coarseDataTools)
devtools::document()
devtools::check() 
devtools::build()
devtools::install()

# check simplest version of EpiEstim works well 

data("Flu2009")
R <- EstimateR(Flu2009$Incidence, T.Start=2:26, T.End=8:32, method="NonParametricSI", 
               SI.Distr=Flu2009$SI.Distr, plot=TRUE, leg.pos=xy.coords(1,3))

Rc <- WT(Flu2009$Incidence, T.Start=2:26, T.End=8:32, method="NonParametricSI", 
               SI.Distr=Flu2009$SI.Distr, plot=TRUE, leg.pos=xy.coords(1,3))

# now check the new SIFromData method using the mock rotavirus data 

data("MockRotavirus")

## estimate the reproduction number (method "SIFromData")
set.seed(1)
R_SIFromData <- EstimateR(MockRotavirus$Incidence, 
                          T.Start=2:47, T.End=8:53, 
                          method="SIFromData", SI.Data=MockRotavirus$SI.Data, 
                          SI.parametricDistr = "G", 
                          MCMC.control = list(burnin = 1000, thin=10), 
                          n1 = 500, n2 = 50,
                          plot=TRUE, leg.pos=xy.coords(1,3))

# check convergence of MCMC
# R_SIFromData$MCMC_converged
# graphically check convergence of MCMC
# par(mfrow=c(2,1))
# plot(R_SIFromData$SI.Moments[,"Mean"], type="l", xlab="Iterations", ylab="Mean SI")
# plot(R_SIFromData$SI.Moments[,"Std"], type="l", xlab="Iterations", ylab="Std SI")

# compare with version with no uncertainty

R_Parametric <- EstimateR(MockRotavirus$Incidence, 
                          T.Start=2:47, T.End=8:53, 
                          method="ParametricSI", 
                          Mean.SI = mean(R_SIFromData$SI.Moments$Mean), Std.SI = mean(R_SIFromData$SI.Moments$Std), 
                          plot=TRUE)

p_uncertainty <- plots(R_SIFromData, "R", ylim=c(0,1.5))
p_no_uncertainty <- plots(R_Parametric, "R", ylim=c(0,1.5))
gridExtra::grid.arrange(p_uncertainty, p_no_uncertainty,ncol=2)

# now check the new SIFromSample method using the mock rotavirus data 

set.seed(1)
SI.fit <- coarseDataTools::dic.fit.mcmc(dat = MockRotavirus$SI.Data, dist="G", burnin = 1000, n.samples = 5000, init.pars=init_MCMC_params(MockRotavirus$SI.Data, "G")) # may take a few minutes
SI.Sample <- coarse2estim(SI.fit, thin=10)$SI.Sample

R_SIFromSample <- EstimateR(MockRotavirus$Incidence, 
                            T.Start=2:47, T.End=8:53, 
                            method="SIFromSample", SI.Sample=SI.Sample,
                            n2 = 50,
                            plot=TRUE, leg.pos=xy.coords(1,3))

# check that R_SIFromSample is the same as R_SIFromData 
# since they were generated using the same MCMC algorithm to generate the SI sample
# (either internally to EpiEstim or externally)
all(R_SIFromSample$R$`Mean(R)` == R_SIFromData$R$`Mean(R)`) ### note this could be turned into a test :-)

# Run this on the Ebola data [for now all countries together]
ebola.SI.Data <- readRDS("~/Dropbox/EbolaDataForEpiEstim/EbolaSIData.RDS")
ebola.SI.Data.early <- readRDS("~/Dropbox/EbolaDataForEpiEstim/EbolaSIData_early.RDS")
ebola.Incid.Data <- readRDS("~/Dropbox/EbolaDataForEpiEstim/EbolaIncidenceData.RDS")
ebola.Incid.Data.early <- ebola.Incid.Data[ebola.Incid.Data$Date<=as.Date("2014-03-31", format="%Y-%m-%d"),]

ebola <- list(Incidence = ebola.Incid.Data$Incid, SI.Data = ebola.SI.Data)
ebola_early <- list(Incidence = ebola.Incid.Data.early$Incid, SI.Data = ebola.SI.Data.early)

## estimates from whole incidence but with early SI estimates

set.seed(1)
R_ebola_SIFromData_early <- EstimateR(ebola$Incidence, 
                          T.Start=2:477, T.End=8:483, 
                          method="SIFromData", SI.Data=ebola_early$SI.Data, 
                          SI.parametricDistr = "G", 
                          MCMC.control = list(burnin = 1000, thin=10), 
                          n1 = 500, n2 = 50,
                          plot=TRUE, leg.pos=xy.coords(1,3))

# what if instead we had neglected the uncertain point? 
tmp <- ebola.SI.Data.early[ebola.SI.Data.early$type==2,]
meanSI_notype1 <- mean(tmp$SL-tmp$EL) # compare with ebola_SI_moments_early['Mean']
stdSI_notype1 <- sd(tmp$SL-tmp$EL) # compare with ebola_SI_moments_early['Std']

R_ebola_ParametricSI_early_notype1 <- EstimateR(ebola$Incidence, 
                                      T.Start=2:477, T.End=8:483, 
                                      method="ParametricSI", Mean.SI = meanSI_notype1, Std.SI = stdSI_notype1, 
                                      plot=TRUE, leg.pos=xy.coords(1,3))

# what if instead we had taken the midpoint of the uncertain points?
tmp <- ebola.SI.Data.early
tmp$Emid <- (tmp$EL+tmp$ER)/2
meanSI_midpoint <- mean(tmp$SL-tmp$Emid) # compare with ebola_SI_moments_early['Mean']
stdSI_midpoint <- sd(tmp$SL-tmp$Emid) # compare with ebola_SI_moments_early['Std']

R_ebola_ParametricSI_early_midpoint <- EstimateR(ebola$Incidence, 
                                                T.Start=2:477, T.End=8:483, 
                                                method="ParametricSI", Mean.SI = meanSI_midpoint, Std.SI = stdSI_midpoint, 
                                                plot=TRUE, leg.pos=xy.coords(1,3))

p1 <- plots(R_ebola_SIFromData_early, "R", ylim=c(0, 11))
p2 <- plots(R_ebola_ParametricSI_early_notype1, "R", ylim=c(0, 11))
p3 <- plots(R_ebola_ParametricSI_early_midpoint, "R", ylim=c(0, 11))
gridExtra::grid.arrange(p1, p2, p3, ncol=1)

rel_error_std_notype1 <- (R_ebola_ParametricSI_early_notype1$R$`Std(R)` / R_ebola_SIFromData_early$R$`Std(R)`)
rel_error_std_midpoint <- (R_ebola_ParametricSI_early_midpoint$R$`Std(R)` / R_ebola_SIFromData_early$R$`Std(R)`)

plot(rel_error_std_notype1, type="l")
lines(rel_error_std_midpoint, col="blue")

plot(R_ebola_SIFromData_early$R$`Quantile.0.975(R)`, type="l", ylim=c(0,11))
lines(R_ebola_SIFromData_early$R$`Quantile.0.025(R)`)
lines(R_ebola_ParametricSI_early_notype1$R$`Quantile.0.975(R)`, col="blue")
lines(R_ebola_ParametricSI_early_notype1$R$`Quantile.0.025(R)`, col="blue")
abline(h=1,col="red")
lines(rel_error_std_notype1, col="green")

## estimates from whole incidence but with late SI estimates

set.seed(1)
R_ebola_SIFromData <- EstimateR(ebola$Incidence, 
                                      T.Start=2:477, T.End=8:483, 
                                      method="SIFromData", SI.Data=ebola$SI.Data, 
                                      SI.parametricDistr = "G", 
                                      MCMC.control = list(burnin = 1000, thin=10), 
                                      n1 = 500, n2 = 50,
                                      plot=TRUE, leg.pos=xy.coords(1,3))

# what if instead we had neglected the uncertain point? 
tmp <- ebola.SI.Data[ebola.SI.Data$type==2,]
meanSI_notype1 <- mean(tmp$SL-tmp$EL) # compare with ebola_SI_moments['Mean']
stdSI_notype1 <- sd(tmp$SL-tmp$EL) # compare with ebola_SI_moments['Std']

R_ebola_ParametricSI_notype1 <- EstimateR(ebola$Incidence, 
                                                T.Start=2:477, T.End=8:483, 
                                                method="ParametricSI", Mean.SI = meanSI_notype1, Std.SI = stdSI_notype1, 
                                                plot=TRUE, leg.pos=xy.coords(1,3))

# what if instead we had taken the midpoint of the uncertain points?
tmp <- ebola.SI.Data
tmp$Emid <- (tmp$EL+tmp$ER)/2
meanSI_midpoint <- mean(tmp$SL-tmp$Emid) # compare with ebola_SI_moments['Mean']
stdSI_midpoint <- sd(tmp$SL-tmp$Emid) # compare with ebola_SI_moments['Std']

R_ebola_ParametricSI_midpoint <- EstimateR(ebola$Incidence, 
                                                 T.Start=2:477, T.End=8:483, 
                                                 method="ParametricSI", Mean.SI = meanSI_midpoint, Std.SI = stdSI_midpoint, 
                                                 plot=TRUE, leg.pos=xy.coords(1,3))

p1 <- plots(R_ebola_SIFromData, "R", ylim=c(0, 11))
p2 <- plots(R_ebola_ParametricSI_notype1, "R", ylim=c(0, 11))
p3 <- plots(R_ebola_ParametricSI_midpoint, "R", ylim=c(0, 11))
gridExtra::grid.arrange(p1, p2, p3, ncol=1)

rel_error_std_notype1 <- (R_ebola_ParametricSI_notype1$R$`Std(R)` / R_ebola_SIFromData$R$`Std(R)`)
rel_error_std_midpoint <- (R_ebola_ParametricSI_midpoint$R$`Std(R)` / R_ebola_SIFromData$R$`Std(R)`)

plot(rel_error_std_notype1, type="l")
lines(rel_error_std_midpoint, col="blue")


p1 <- plots(R_ebola_SIFromData_early, "R", ylim=c(0, 5))
p2 <- plots(R_ebola_SIFromData, "R", ylim=c(0, 5))
gridExtra::grid.arrange(p1, p2,ncol=1)



