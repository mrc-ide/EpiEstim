rm(list=ls())

install_github("nickreich/coarseDataTools", ref = "hackout3")

devtools::document()
devtools::check() 
devtools::build()
devtools::install()

# check simplest version of EpiEstim works well 

data("Flu2009")
R <- EstimateR(Flu2009$Incidence, T.Start=2:26, T.End=8:32, method="NonParametricSI", 
               SI.Distr=Flu2009$SI.Distr, plot=TRUE, leg.pos=xy.coords(1,3))
R

# now check the new SIFromData method using the mock rotavirus data 

data("MockRotavirus")

## estimate the reproduction number (method "SIFromData")
set.seed(1)
R_SIFromData <- EstimateR(MockRotavirus$Incidence, 
                          T.Start=2:47, T.End=8:53, 
                          method="SIFromData", SI.Data=MockRotavirus$SI.Data, SI.parametricDistr = "G", MCMC.burnin = 1000, 
                          n1 = 1000, n2 = 50,
                          plot=TRUE, leg.pos=xy.coords(1,3))

SI.Data <- MockRotavirus$SI.Data
SI.parametricDistr <- "off1G"
dic.fit.mcmc(dat = SI.Data, dist=SI.parametricDistr, burnin = MCMC.burnin, n.samples = n1, init.pars=init_MCMC_params(SI.Data, SI.parametricDistr))

### WHERE IS THIS FUNCTION???
### dist.optim.untransform

# check convergence of MCMC
# par(mfrow=c(2,1))
# plot(R_SIFromData$SIDistr[,"Mean.SI.sample"], type="l", xlab="Iterations", ylab="Mean SI")
# plot(R_SIFromData$SIDistr[,"Std.SI.sample"], type="l", xlab="Iterations", ylab="Std SI")

# compare with version with no uncertainty

R_Parametric <- EstimateR(MockRotavirus$Incidence, 
                          T.Start=2:47, T.End=8:53, 
                          method="ParametricSI", 
                          Mean.SI = mean(R_SIFromData$SIDistr$Mean.SI.sample), Std.SI = mean(R_SIFromData$SIDistr$Std.SI.sample), 
                          plot=TRUE)

p_uncertainty <- plots(R_SIFromData, "R")
p_no_uncertainty <- plots(R_Parametric, "R")
gridExtra::grid.arrange(p_uncertainty, p_no_uncertainty,ncol=2)

# now check the new SIFromSample method using the mock rotavirus data 

set.seed(1)

dist <- "G"
SI.Data <- MockRotavirus$SI.Data

SI.fit <- coarseDataTools::dic.fit.mcmc(dat = MockRotavirus$SI.Data, dist="G", burnin = 1, n.samples = 1000) # may take a few minutes
SI.fit2 <- coarseDataTools::dic.fit.mcmc(dat = MockRotavirus$SI.Data, dist="G",init.pars=init_MCMC_params(MockRotavirus$SI.Data, "G"), burnin = 1, n.samples = 1000) # may take a few minutes



SI.Sample <- coarse2estim(SI.fit, nrow(SI.fit@samples))$prob_matrix

R_SIFromSample <- EstimateR(MockRotavirus$Incidence, 
                            T.Start=2:47, T.End=8:53, 
                            method="SIFromSample", SI.Sample=SI.Sample,
                            n2 = 50,
                            plot=TRUE, leg.pos=xy.coords(1,3))

# check that R_SIFromSample is the same as R_SIFromData 
# since they were generated using the same MCMC algorithm to generate the SI sample
# (either internally to EpiEstim or externally)
all(R_SIFromSample$R$`Mean(R)` == R_SIFromData$R$`Mean(R)`) ### note this could be turned into a test :-)
