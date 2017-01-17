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

# now check the new version using the mock rotavirus data 

data("MockRotavirus")

## estimate the reproduction number (method "NonParametricUncertainSI")
R_NonParametricUncertainSI <- EstimateR(MockRotavirus$Incidence, 
          T.Start=2:47, T.End=8:53, 
          method="NonParametricUncertainSI", SI.Data=MockRotavirus$SI.Data, SI.parametricDistr = "G", MCMC.burnin = 5000, 
          n1 = 1000, n2 = 50,
          plot=TRUE, leg.pos=xy.coords(1,3))

#SI.fit <- dic.fit.mcmc(dat = MockRotavirus$SI.Data, dist="G", burnin = 5000, n.samples = 1000) # may take a few minutes

# compare with version with no uncertainty

R_Parametric <- EstimateR(MockRotavirus$Incidence, 
                          T.Start=2:47, T.End=8:53, 
                          method="ParametricSI", 
                          Mean.SI = mean(R_NonParametricUncertainSI$SIDistr$Mean.SI.sample), Std.SI = mean(R_NonParametricUncertainSI$SIDistr$Std.SI.sample), 
                          plot=TRUE)
  
R_NonParametricUncertainSI$R$`Mean(R)`
R_Parametric$R$`Mean(R)`

R_NonParametricUncertainSI$R$`Quantile.0.975(R)`
R_Parametric$R$`Quantile.0.975(R)`

p_uncertainty <- plots(R_NonParametricUncertainSI, "R")
p_no_uncertainty <- plots(R_Parametric, "R")
gridExtra::grid.arrange(p_uncertainty, p_no_uncertainty,ncol=2)

