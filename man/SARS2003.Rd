\name{SARS2003}
     \docType{data}
     \alias{SARS2003}
     \title{Data on the 2003 SARS epidemic in Hong Kong}
     \description{
       This data set gives 
	   
	   1/ the daily incidence of onset of symptoms in Hong Kong during the 2003 severe acute respiratory syndrome (SARS) epidemic (see source and references),
	   
	   2/ the discrete daily distribution of the serial interval for SARS, assuming a shifted Gamma distribution with mean 8.4 days, standard deviation 3.8 days and shift 1 day (see references).
	   
     }
     \usage{SARS2003}
     \format{
		A list of two elements: 
	 
		1/ "Incidence": a vector containing 107 days of observation,
		
		2/ "SI.Distr": a vector containing a set of 25 probabilities.
		}
		\source{Cori A. et al. (2009) Temporal variability and social heterogeneity in disease transmission: the case of SARS in Hong Kong. PLoS Comput Biol 5(8): e1000471.}
     \references{
       Cori A. et al. (2009) Temporal variability and social heterogeneity in disease transmission: the case of SARS in Hong Kong. PLoS Comput Biol 5(8): e1000471.
	   
	   Lipsitch M. et al. (2003) Transmission dynamics and control of severe acute respiratory syndrome. Science 300(5627): 1966-1970.
     }
	 \examples{
		## load data on SARS in Hong Kong in 2003
		data("SARS2003")

		## estimate the reproduction number (method "NonParametricSI")
		EstimateR(SARS2003$Incidence, T.Start=14:101, T.End=20:107, method="NonParametricSI", 
			SI.Distr=SARS2003$SI.Distr, plot=TRUE, leg.pos=xy.coords(1,7))
		# the second plot produced shows, at each each day, 
		# the estimate of the reproduction number over the 7-day window finishing on that day.

	}
     \keyword{datasets}