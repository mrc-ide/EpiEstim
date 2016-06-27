\name{Measles1861}
     \docType{data}
     \alias{Measles1861}
     \title{Data on the 1861 measles epidemic in Hagelloch, Germany}
     \description{
       This data set gives 
	   
	   1/ the daily incidence of onset of symptoms in Hallegoch (Germany) during the 1861 measles epidemic (see source and references),
	   
	   2/ the discrete daily distribution of the serial interval for measles, assuming a shifted Gamma distribution with mean 14.9 days, standard deviation 3.9 days and shift 1 day (see references).
	   
     }
     \usage{Measles1861}
     \format{
		A list of two elements: 
	 
		1/ "Incidence": a vector containing 48 days of observation,
		
		2/ "SI.Distr": a vector containing a set of 38 probabilities.
		}
		\source{Groendyke C. et al. (2011) Bayesian Inference for Contact Networks Given Epidemic Data. Scandinavian Journal of Statistics 38(3): 600-616.}
     \references{
       Groendyke C. et al. (2011) Bayesian Inference for Contact Networks Given Epidemic Data. Scandinavian Journal of Statistics 38(3): 600-616.
     }
	 \examples{
		## load data on measles in Hallegoch in 1861
		data("Measles1861")

		## estimate the reproduction number (method "NonParametricSI")
		EstimateR(Measles1861$Incidence, 
		  T.Start=17:42, T.End=23:48, 
		  method="NonParametricSI", SI.Distr=Measles1861$SI.Distr, 
		  plot=TRUE, leg.pos=xy.coords(1,7))
		# the second plot produced shows, at each each day, 
		# the estimate of the reproduction number 
		# over the 7-day window finishing on that day.

	}
     \keyword{datasets}