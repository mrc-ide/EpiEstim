\name{Smallpox1972}
     \docType{data}
     \alias{Smallpox1972}
     \title{Data on the 1972 smallpox epidemic in Kosovo}
     \description{
       This data set gives 
	   
	   1/ the daily incidence of onset of symptoms in Kosovo during the 1972 smallpox epidemic (see source and references),
	   
	   2/ the discrete daily distribution of the serial interval for smallpox, assuming a shifted Gamma distribution with mean 22.4 days, standard deviation 6.1 days and shift 1 day (see references).
	   
     }
     \usage{Smallpox1972}
     \format{
		A list of two elements: 
	 
		1/ "Incidence": a vector containing 57 days of observation,
		
		2/ "SI.Distr": a vector containing a set of 46 probabilities.
		}
		\source{Fenner F. et al. (1988) Smallpox and its Eradication. Geneva, World Health Organization.}
     \references{
       Fenner F. et al. (1988) Smallpox and its Eradication. Geneva, World Health Organization.
	   
	   Gani R. and S. Leach (2001) Transmission potential of smallpox in contemporary populations. Nature 414(6865): 748-751.
	   
	   Riley S. and N. M. Ferguson (2006) Smallpox transmission and control: spatial dynamics in Great Britain. Proc Natl Acad Sci U S A 103(33): 12637-12642.
     }
	 \examples{
		## load data on smallpox in Kosovo in 1972
		data("Smallpox1972")

		## estimate the reproduction number (method "NonParametricSI")
		EstimateR(Smallpox1972$Incidence, T.Start=27:51, T.End=33:57, method="NonParametricSI", 
			SI.Distr=Smallpox1972$SI.Distr, plot=TRUE, leg.pos=xy.coords(1,15))
		# the second plot produced shows, at each each day, 
		# the estimate of the reproduction number over the 7-day window finishing on that day.

	}
     \keyword{datasets}