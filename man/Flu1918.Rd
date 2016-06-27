\name{Flu1918}
     \docType{data}
     \alias{Flu1918}
     \title{Data on the 1918 H1N1 influenza pandemic in Baltimore}
     \description{
       This data set gives 
	   
	   1/ the daily incidence of onset of disease in Baltimore during the 1918 H1N1 influenza pandemic (see source and references),
	   
	   2/ the discrete daily distribution of the serial interval for influenza, assuming a shifted Gamma distribution with mean 2.6 days, standard deviation 1.5 days and shift 1 day (see references).
	   
     }
     \usage{Flu1918}
     \format{
		A list of two elements: 
	 
		1/ "Incidence": a vector containing 92 days of observation,
		
		2/ "SI.Distr": a vector containing a set of 12 probabilities.
		}
		\source{Frost W. and E. Sydenstricker (1919) Influenza in Maryland: preliminary statistics of certain localities. Public Health Rep.(34): 491-504.}
     \references{
       Cauchemez S. et al. (2011) Role of social networks in shaping disease transmission during a community outbreak of 2009 H1N1 pandemic influenza. Proc Natl Acad Sci U S A 108(7), 2825-2830.
	   
	   Ferguson N.M. et al. (2005) Strategies for containing an emerging influenza pandemic in Southeast Asia. Nature 437(7056), 209-214.
	   
	   Fraser C. et al. (2011) Influenza Transmission in Households During the 1918 Pandemic. Am J Epidemiol 174(5): 505-514.
	   
	   Frost W. and E. Sydenstricker (1919) Influenza in Maryland: preliminary statistics of certain localities. Public Health Rep.(34): 491-504.
	   
	   Vynnycky E. et al. (2007) Estimates of the reproduction numbers of Spanish influenza using morbidity data. Int J Epidemiol 36(4): 881-889.
	   
	   White L.F. and M. Pagano (2008) Transmissibility of the influenza virus in the 1918 pandemic. PLoS One 3(1): e1498.
     }
	 \examples{
		## load data on pandemic flu in Baltimore in 1918
		data("Flu1918")

		## estimate the reproduction number (method "NonParametricSI")
		EstimateR(Flu1918$Incidence, 
		T.Start=2:86, T.End=8:92, 
		  method="NonParametricSI", SI.Distr=Flu1918$SI.Distr, 
		  plot=TRUE, leg.pos=xy.coords(1,3))
		# the second plot produced shows, at each each day, 
		# the estimate of the reproduction number 
		# over the 7-day window finishing on that day.

	}
     \keyword{datasets}