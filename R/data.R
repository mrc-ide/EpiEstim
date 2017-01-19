#' Data on the 1918 H1N1 influenza pandemic in Baltimore.
#'
#' This data set gives:
#' 1/ the daily incidence of onset of disease in Baltimore during the 1918 H1N1 influenza pandemic (see source and references),
#' 2/ the discrete daily distribution of the serial interval for influenza, assuming a shifted Gamma distribution with mean 2.6 days, standard deviation 1.5 days and shift 1 day (see references).
#' 
#' @format A list of two elements: 
#' \describe{
#'   \item{Incidence}{a vector containing 92 days of observation,}
#'   \item{SI.Distr}{a vector containing a set of 12 probabilities.}
#' }
#' @source Frost W. and E. Sydenstricker (1919) Influenza in Maryland: preliminary statistics of certain localities. Public Health Rep.(34): 491-504.
#' @references 
#' {
#' Cauchemez S. et al. (2011) Role of social networks in shaping disease transmission during a community outbreak of 2009 H1N1 pandemic influenza. Proc Natl Acad Sci U S A 108(7), 2825-2830.
#' 
#' Ferguson N.M. et al. (2005) Strategies for containing an emerging influenza pandemic in Southeast Asia. Nature 437(7056), 209-214.
#' 
#' Fraser C. et al. (2011) Influenza Transmission in Households During the 1918 Pandemic. Am J Epidemiol 174(5): 505-514.
#' 
#' Frost W. and E. Sydenstricker (1919) Influenza in Maryland: preliminary statistics of certain localities. Public Health Rep.(34): 491-504.
#' 
#' Vynnycky E. et al. (2007) Estimates of the reproduction numbers of Spanish influenza using morbidity data. Int J Epidemiol 36(4): 881-889.
#' 
#' White L.F. and M. Pagano (2008) Transmissibility of the influenza virus in the 1918 pandemic. PLoS One 3(1): e1498.
#' }
#' @examples 
#' ## load data on pandemic flu in Baltimore in 1918
#' data("Flu1918")
#' 
#' ## estimate the reproduction number (method "NonParametricSI")
#' EstimateR(Flu1918$Incidence, 
#'           T.Start=2:86, T.End=8:92, 
#'           method="NonParametricSI", SI.Distr=Flu1918$SI.Distr, 
#'           plot=TRUE, leg.pos=xy.coords(1,3))
#' # the second plot produced shows, at each each day, 
#' # the estimate of the reproduction number 
#' # over the 7-day window finishing on that day.
"Flu1918"

##################################################################################################

#' Data on the 2009 H1N1 influenza pandemic in a school in Pennsylvania.
#'
#' This data set gives:
#' 1/ the daily incidence of onset of acute respiratory illness 
#' (ARI, defined as at least two symptoms among fever, cough, sore throat, and runny nose)
#' amongst children in a school in Pennsylvania during the 2009 H1N1 influenza pandemic (see source and references),
#' 2/ the discrete daily distribution of the serial interval for influenza, assuming a shifted Gamma distribution with mean 2.6 days, standard deviation 1.5 days and shift 1 day (see references).
#' @format A list of two elements: 
#' \describe{
#'   \item{Incidence}{a vector containing 32 days of observation,}
#'   \item{SI.Distr}{a vector containing a set of 12 probabilities.}
#' }
#' @source Cauchemez S. et al. (2011) Role of social networks in shaping disease transmission during a community outbreak of 2009 H1N1 pandemic influenza. Proc Natl Acad Sci U S A 108(7), 2825-2830.
#' @references 
#' {
#' Cauchemez S. et al. (2011) Role of social networks in shaping disease transmission during a community outbreak of 2009 H1N1 pandemic influenza. Proc Natl Acad Sci U S A 108(7), 2825-2830.
#' 
#' Ferguson N.M. et al. (2005) Strategies for containing an emerging influenza pandemic in Southeast Asia. Nature 437(7056), 209-214.
#' }
#' @examples 
#' ## load data on pandemic flu in a school in 2009
#' data("Flu2009")
#' 
#' ## estimate the reproduction number (method "NonParametricSI")
#' EstimateR(Flu2009$Incidence, 
#'           T.Start=2:26, T.End=8:32, 
#'           method="NonParametricSI", SI.Distr=Flu2009$SI.Distr, 
#'           plot=TRUE, leg.pos=xy.coords(1,3))
#' # the second plot produced shows, at each each day, 
#' # the estimate of the reproduction number 
#' # over the 7-day window finishing on that day.
"Flu2009"


##################################################################################################

#' Data on the 1861 measles epidemic in Hagelloch, Germany.
#'
#' This data set gives:
#' 1/ the daily incidence of onset of symptoms in Hallegoch (Germany) during the 1861 measles epidemic (see source and references),
#' 2/ the discrete daily distribution of the serial interval for measles, assuming a shifted Gamma distribution with mean 14.9 days, standard deviation 3.9 days and shift 1 day (see references).
#' @format A list of two elements: 
#' \describe{
#'   \item{Incidence}{a vector containing 48 days of observation,}
#'   \item{SI.Distr}{a vector containing a set of 38 probabilities.}
#' }
#' @source Groendyke C. et al. (2011) Bayesian Inference for Contact Networks Given Epidemic Data. Scandinavian Journal of Statistics 38(3): 600-616.
#' @references Groendyke C. et al. (2011) Bayesian Inference for Contact Networks Given Epidemic Data. Scandinavian Journal of Statistics 38(3): 600-616.
#' @examples 
#' ## load data on measles in Hallegoch in 1861
#' data("Measles1861")
#' 
#' ## estimate the reproduction number (method "NonParametricSI")
#' EstimateR(Measles1861$Incidence, 
#'           T.Start=17:42, T.End=23:48, 
#'           method="NonParametricSI", SI.Distr=Measles1861$SI.Distr, 
#'           plot=TRUE, leg.pos=xy.coords(1,7))
#' # the second plot produced shows, at each each day, 
#' # the estimate of the reproduction number 
#' # over the 7-day window finishing on that day.
"Measles1861"

##################################################################################################

#' Data on the 2003 SARS epidemic in Hong Kong.
#'
#' This data set gives:
#' 1/ the daily incidence of onset of symptoms in Hong Kong during the 2003 severe acute respiratory syndrome (SARS) epidemic (see source and references),
#' 2/ the discrete daily distribution of the serial interval for SARS, assuming a shifted Gamma distribution with mean 8.4 days, standard deviation 3.8 days and shift 1 day (see references).
#' @format A list of two elements: 
#' \describe{
#'   \item{Incidence}{a vector containing 107 days of observation,}
#'   \item{SI.Distr}{a vector containing a set of 25 probabilities.}
#' }
#' @source Cori A. et al. (2009) Temporal variability and social heterogeneity in disease transmission: the case of SARS in Hong Kong. PLoS Comput Biol 5(8): e1000471.
#' @references 
#' {
#' Cori A. et al. (2009) Temporal variability and social heterogeneity in disease transmission: the case of SARS in Hong Kong. PLoS Comput Biol 5(8): e1000471.
#' 
#' Lipsitch M. et al. (2003) Transmission dynamics and control of severe acute respiratory syndrome. Science 300(5627): 1966-1970.
#' }
#' @examples 
#' ## load data on SARS in Hong Kong in 2003
#' data("SARS2003")
#' 
#' ## estimate the reproduction number (method "NonParametricSI")
#' EstimateR(SARS2003$Incidence, 
#'           T.Start=14:101, T.End=20:107, 
#'           method="NonParametricSI", SI.Distr=SARS2003$SI.Distr, 
#'           plot=TRUE, leg.pos=xy.coords(1,7))
#' # the second plot produced shows, at each each day, 
#' # the estimate of the reproduction number 
#' # over the 7-day window finishing on that day.
"SARS2003"

##################################################################################################

#' Data on the 1972 smallpox epidemic in Kosovo
#'
#' This data set gives:
#' 1/ the daily incidence of onset of symptoms in Kosovo during the 1972 smallpox epidemic (see source and references),
#' 2/ the discrete daily distribution of the serial interval for smallpox, assuming a shifted Gamma distribution with mean 22.4 days, standard deviation 6.1 days and shift 1 day (see references).
#' @format A list of two elements: 
#' \describe{
#'   \item{Incidence}{a vector containing 57 days of observation,}
#'   \item{SI.Distr}{a vector containing a set of 46 probabilities.}
#' }
#' @source Fenner F. et al. (1988) Smallpox and its Eradication. Geneva, World Health Organization.
#' @references 
#' {
#' Fenner F. et al. (1988) Smallpox and its Eradication. Geneva, World Health Organization.
#' 
#' Gani R. and S. Leach (2001) Transmission potential of smallpox in contemporary populations. Nature 414(6865): 748-751.
#'
#' Riley S. and N. M. Ferguson (2006) Smallpox transmission and control: spatial dynamics in Great Britain. Proc Natl Acad Sci U S A 103(33): 12637-12642.
#' }
#' @examples 
#' ## load data on smallpox in Kosovo in 1972
#' data("Smallpox1972")
#' 
#' ## estimate the reproduction number (method "NonParametricSI")
#' EstimateR(Smallpox1972$Incidence, 
#'           T.Start=27:51, T.End=33:57, 
#'           method="NonParametricSI", SI.Distr=Smallpox1972$SI.Distr, 
#'           plot=TRUE, leg.pos=xy.coords(1,15))
#' # the second plot produced shows, at each each day, 
#' # the estimate of the reproduction number 
#' # over the 7-day window finishing on that day.
"Smallpox1972"

##################################################################################################

#' Mock data on a rotavirus epidemic.
#'
#' This data set gives:
#' 1/ the daily incidence of onset of symptoms in a mock outbreak of rotavirus ,
#' 2/ mock observations of symptom onset dates for 19 pairs of infector/infected individuals.
#' @format A list of two elements: 
#' \describe{
#'   \item{Incidence}{a vector containing 53 days of observation,}
#'   \item{SI.Distr}{a dataframe containing a set of 19 observations; each observation corresponds to a pair of infector/infected individuals. EL and ER columns contain the lower an upper bounds of the dates of symptoms onset in the infectors. SL and SR columns contain the lower an upper bounds of the dates of symptoms onset in the infected indiviuals. The type column corresponds to XXX TO BE COMPLETED XXX}
#' }
#' @source XXX TO BE COMPLETED XXX
#' @references 
#' {
#' XXX TO BE COMPLETED XXX
#' }
#' @examples 
#' \dontrun{
#' ## Note the following example uses an MCMC routine 
#' ## to estimate the serial interval distribution from data, 
#' ## so may take a few minutes to run
#' 
#' ## load data 
#' data("MockRotavirus")
#' 
#' ## estimate the reproduction number (method "SIFromData")
#' EstimateR(MockRotavirus$Incidence, 
#'           T.Start=2:47, T.End=8:53, 
#'           method="SIFromData", 
#'           SI.Data=MockRotavirus$SI.Data, 
#'           SI.parametricDistr = "G", 
#'           MCMC.control = list(burnin = 3000), 
#'           n1 = 1000, n2 = 50,
#'           plot=TRUE, leg.pos=xy.coords(1,3))
#' # the second plot produced shows, at each each day, 
#' # the estimate of the reproduction number 
#' # over the 7-day window finishing on that day.
#' }
"MockRotavirus"



