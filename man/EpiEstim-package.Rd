\encoding{UTF-8}
\name{EpiEstim-package}
\alias{EpiEstim-package}
\alias{EpiEstim}
\docType{package}
\title{The EpiEstim package}
\description{
  Quantifying transmissibility during epidemics is essential to calibrate and adjust public health responses. 
  Transmissibility can be measured by the reproduction number R, the average number of secondary cases caused by an infected individual. 
  Several methods have been proposed to estimate R over the course of an epidemic. However, they are usually difficult to implement and often specific to a particular epidemic. 
  In this package, we propose a ready-to-use tool for estimating R from the incidence time-series. 
  This tool uses novel analytical estimates of R, and also incorporates uncertainty on the serial interval distribution. 
  
  The functions implemented in this package are:
  \itemize{
         \item \code{\link{DiscrSI}}, which computes the discrete distribution of the serial interval with a given mean and standard deviation;
         \item \code{\link{OverallInfectivity}}, which calculates, at each time step, the overall infectivity due to previously infected individuals;
         \item \code{\link{EstimateR}}, which estimates, for each time step, the instantaneous reproduction number, given the incidence time series and the serial interval distribution. 
         \item \code{\link{WT}}, which estimates, for each time step, the case reproduction number, given the incidence time series and the serial interval distribution. 
       }
  
  Together with the package a few datasets are provided that contain time series of incidence and serial interval distribution for the following epidemics:
  \itemize{
         \item \code{\link{Measles1861}}: 1861 measles epidemic in Hagelloch, Germany;
		 \item \code{\link{Flu1918}}: 1918 H1N1 influenza pandemic in Baltimore;
		 \item \code{\link{Smallpox1972}}: 1972 smallpox epidemic in Kosovo;
		 \item \code{\link{SARS2003}}: 2003 severe acute respiratory syndrome (SARS) epidemic in Hong Kong;
         \item \code{\link{Flu2009}}: 2009 H1N1 influenza pandemic in a school in Pennsylvania.
       }
  
}
\references{
  Cori, A. et al. A new framework and software to estimate time-varying reproduction numbers during epidemics (AJE 2013).
}
\details{
  Package: EpiEstim
  
  Type: Package
  
  Version: 1.0-0
  
  Date: 2012-06-06
  
  License: GPL (>=2)
}
\author{ 
  Anne Cori \email{a.cori@imperial.ac.uk}
  
  with contributions from Simon Cauchemez, Neil Ferguson, Christophe Fraser and Thibaut Jombart.
}
