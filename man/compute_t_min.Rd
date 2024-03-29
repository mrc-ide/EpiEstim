% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gibbs_draws.R
\name{compute_t_min}
\alias{compute_t_min}
\title{Compute the smallest index at which joint estimation
should start}
\usage{
compute_t_min(incid, si_distr, miss_at_most)
}
\arguments{
\item{incid}{a multidimensional array containing values of the incidence
for each time step (1st dimension), location (2nd dimension) and
pathogen/strain/variant (3rd dimension)}

\item{si_distr}{a matrix with two columns, each containing the probability mass
function for the discrete serial interval for each of the two
pathogen/strain/variants, starting with the probability mass function
for day 0 in the first row, which should be 0. each column in the matrix
should sum to 1}

\item{miss_at_most}{numeric. probability mass in the tail of the SI distribution}
}
\value{
integer
}
\description{
Unless specified by the user, t_min in \code{estimate_advantage}
is computed as the sum of two indices:
(i) the first day of non-zero incidence across all locations,
computed using \code{first_nonzero_incid}
and (ii) the 95th percentile of the probability mass function of the
SI distribution across all variants computed using \code{compute_si_cutoff}
}
\author{
Sangeeta Bhatia
}
