% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gibbs_draws.R
\name{compute_si_cutoff}
\alias{compute_si_cutoff}
\title{Index before which at most a given probability
mass is captured}
\usage{
compute_si_cutoff(si_distr, miss_at_most = 0.05)
}
\arguments{
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
Across a matrix of discretised probability distributions
(see \code{estimate_advantage}
this function returns the largest index
(across all columns) such that the
cumulative probability mass before index is
1 - \code{miss_at_most}.
}
\author{
Sangeeta Bhatia
}
