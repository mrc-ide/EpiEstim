# Compute the overall infectivity

Compute the overall infectivity

## Usage

``` r
compute_lambda(incid, si_distr)
```

## Arguments

- incid:

  a list (as obtained from function
  [`process_I_multivariant()`](https://mrc-ide.github.io/EpiEstim/reference/process_I_multivariant.md))
  of two multidimensional arrays (`local` and `imported`) containing
  values of the incidence for each time step (1st dimension), location
  (2nd dimension) and pathogen/strain/variant (3rd dimension)

- si_distr:

  a matrix where each column contains the probability mass function for
  the discrete serial interval for each of the pathogen/strain/variants,
  starting with the probability mass function for day 0 in the first
  row, which should be 0. Each column in the matrix should sum to 1

## Value

a multidimensional array containing values of the overall infectivity
for each time step (1st dimension), location (2nd dimension) and
pathogen/strain/variant (3rd dimension). The overall infectivity for a
given location and pathogen/strain/variant represents the sum of the
incidence for that location and that pathogen/strain/variant at all
previous time steps, weighted by the current infectivity of those past
incident cases. Pre-calculating the overall infectivity makes the
algorithm much faster.

## Examples

``` r
n_v <- 2
n_loc <- 3 # 3 locations
T <- 100 # 100 time steps
priors <- default_priors()
# constant incidence 10 per day everywhere
incid <- array(10, dim = c(T, n_loc, n_v))
incid <- process_I_multivariant(incid)
# arbitrary serial interval, same for both variants
w_v <- c(0, 0.2, 0.5, 0.3)
si_distr <- cbind(w_v, w_v)
lambda <- compute_lambda(incid, si_distr)
```
