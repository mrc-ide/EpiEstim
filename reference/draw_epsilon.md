# Draw epsilon from marginal posterior distribution

Draw epsilon from marginal posterior distribution

## Usage

``` r
draw_epsilon(
  R,
  incid,
  lambda,
  priors,
  shape_epsilon = NULL,
  t_min = 2L,
  t_max = nrow(incid),
  seed = NULL
)
```

## Arguments

- R:

  a matrix with dimensions containing values of the instantaneous
  reproduction number for each time step (row) and location (column),
  for the reference pathogen/strain/variant

- incid:

  a multidimensional array containing values of the (local) incidence
  for each time step (1st dimension), location (2nd dimension) and
  pathogen/strain/variant (3rd dimension)

- lambda:

  a multidimensional array containing values of the overall infectivity
  for each time step (1st dimension), location (2nd dimension) and
  pathogen/strain/variant (3rd dimension). The overall infectivity for a
  given location and pathogen/strain/variant represents the sum of the
  incidence for that location and that pathogen/strain/variant at all
  previous time steps, weighted by the current infectivity of those past
  incident cases. It can be calculated from the incidence `incid` and
  the distribution of the serial interval using function
  [`compute_lambda()`](https://mrc-ide.github.io/EpiEstim/reference/compute_lambda.md)

- priors:

  a list of prior parameters (shape and scale of a gamma distribution)
  for epsilon and R; can be obtained from the function
  [`default_priors()`](https://mrc-ide.github.io/EpiEstim/reference/default_priors.md).
  The prior for R is assumed to be the same for all time steps and all
  locations

- shape_epsilon:

  a value or vector of values of the shape of the posterior distribution
  of epsilon for each of the non-reference variants, as returned by
  function
  [`get_shape_epsilon()`](https://mrc-ide.github.io/EpiEstim/reference/get_shape_epsilon.md)

- t_min:

  an integer \> 1 giving the minimum time step to consider in the
  estimation. Default value is 2 (as the estimation is conditional on
  observations at time step 1 and can therefore only start at time step
  2).

- t_max:

  an integer \> `t_min` and \<= `nrow(incid)` giving the maximum time
  step to consider in the estimation. Default value is `nrow(incid)`.

- seed:

  a numeric value used to fix the random seed

## Value

A value or vector of values for epsilon for each non reference
pathogen/strain/variant, drawn from the marginal posterior distribution

## Examples

``` r
n_loc <- 4 # 4 locations
n_v <- 3 # 3 strains
T <- 100 # 100 time steps
priors <- default_priors()
# constant incidence 10 per day everywhere
incid <- array(10, dim = c(T, n_loc, n_v))
incid <- process_I_multivariant(incid)
# arbitrary serial interval, same for both variants
w_v <- c(0, 0.2, 0.5, 0.3)
si_distr <- cbind(w_v, w_v, w_v)
lambda <- compute_lambda(incid, si_distr)
# Constant reproduction number of 1
R <- matrix(1, nrow = T, ncol = n_loc)
R[1, ] <- NA # no estimates of R on first time step
draw_epsilon(R, incid$local, lambda, priors, seed = 1)
#> [1] 1.001066 1.032584
```
