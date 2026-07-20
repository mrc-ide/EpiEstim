# Precompute shape of posterior distribution for R

Precompute shape of posterior distribution for R

## Usage

``` r
get_shape_R_flat(incid, priors, t_min = 2L, t_max = nrow(incid))
```

## Arguments

- incid:

  a multidimensional array containing values of the (local) incidence
  for each time step (1st dimension), location (2nd dimension) and
  pathogen/strain/variant (3rd dimension)

- priors:

  a list of prior parameters (shape and scale of a gamma distribution)
  for epsilon and R; can be obtained from the function `default_priors`.
  The prior for R is assumed to be the same for all time steps and all
  locations

- t_min:

  an integer \> 1 giving the minimum time step to consider in the
  estimation. Default value is 2 (as the estimation is conditional on
  observations at time step 1 and can therefore only start at time step
  2).

- t_max:

  an integer \> `t_min` and \<= `nrow(incid)` giving the maximum time
  step to consider in the estimation. Default value is `nrow(incid)`.

## Value

a vector of the shape of the posterior distribution of R for each time
step t and each location l (stored in element
`(l-1)*(t_max - t_min + 1) + t` of the vector)

## Examples

``` r
n_v <- 2
n_loc <- 3 # 3 locations
T <- 100 # 100 time steps
priors <- default_priors()
# constant incidence 10 per day everywhere
incid <- array(10, dim = c(T, n_loc, n_v))
get_shape_R_flat(incid, priors)
#>   [1] 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04
#>  [13] 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04
#>  [25] 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04
#>  [37] 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04
#>  [49] 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04
#>  [61] 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04
#>  [73] 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04
#>  [85] 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04
#>  [97] 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04
#> [109] 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04
#> [121] 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04
#> [133] 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04
#> [145] 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04
#> [157] 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04
#> [169] 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04
#> [181] 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04
#> [193] 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04
#> [205] 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04
#> [217] 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04
#> [229] 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04
#> [241] 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04
#> [253] 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04
#> [265] 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04
#> [277] 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04
#> [289] 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04 20.04
```
