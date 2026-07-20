# Compute the smallest index at which joint estimation should start

Unless specified by the user, `t_min` in
[`estimate_advantage()`](https://mrc-ide.github.io/EpiEstim/reference/estimate_advantage.md)
is computed as the sum of two indices:

## Usage

``` r
compute_t_min(incid, si_distr, miss_at_most)
```

## Arguments

- incid:

  a multidimensional array containing values of the incidence for each
  time step (1st dimension), location (2nd dimension) and
  pathogen/strain/variant (3rd dimension)

- si_distr:

  a matrix with two columns, each containing the probability mass
  function for the discrete serial interval for each of the two
  pathogen/strain/variants, starting with the probability mass function
  for day 0 in the first row, which should be 0. Each column in the
  matrix should sum to 1.

- miss_at_most:

  numeric. Probability mass in the tail of the SI distribution

## Value

integer

## Details

- the first day of non-zero incidence across all locations, computed
  using
  [`first_nonzero_incid()`](https://mrc-ide.github.io/EpiEstim/reference/first_nonzero_incid.md)

- the 95th percentile of the probability mass function of the SI
  distribution across all variants computed using
  [`compute_si_cutoff()`](https://mrc-ide.github.io/EpiEstim/reference/compute_si_cutoff.md)

## Author

Sangeeta Bhatia
