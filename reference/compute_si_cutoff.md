# Index before which at most a given probability mass is captured

Across a matrix of discretised probability distributions (see
[`estimate_advantage()`](https://mrc-ide.github.io/EpiEstim/reference/estimate_advantage.md))
this function returns the largest index (across all columns) such that
the cumulative probability mass before index is `1 - miss_at_most`.

## Usage

``` r
compute_si_cutoff(si_distr, miss_at_most = 0.05)
```

## Arguments

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

## Author

Sangeeta Bhatia
