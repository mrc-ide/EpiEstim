# Impute unobserved generations of infection

This function imputes incidence prior the first date of reported cases
to address early bias in `R_t` estimates. A simple linear model is
fitted on shifted, logged-incidence cases, based on an initial
observation window. Log-incidence is computed as `log(local + 0.5)` to
avoid -infinite logs. Currently, no cases are assumed to be imported.

## Usage

``` r
backimpute_I(incid, window_b)
```

## Arguments

- incid:

  the raw, reported incidence cases.

- window_b:

  length of the observation window to fit the exponential growth model
  for back-imputation

## Value

an incidence `data.frame`, combining back-imputed cases for a maximum of
100 time points (with rows indexed by a negative integer rowname) and
cases (with rows indexed by a non-negative integer)

## Examples

``` r
incid_all <- ceiling(exp(0.3 * 0:20))
incid_trunc <- tail(incid_all, 10)
x <- backimpute_I(incid = incid_trunc, window_b = 6)
idx <- as.integer(rownames(x)) > -10
x[idx, ]$local
#>  [1]   1.005090   1.517613   2.204664   3.125674   4.360311   6.015376
#>  [7]   8.234033  11.208202  15.195154  20.539768  28.000000  37.000000
#> [13]  50.000000  67.000000  91.000000 122.000000 165.000000 222.000000
#> [19] 299.000000 404.000000
incid_all
#>  [1]   1   2   2   3   4   5   7   9  12  15  21  28  37  50  67  91 122 165 222
#> [20] 299 404
```
