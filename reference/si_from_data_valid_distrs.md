# Distribution names valid when using MCMC to estimate SI from data

When using si_from_data method, the package will use
[`dic.fit.mcmc`](http://nickreich.github.io/coarseDataTools/reference/dic.fit.mcmc.md)
to fit the serial interval distribution. This method supports only a
limited set of distributions. This function returns the valid
distribution names as used by EpiEstim. The names are internally
converted to those used by coarsedatatools by
[`convert_distr_name_for_mcmc`](https://mrc-ide.github.io/EpiEstim/reference/convert_distr_name_for_mcmc.md)
function.

## Usage

``` r
si_from_data_valid_distrs(dist)
```

## Arguments

- dist:

  The parametric distribution used when estimating the serial interval.
  Should be one of "gamma", "weibull", "lognormal", "gamma_offset_1",
  "weibull_offset_1", or "lognormal_offset_1". Note the different naming
  convention compared to
  [`dic.fit.mcmc`](http://nickreich.github.io/coarseDataTools/reference/dic.fit.mcmc.md).
  The distribution may also be specified using the abbreviated forms
  "G", "W", "L", "off1G", "off1W", and "off1L" as used in
  [`dic.fit.mcmc`](http://nickreich.github.io/coarseDataTools/reference/dic.fit.mcmc.md).
  However, we recommend using the full names to avoid confusion, and a
  warning will be issued if the abbreviated forms are used. If not
  present, computed automatically from `x`.

## Value

A two element list - the first element is a flag `is_dist_valid`
indicating whether the passed distribution is valid. Te second element
`all_valid_distrs` is a character vector with the valid distribution
names.

## Author

Sangeeta Bhatia
