# Set default for Gamma priors

Set default for Gamma priors

## Usage

``` r
default_priors()
```

## Value

a list of default parameters for the priors. Values can then be manually
edited as in the examples below. Users could use functions
[`epitrix::gamma_shapescale2mucv()`](http://www.repidemicsconsortium.org/epitrix/reference/gamma_tools.md)
and
[`epitrix::gamma_mucv2shapescale()`](http://www.repidemicsconsortium.org/epitrix/reference/gamma_tools.md)
to set the shape and scale corresponding to the desired prior mean and
coefficient of variation.

## Examples

``` r
priors <- default_priors()
# change the prior for R to have a mean of 3
priors$R$shape <- 3
```
