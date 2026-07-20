# Set default for MCMC control

Set default for MCMC control

## Usage

``` r
default_mcmc_controls()
```

## Value

a list of default MCMC control parameters, containing:

- `n_iter`: the number of iterations for the MCMC to perform

- `burnin`: the burnin to use; MCMC iterations will only be recorded
  after the burnin

- `thin`: MCMC iterations will only be recorded after the burnin and
  every `thin` iteration

Values can then be manually edited as in the examples below.

## Examples

``` r
mcmc_control <- default_mcmc_controls()
# change to run for 10 times longer
mcmc_control$n_iter <- mcmc_control$n_iter * 10
```
