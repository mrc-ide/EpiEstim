# Check MCMC chain convergence using the Gelman-Rubin algorithm

This function splits an MCMC chain in two halves and uses the
Gelman-Rubin algorithm to assess convergence of the chain by comparing
its two halves.

## Usage

``` r
check_cdt_samples_convergence(cdt_samples)
```

## Arguments

- cdt_samples:

  the `@sample` slot of a `cd.fit.mcmc` S4 object (see package
  `coarseDataTools`)

## Value

TRUE if the Gelman Rubin test for convergence was successful, FALSE
otherwise

## See also

[`estimate_R()`](https://mrc-ide.github.io/EpiEstim/reference/estimate_R.md)

## Author

Anne Cori

## Examples

``` r
if (FALSE) { # \dontrun{
## Note the following examples use an MCMC routine
## to estimate the serial interval distribution from data,
## so they may take a few minutes to run

## load data on rotavirus
data("MockRotavirus")

## estimate the serial interval from data
SI_fit <- coarseDataTools::dic.fit.mcmc(dat = MockRotavirus$si_data,
                     dist = "G",
                     init_pars = init_mcmc_params(MockRotavirus$si_data, "G"),
                     burnin = 1000,
                     n.samples = 5000)

## use check_cdt_samples_convergence to check convergence
converg_diag <- check_cdt_samples_convergence(SI_fit@samples)
converg_diag

} # }
```
