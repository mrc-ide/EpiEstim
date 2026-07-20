# Convert EpiEstim distribution names to those used by coarsedatatools

coarseDataTools uses abberviated names for distributions e.g. "G" for
gamma etc To provide a smooth user experience, we convert the more
descriptive names. This function performs the user-provided names to the
abberviated ones.

## Usage

``` r
convert_distr_name_for_mcmc(distr)
```

## Arguments

- distr:

  A string with the name of the distribution as provided by the user.

## Value

A string with the converted distribution name.

## Author

Sangeeta Bhatia
