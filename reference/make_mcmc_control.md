# Create list of MCMC control parameters

Creates a list of MCMC control parameters to be used in
`config$mcmc_control`, where `config` is an argument of the
[`estimate_R()`](https://mrc-ide.github.io/EpiEstim/reference/estimate_R.md)
function. This is used to configure the MCMC chain used to estimate the
serial interval within
[`estimate_R()`](https://mrc-ide.github.io/EpiEstim/reference/estimate_R.md)
(with method "si_from_data").

## Usage

``` r
make_mcmc_control(
  burnin = 3000,
  thin = 10,
  seed = as.integer(Sys.time()),
  init_pars = NULL
)
```

## Arguments

- burnin:

  A positive integer giving the burnin used in the MCMC when estimating
  the serial interval distribution.

- thin:

  A positive integer corresponding to thinning parameter; the MCMC will
  be run for `burnin + n1*thin` iterations; 1 in `thin` iterations will
  be recorded, after the burnin phase, so the posterior sample size is
  `n1`.

- seed:

  An integer used as the seed for the random number generator at the
  start of the MCMC estimation; useful to get reproducible results.

- init_pars:

  vector of size 2 corresponding to the initial values of parameters to
  use for the SI distribution. This is the shape and scale for all but
  the lognormal distribution, for which it is the meanlog and sdlog.

## Value

An object of class `estimate_R_mcmc_control` with components `burnin`,
`thin`, `seed`, `init_pars`. This can be used as an argument of function
[`make_config()`](https://mrc-ide.github.io/EpiEstim/reference/make_config.md).

## Details

The argument `si_data`, should be a dataframe with 5 columns:

- `EL`: the lower bound of the symptom onset date of the infector (given
  as an integer)

- `ER`: the upper bound of the symptom onset date of the infector (given
  as an integer). Should be such that `ER >= EL`

- `SL`: the lower bound of the symptom onset date of the infected
  individual (given as an integer)

- `SR`: the upper bound of the symptom onset date of the infected
  individual (given as an integer). Should be such that `SR >= SL`

- `type` (optional): can have entries 0, 1, or 2, corresponding to
  doubly interval-censored, single interval-censored or exact
  observations, respectively, see Reich et al. Statist. Med. 2009. If
  not specified, this will be automatically computed from the dates

Assuming a given parametric distribution for the serial interval
distribution (specified in `si_parametric_distr`), the posterior
distribution of the serial interval is estimated directly from these
data using MCMC methods implemented in the package

## Examples
