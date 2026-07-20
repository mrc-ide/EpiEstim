# Multivariant EpiEstim

The `EpiEstim` package has been extended to allow users to estimate the
effective transmission advantage of new pathogens/variants/strains
compared to a reference pathogen/variant/strain in real-time. This
vignette will take you through the key stages of estimating this
advantage and also provide an example of how this method can be applied
in a typical outbreak analysis workflow. Please also see the ‘FAQs’
section.

### Input data

The **incidence data** needs to be supplied in the form of a
multidimensional array, containing the incidence for each day of the
outbreak (row), location (column) and variant/strain under investigation
(third dimension), e.g.:

[TABLE]

The **serial interval distributions** must be supplied as a matrix with
the number of columns matching the number of variants being considered.
Each column should contain the probability mass function (PMF) for the
discrete serial interval of each variant, starting with the PMF for day
zero in the first row (which should be 0) and each column should sum to
1, e.g.:

|        | Variant_1 | Variant_2 |
|:-------|:----------|:----------|
| \[1,\] | 0         | 0         |
| \[2,\] | 0.1       | 0.2       |
| \[3,\] | 0.3       | 0.4       |
| \[4,\] | 0.5       | 0.6       |
| …      | …         | …         |

\##Estimate advantage

To estimate the effective transmission advantage of new variants
compared to a reference variant we use the
[`estimate_advantage()`](https://mrc-ide.github.io/EpiEstim/reference/estimate_advantage.md)
function in `EpiEstim`.

\###[`estimate_advantage()`](https://mrc-ide.github.io/EpiEstim/reference/estimate_advantage.md)

``` r

estimate_advantage(incid = incidence_object, 
                   si_distr = si_matrix, 
                   priors = default_priors(),
                   mcmc_control = default_mcmc_controls(),
                   t_min = NULL, t_max = nrow(incid),
                   incid_imported = NULL,
                   precompute = TRUE,
                   reorder_incid = TRUE)
```

As described above, the `incid=` and `si_distr=` arguments should be
supplied with a multidimensional array of incidence and a matrix
containing the serial interval for each variant respectively. In
addition, the user can supply:

- `priors=` with a list of the prior parameters (shape and scale) for
  the reproduction number and the transmission advantage. The default
  priors can be obtained using the
  [`default_priors()`](https://mrc-ide.github.io/EpiEstim/reference/default_priors.md)
  function.

- `mcmc_control=` with a list of the default properties of the MCMC,
  i.e. the number of iterations, burn-in and the thinning of the MCMC
  chains. Burn-in is the number of iterations to be discarded as
  “warm-up”, for instance, if burnin=10 the output is only recorded
  after 10 iterations have run. Thinning determines how much the MCMC
  chains should be thinned out, if thin = 10 then 1 in every 10
  iterations will be kept. The default parameters can be obtained using
  [`default_mcmc_controls()`](https://mrc-ide.github.io/EpiEstim/reference/default_mcmc_controls.md).

- `t_min=` and `t_max=` with integers that give the minimum and maximum
  time step over which the transmission advantage will be estimated.
  Note that cases before `t_min` will still be accounted for as
  potential infectors in the likelihood. Indeed, cases before `t_min`
  can contribute to the overall infectiousness from t_min to t_max
  through the serial interval distribution. `t_min` must be \>1 and if
  `t_min = NULL` then this is automatically calculated as the maximum of
  the 95th percentile of the SI distribution. `t_max` must be \> `t_min`
  and \<= the number of rows of the incidence data. The default value of
  `t_max` is `nrow(incid)`.

Other optional parameters include:

- `incid_imported=` where a multidimensional array can be supplied (row
  = time, column = location, third dimension = variant) containing the
  incidence of imported cases, such that `incid - incid_imported` would
  be the incidence of locally infected cases. If `incid_imported = NULL`
  it will be assumed that (other than in the first time step) all cases
  arose from local transmission.

- `precompute=` which can be `TRUE` or `FALSE`, but defaults to `TRUE`.
  This determines whether the shape of the posterior distributions for R
  and the transmission advantage (for the non-reference variant(s)) for
  each time step and location should be precalculated. This makes
  [`estimate_advantage()`](https://mrc-ide.github.io/EpiEstim/reference/estimate_advantage.md)
  faster, only use `precompute=FALSE` for debugging.

- `reorder_incid=` which can be `TRUE` or `FALSE`, but defaults to
  `TRUE`. If `TRUE` then during the estimation of the transmission
  advantage the incidence array will be temporarily re-ordered so that
  all variants are compared to the most transmissible variant
  (temporarily considered the reference) regardless of the order in
  which the user supplied the data. Note that the results will always be
  returned in the order supplied by the user. We recommend that this is
  set to `TRUE`.

## Examples

### SARS-CoV-2 variants

This example will take you through a workflow using incidence data for
two variants of SARS-CoV-2 (wildtype and alpha) across 7 NHS regions in
England. (For detailed description of the data see Bhatia et al. 2021)

``` r

incid <- readRDS("./data_mv_vignette/incid.rds")
```

**Incidence**

Let us say we have incidence data in the format below, where each row
corresponds to a time step in the outbreak, each column corresponds to a
region in England, and each third dimension corresponds to a variant.

``` r

head(incid)
#> , , wild
#> 
#>      East of England London Midlands North East and Yorkshire North West
#> [1,]               7     19       28                       21         37
#> [2,]               8     20       44                       30         48
#> [3,]              10     19       43                       30         48
#> [4,]              10     23       52                       30         58
#> [5,]               8     23       46                       36         53
#> [6,]               6     19       32                       14         37
#>      South East South West
#> [1,]          7          8
#> [2,]         10          7
#> [3,]         11          6
#> [4,]         10          6
#> [5,]         12          5
#> [6,]         13          2
#> 
#> , , alpha
#> 
#>      East of England London Midlands North East and Yorkshire North West
#> [1,]               0      0        1                        0          0
#> [2,]               0      1        1                        0          1
#> [3,]               0      0        1                        0          0
#> [4,]               1      0        1                        0          0
#> [5,]               0      0        0                        0          0
#> [6,]               0      0        0                        0          1
#>      South East South West
#> [1,]          0          0
#> [2,]          0          0
#> [3,]          1          1
#> [4,]          0          0
#> [5,]          0          0
#> [6,]          0          0
```

We also assume that the SI is the same for both variants, with a mean of
5.4 days and a standard deviation of 1.5 days (Rai, Shukla, and Dwivedi
2021) and produce a matrix of SI’s with the number of columns equal to
the number of variants:

``` r

mean_SI <- 5.4
sd_SI <- 1.5
SI <- EpiEstim::discr_si(seq(0, 20), mean_SI, sd_SI)
si_matrix <- cbind(SI,SI)
```

``` r

head(si_matrix)
#>                SI           SI
#> [1,] 0.000000e+00 0.000000e+00
#> [2,] 4.674689e-05 4.674689e-05
#> [3,] 8.034951e-03 8.034951e-03
#> [4,] 7.964056e-02 7.964056e-02
#> [5,] 2.125381e-01 2.125381e-01
#> [6,] 2.683362e-01 2.683362e-01
```

**Estimate transmission advantage**

Now that we have our incidence and SI distributions we can supply these
to the
[`estimate_advantage()`](https://mrc-ide.github.io/EpiEstim/reference/estimate_advantage.md)
function. We use the default priors and MCMC controls.

``` r

output <- EpiEstim::estimate_advantage(incid = incid,
                             si_distr = si_matrix,
                             priors = default_priors(),
                             mcmc_control = default_mcmc_controls())
```

The output is a list containing 4 elements:

``` r

output$epsilon
```

- epsilon - the estimated effective transmission advantage. This is a
  matrix containing the MCMC chain (thinned and after burnin) for the
  relative transmissibility of the “new” variant (alpha) compared to the
  reference variant (wildtype) across all locations. Each row in the
  matrix is a “new” variant (in this case only 1) and each column an
  iteration of the MCMC. If epsilon \>1 then the corresponding variant
  is estimated to have a transmission advantage over the reference
  variant.

``` r

output$R
```

- R - the reproduction number estimate for the reference variant. This
  is an array containing the MCMC chain (thinned and after burnin) for
  the reproduction number of the reference variant. The first dimension
  of the array is time, the second is the location, and the third is the
  iteration of the MCMC.

In this example, we did not specify `t_min=` or `t_max=`. You will
notice that the R estimates do not start until day 25, this is
because: 1) One of the regions doesn’t have cases of the alpha variant
until day 16, 2) The 95th percentile of the SI is 9 days.

``` r

output$convergence
```

- convergence - result of the Gelman-Rubin convergence diagnostic. This
  is either ‘TRUE’ or ‘FALSE’ and tells you whether the corresponding
  epsilon MCMC chain for each “new” variant has converged within the
  number of iterations specified.

``` r

output$diag
```

- diag - nested list of the point estimate and upper confidence limits
  of the Gelman-Rubin convergence diagnostics (as implemented in the R
  package `coda`). The length of diag is equal to the number of rows in
  epsilon (i.e. the number of non-reference variants, in this case only
  1). Each element of diag contains a named list of the point estimate
  and upper confidence limits.

**Summarise results**

To summarise the estimated transmission advantage one can use the
posterior median and 95% CrI as follows.

``` r

quantile(output$epsilon[1,], c(0.5,0.025,0.975))
#>      50%     2.5%    97.5% 
#> 1.491769 1.470332 1.514741
```

## FAQs

- **Why should the estimates of the reproduction number produced by
  [`estimate_advantage()`](https://mrc-ide.github.io/EpiEstim/reference/estimate_advantage.md)
  be interpreted with caution?**

The inference framework jointly estimates the instantaneous reproduction
number of the reference variant and the effective transmission advantage
of new variants. R estimates here should be interpreted with caution
because they represent estimates from the joint distribution of the
reproduction number and transmission advantage, therefore depending on
the incidence of the reference *as well as* the new variant(s).
Moreover, temporal variation in estimates of the instantaneous
reproduction number can arise from a number of factors, such as the
implementation of control measures, changes in population behaviour, or
variability in the reporting of cases over time. If only interested in
estimating the reference R_(t), we recommend using
[`estimate_R()`](https://mrc-ide.github.io/EpiEstim/reference/estimate_R.md).

- **What do some of the common warning/error messages mean?**

&nbsp;

    "The Gelman-Rubin algorithm suggests the MCMC may not have converged within the number of iterations specified."

This means that the MCMC chains for the estimation have not converged.
To avoid this you may need to run
[`estimate_advantage()`](https://mrc-ide.github.io/EpiEstim/reference/estimate_advantage.md)
with more MCMC iterations. E.g. by altering `mcmc_control()`.

    "Input SI distributions should sum to 1. Normalising now"

The SI distributions supplied have been re-normalised so that they sum
to 1.
