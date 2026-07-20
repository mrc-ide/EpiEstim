# EpiEstim for aggregated incidence data

The `EpiEstim` package has been extended to allow users to estimate the
time varying reproduction number (R_(t)) from temporally aggregated
incidence data (Nash et al.). This approach reconstructs daily incidence
from data supplied at any timescale. This vignette will take you through
the different ways aggregated incidence data can be supplied, the
additional parameters needed, and an example using aggregated UK
COVID-19 data. Please also see the FAQs section.

### Input data

**Incidence data**

Many diseases are not reported on a daily basis. EpiEstim can now use
incidence data that has been aggregated in multiple ways, e.g.:

- Constant aggregations, such as weekly data
- Repeating patterns of aggregations, such as regular reporting 3x per
  week
- Variable aggregations over the course of an outbreak

Daily incidence data is reconstructed from aggregated data using an
expectation-maximisation (EM) algorithm. There are three stages of the
EM algorithm:

- **Initialisation.** The EM algorithm is initialised with a naive
  disaggregation of the incidence data. For example, if there were 70
  cases over the course of a week, this would be naively split into 10
  cases per day.
- **Expectation.** The reproduction number is estimated for each
  aggregation window, except for the first aggregation window (as there
  is no past incidence data). This means that the earliest the incidence
  reconstruction can start is at least the first day of the second
  aggregation window. Additionally, if the disaggregated incidence in
  subsequent aggregation windows is too low to estimate the reproduction
  number, this will mean that the reconstruction will not start until
  case numbers are sufficiently high.
- **Maximisation.** The reproduction number estimates are then
  translated into growth rates for each aggregation window (Wallinga &
  Lipsitch, 2007) and used to reconstruct daily incidence data assuming
  exponential growth. The daily incidence is adjusted by a constant to
  ensure that if the daily incidence were to be re-aggregated, it would
  still sum to the original aggregated totals. The expectation and
  maximisation steps repeat iteratively until convergence.

The daily incidence that is reconstructed after the final iteration of
the EM algorithm is then used to estimate R_(t) using the same process
as the original
[`estimate_R()`](https://mrc-ide.github.io/EpiEstim/reference/estimate_R.md)
function, with sliding weekly time windows used as the default.

Note that we do not support the distinction between local and imported
cases when using temporally aggregated incidence as we assume that this
level of data would not be available.

**Aggregation windows**

Aggregation windows can be specified using the parameter `dt`, and can
be provided in one of three ways:

- A single integer - for constant aggregation windows, e.g. `dt = 7L`
  for weekly data
- A repeating vector of integers - for repeating aggregation patterns
  e.g. reporting 3x per week on the same day of the week could be
  `dt = c(2L,2L,3L)`
- A vector of aggregations matching the length of incidence data
  supplied

**Serial interval distribution**

The serial interval can be provided on a daily timescale (as usual),
either as the mean and standard deviation (parametric distribution) or
the full distribution (non-parametric distribution). See
‘full_EpiEstim_vignette’ for more details.

### Estimate R_(t) from temporally aggregated incidence data

To estimate R_(t) from temporally aggregated incidence data, we simply
use the
[`estimate_R()`](https://mrc-ide.github.io/EpiEstim/reference/estimate_R.md)
function with two additional parameters required, `dt` and `dt_out`, and
some optional parameters, `recon_opt`, `iter`, `tol`, and `grid`.

#### `estimate_R()` for aggregated data

``` r

estimate_R(incid = aggregated_incidence,
           dt = 7L,
           dt_out = 7L,
           recon_opt = "naive",
           iter = 10L,
           tol = 1e-6,
           grid = list(precision = 0.001, min = -1, max = 1),
           config = config,
           method = method)
```

- As described above, `dt` can be supplied as a single integer, a vector
  of repeating integers, or a full vector of integers matching the
  length of the incidence data.

- `dt_out` is the length of the sliding windows used to estimate R_(t)
  from the reconstructed daily incidence data, this is `dt = 7L` (weekly
  sliding windows) by default. We recommend that `dt_out` is at least
  equal to the length of the longest aggregation window (`dt`) in the
  data.

- `recon_opt` can be one of two options: `"naive"` or `"match"`. This
  specifies how to handle the initial incidence data that cannot be
  reconstructed by the EM algorithm (e.g. the incidence data for the
  aggregation window that precedes the first aggregation window that R
  can be estimated for). If `"naive"` is chosen, the naive
  disaggregation of the incidence data will be kept. If `"match"` is
  chosen, the incidence in the preceding aggregation window will be
  reconstructed by assuming that the growth rate matches that of the
  first estimation window. This is `"naive"` by default.

There are three other optional parameters that can be modified, however,
we recommend that the default values are used:

- `iter` is the number of iterations of the EM algorithm used to
  reconstruct the daily incidence data. This is `iter = 10L` by default.

- `tol` is the tolerance value used for the convergence check. The
  tolerance is how much the final iteration of the reconstructed daily
  incidence is allowed to differ from the reconstructed incidence
  produced in the previous iteration without returning a warning. This
  is `tol = 1e-6` by default.

- `grid` is a list of “precision”, “min”, and “max” values to define a
  grid of growth rate parameters used inside the EM algorithm. The grid
  is used to convert reproduction number estimates for each aggregation
  of incidence data into growth rates, which are then used to
  reconstruct the daily incidence data assuming exponential growth. The
  grid will auto-adjust if it is not large enough, so we recommend using
  the default values.

The SI distibution can be specified as normal using the `method` and
`config` parameters (see the full_EpiEstim_vignette for more details).

## Examples

### Estimate R_(t) from weekly COVID-19 data

This example will take you through a workflow using weekly incidence
data for UK COVID-19 cases. (For detailed description of the data see
Nash et al.)

``` r

incid <- readRDS("./aggregated_data/UK_covid_cases.rds")
```

**Incidence**

Let us say we have a vector of weekly incidence data for COVID-19 cases.

``` r

incid
#>  [1]      21     241    1503    4714   14294   27408   33124   30126   33288
#> [10]   31944   25581   20435   17320   12195    9259    7102    6948    5596
#> [19]    4376    4258    4287    4704    5507    5964    7408    7500    8244
#> [28]   13823   22060   23924   41711   66725  106925  116472  146015  150832
#> [37]  159310  172001  146202  107282  102508  128194  202680  260353  352685
#> [46]  397790  316549  250624  175157  133220   94886   79373   63185   41334
#> [55]   38741   37560   37213   27387   19128   16934   16182   14287   13422
#> [64]   14743   15882   20430   30671   46961   62642   97320  163612  214681
#> [73]  294984  279099  187392  186748  198318  225124  239123  249063  251378
#> [82]  199736  235620  234391  244202  292109  326688  280646  255388  253735
#> [91]  287826  295546  322535  345941  525261  800717 1034989
```

We need to specify how the data is aggregated, which in this case, is by
constant weekly aggregation windows. We do this by supplying `dt` with a
single integer (7L).

``` r

dt <- 7L
```

We can take an estimate from the literature to specify a parametric SI
with a mean of 6.3 days and a standard deviation of 4.2 days (Bi et al
2020).

``` r

mean_si <- 6.3
std_si <- 4.2
method <- "parametric_si"
config <- make_config(list(mean_si = mean_si,
                           std_si = std_si))
```

**Estimate R_(t)**

Now that we have our aggregated incidence, our aggregation time window,
and SI distribution, we can supply these to the
[`estimate_R()`](https://mrc-ide.github.io/EpiEstim/reference/estimate_R.md)
function. We do not need to specify `dt_out`, `iter`, `tol`, or `grid`
because we are going to use the default values.

``` r

output <- EpiEstim::estimate_R(incid = incid,
                               dt = dt,
                               recon_opt = "match",
                               method = method,
                               config = config)
```

The output consists of multiple elements, including the reconstructed
daily incidence data:

``` r

head(output$I)
#> [1] 0.817031 1.164069 1.658514 2.362976 3.366662 4.796669
```

And the R_(t) estimates (in this case, using the default weekly sliding
windows):

``` r

head(output$R)
#>   t_start t_end  Mean(R)    Std(R) Quantile.0.025(R) Quantile.0.05(R)
#> 1       8    14 4.852244 0.3119140          4.260110         4.350840
#> 2       9    15 4.491777 0.2516604          4.012019         4.086002
#> 3      10    16 4.169784 0.2042461          3.779027         3.839612
#> 4      11    17 3.934157 0.1688907          3.610057         3.660539
#> 5      12    18 3.766547 0.1420022          3.493333         3.536057
#> 6      13    19 3.645352 0.1209050          3.412203         3.448787
#>   Quantile.0.25(R) Median(R) Quantile.0.75(R) Quantile.0.95(R)
#> 1         4.638377  4.845562         5.058827         5.376440
#> 2         4.319571  4.487078         4.658862         4.913581
#> 3         4.030265  4.166450         4.305669         4.511330
#> 4         3.818963  3.931740         4.046717         4.216018
#> 5         3.669820  3.764763         3.861329         4.003124
#> 6         3.563091  3.644016         3.726157         3.846477
#>   Quantile.0.975(R)
#> 1          5.482348
#> 2          4.998238
#> 3          4.579489
#> 4          4.271989
#> 5          4.049902
#> 6          3.886097
```

In this example, you will notice that R_(t) estimation does not start
until day 8. This is because the daily incidence data cannot be
reconstructed, and R_(t) estimation cannot start, until the first day of
the second aggregation window. The start of R_(t) estimation may also be
delayed if incidence is too low, but this was not the case here.

**Plot results**

As normal, simply plot the full or partial output.

``` r

plot(output) # full output

plot(output, "incid") # Reconstructed daily incidence only
plot(output, "R") # Rt estimates only
plot(output, "SI") # SI estimates only
```

**Convergence**

Convergence is checked automatically to ensure that the final iteration
of the reconstructed daily incidence does not differ from the previous
iteration beyond a tolerance of 10$`^{-6}`$ by default. If convergence
is not reached, a warning will be returned and the algorithm can be run
with more iterations by modifying `iter =`. The tolerance threshold can
also be modified using `tol =`.

## FAQs

- **Why does R_(t) estimation start later than the start date of the
  incidence data supplied?**

In order to reconstruct daily incidence data, the method requires that
R_(t) is estimated for each aggregation window in turn, which is
translated into a growth rate and used to reconstruct daily incidence
assuming exponential growth. As there is no past incidence data beyond
the first aggregation window, R_(t) cannot be estimated and the daily
incidence cannot be reconstructed until the first day of the second
aggregation window.

Additionally, R_(t) estimation will not start until case numbers are
sufficiently high.

**Please also see the FAQ section in the main “EpiEstim Vignette”.**

### References

Nash RK, Cori A, Nouvellet P. Estimating the epidemic reproduction
number from temporally aggregated incidence data: a statistical
modelling approach and software tool. medRxiv pre-print.

Bi Q, et al. Epidemiology and transmission of COVID-19 in 391 cases and
1286 of their close contacts in Shenzhen, China: a retrospective cohort
study. Lancet. 2020.

Wallinga J, Lipsitch M. How generation intervals shape the relationship
between growth rates and reproductive numbers. Proceedings of the Royal
Society B: Biological Sciences. 2007 Feb 22;274(1609):599–604.
