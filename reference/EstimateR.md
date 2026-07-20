# Function to ensure compatibility with EpiEstim versions \<2.0

Please only use for compatibility; Prefer the new
[`estimate_R()`](https://mrc-ide.github.io/EpiEstim/reference/estimate_R.md)
function instead.

## Usage

``` r
EstimateR(
  I,
  T.Start,
  T.End,
  method = c("NonParametricSI", "ParametricSI", "UncertainSI"),
  n1 = NULL,
  n2 = NULL,
  Mean.SI = NULL,
  Std.SI = NULL,
  Std.Mean.SI = NULL,
  Min.Mean.SI = NULL,
  Max.Mean.SI = NULL,
  Std.Std.SI = NULL,
  Min.Std.SI = NULL,
  Max.Std.SI = NULL,
  SI.Distr = NULL,
  Mean.Prior = 5,
  Std.Prior = 5,
  CV.Posterior = 0.3,
  plot = FALSE,
  leg.pos = "topright"
)
```

## Arguments

- I:

  see `incid` in
  [`estimate_R()`](https://mrc-ide.github.io/EpiEstim/reference/estimate_R.md)

- T.Start:

  see `config$t_start` from
  [`make_config()`](https://mrc-ide.github.io/EpiEstim/reference/make_config.md)
  used in
  [`estimate_R()`](https://mrc-ide.github.io/EpiEstim/reference/estimate_R.md)

- T.End:

  see `config$t_end` from
  [`make_config()`](https://mrc-ide.github.io/EpiEstim/reference/make_config.md)
  used in
  [`estimate_R()`](https://mrc-ide.github.io/EpiEstim/reference/estimate_R.md)

- method:

  see `method` in
  [`estimate_R()`](https://mrc-ide.github.io/EpiEstim/reference/estimate_R.md)
  (but `EstimateR()` uses CamelCase where
  [`estimate_R()`](https://mrc-ide.github.io/EpiEstim/reference/estimate_R.md)
  uses snake_case for the method names)

- n1:

  see `n1` in
  [`estimate_R()`](https://mrc-ide.github.io/EpiEstim/reference/estimate_R.md)

- n2:

  see `n2` in
  [`estimate_R()`](https://mrc-ide.github.io/EpiEstim/reference/estimate_R.md)

- Mean.SI:

  see `config$mean_si` from
  [`make_config()`](https://mrc-ide.github.io/EpiEstim/reference/make_config.md)
  used in
  [`estimate_R()`](https://mrc-ide.github.io/EpiEstim/reference/estimate_R.md)

- Std.SI:

  see `config$std_si` from
  [`make_config()`](https://mrc-ide.github.io/EpiEstim/reference/make_config.md)
  used in
  [`estimate_R()`](https://mrc-ide.github.io/EpiEstim/reference/estimate_R.md)

- Std.Mean.SI:

  see `config$std_mean_si` from
  [`make_config()`](https://mrc-ide.github.io/EpiEstim/reference/make_config.md)
  used in
  [`estimate_R()`](https://mrc-ide.github.io/EpiEstim/reference/estimate_R.md)

- Min.Mean.SI:

  see `config$min_mean_si` from
  [`make_config()`](https://mrc-ide.github.io/EpiEstim/reference/make_config.md)
  used in
  [`estimate_R()`](https://mrc-ide.github.io/EpiEstim/reference/estimate_R.md)

- Max.Mean.SI:

  see `config$max_mean_si` from
  [`make_config()`](https://mrc-ide.github.io/EpiEstim/reference/make_config.md)
  used in
  [`estimate_R()`](https://mrc-ide.github.io/EpiEstim/reference/estimate_R.md)

- Std.Std.SI:

  see `config$std_std_si` from
  [`make_config()`](https://mrc-ide.github.io/EpiEstim/reference/make_config.md)
  used in
  [`estimate_R()`](https://mrc-ide.github.io/EpiEstim/reference/estimate_R.md)

- Min.Std.SI:

  see `config$min_std_si` from
  [`make_config()`](https://mrc-ide.github.io/EpiEstim/reference/make_config.md)
  used in
  [`estimate_R()`](https://mrc-ide.github.io/EpiEstim/reference/estimate_R.md)

- Max.Std.SI:

  see `config$max_std_si` from
  [`make_config()`](https://mrc-ide.github.io/EpiEstim/reference/make_config.md)
  used in
  [`estimate_R()`](https://mrc-ide.github.io/EpiEstim/reference/estimate_R.md)

- SI.Distr:

  see `config$si_distr` from
  [`make_config()`](https://mrc-ide.github.io/EpiEstim/reference/make_config.md)
  used in
  [`estimate_R()`](https://mrc-ide.github.io/EpiEstim/reference/estimate_R.md)

- Mean.Prior:

  see `config$mean_prior` from
  [`make_config()`](https://mrc-ide.github.io/EpiEstim/reference/make_config.md)
  used in
  [`estimate_R()`](https://mrc-ide.github.io/EpiEstim/reference/estimate_R.md)

- Std.Prior:

  see `config$std_prior` from
  [`make_config()`](https://mrc-ide.github.io/EpiEstim/reference/make_config.md)
  used in
  [`estimate_R()`](https://mrc-ide.github.io/EpiEstim/reference/estimate_R.md)

- CV.Posterior:

  see `config$cv_posterior` from
  [`make_config()`](https://mrc-ide.github.io/EpiEstim/reference/make_config.md)
  used in
  [`estimate_R()`](https://mrc-ide.github.io/EpiEstim/reference/estimate_R.md)

- plot:

  Not used anymore, only there for compatibility

- leg.pos:

  Not used anymore, only there for compatibility
