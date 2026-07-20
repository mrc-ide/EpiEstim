# Function to ensure compatibility with EpiEstim versions \<2.0

Please only use for compatibility; Prefer the new
[`wallinga_teunis()`](https://mrc-ide.github.io/EpiEstim/reference/wallinga_teunis.md)
function instead

## Usage

``` r
WT(
  I,
  T.Start,
  T.End,
  method = c("NonParametricSI", "ParametricSI"),
  Mean.SI = NULL,
  Std.SI = NULL,
  SI.Distr = NULL,
  nSim = 10,
  plot = FALSE,
  leg.pos = "topright"
)
```

## Arguments

- I:

  see `incid` in
  [`wallinga_teunis()`](https://mrc-ide.github.io/EpiEstim/reference/wallinga_teunis.md)

- T.Start:

  see `config$t_start` in
  [`wallinga_teunis()`](https://mrc-ide.github.io/EpiEstim/reference/wallinga_teunis.md)

- T.End:

  see `config$t_end` in
  [`wallinga_teunis()`](https://mrc-ide.github.io/EpiEstim/reference/wallinga_teunis.md)

- method:

  see `method` in
  [`wallinga_teunis()`](https://mrc-ide.github.io/EpiEstim/reference/wallinga_teunis.md)
  (but WT uses CamelCase where wallinga_teunis uses snake_case for the
  method names)

- Mean.SI:

  see `config$mean_si` in
  [`wallinga_teunis()`](https://mrc-ide.github.io/EpiEstim/reference/wallinga_teunis.md)

- Std.SI:

  see `config$std_si` in
  [`wallinga_teunis()`](https://mrc-ide.github.io/EpiEstim/reference/wallinga_teunis.md)

- SI.Distr:

  see `config$si_distr` in
  [`wallinga_teunis()`](https://mrc-ide.github.io/EpiEstim/reference/wallinga_teunis.md)

- nSim:

  see `config$n_sim` in
  [`wallinga_teunis()`](https://mrc-ide.github.io/EpiEstim/reference/wallinga_teunis.md)

- plot:

  Not used anymore, only there for compatibility

- leg.pos:

  Not used anymore, only there for compatibility
