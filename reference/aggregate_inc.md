# Aggregating daily incidence to longer time windows

Aggregating daily incidence to longer time windows

## Usage

``` r
aggregate_inc(incid, dt = 7L)
```

## Arguments

- incid:

  a vector of daily incidence values

- dt:

  a positive integer, or vector thereof, indicating the length(s) of the
  desired aggregation window(s). If a vector, this will be recycled. For
  example, `dt = c(3L, 4L)` would correspond to alternating aggregation
  windows of 3 and 4 days

## Value

a vector of incidence values, aggregated to dt

## Examples

``` r
## Constant aggregation e.g. weekly reporting
data("SARS2003")
incid <- SARS2003$incidence
dt <- 7L
aggregate_inc(incid, dt)
#> Incidence aggregated up to day 105 of 107
#>  [1]   4   7  24 149 158 430 225 180 105  63  47  33  23  11   6

## Non-constant aggregation e.g. reporting 3x week
data("SARS2003")
incid <- SARS2003$incidence
dt <- c(2L,2L,3L)
aggregate_inc(incid, dt)
#> Incidence aggregated up to day 107 of 107
#>  [1]   1   1   2   2   3   2   0   5  19  58  38  53  51  54  53  56 199 175  58
#> [20]  80  87  68  56  56  33  30  42  20  12  31  17  16  14  10   9  14   5   8
#> [39]  10   3   1   7   2   2   2   2
```
