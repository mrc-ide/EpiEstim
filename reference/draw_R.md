# Draw R from marginal posterior distribution

Draw R from marginal posterior distribution

## Usage

``` r
draw_R(
  epsilon,
  incid,
  lambda,
  priors,
  shape_R_flat = NULL,
  t_min = NULL,
  t_max = nrow(incid),
  seed = NULL
)
```

## Arguments

- epsilon:

  a value or vector of values for the relative transmissibility of the
  "new" pathogen/strain/variant(s) compared to the reference
  pathogen/strain/variant

- incid:

  a multidimensional array containing values of the (local) incidence
  for each time step (1st dimension), location (2nd dimension) and
  pathogen/strain/variant (3rd dimension)

- lambda:

  a multidimensional array containing values of the overall infectivity
  for each time step (1st dimension), location (2nd dimension) and
  pathogen/strain/variant (3rd dimension). The overall infectivity for a
  given location and pathogen/strain/variant represents the sum of the
  incidence for that location and that pathogen/strain/variant at all
  previous time steps, weighted by the current infectivity of those past
  incident cases. It can be calculated from the incidence `incid` and
  the distribution of the serial interval using function
  [`compute_lambda()`](https://mrc-ide.github.io/EpiEstim/reference/compute_lambda.md)

- priors:

  a list of prior parameters (shape and scale of a gamma distribution)
  for epsilon and R; can be obtained from the function
  [`default_priors()`](https://mrc-ide.github.io/EpiEstim/reference/default_priors.md).
  The prior for R is assumed to be the same for all time steps and all
  locations

- shape_R_flat:

  a vector of the shape of the posterior distribution of R for each time
  step t and each location l (stored in element
  `(l-1)*(t_max - t_min + 1) + t` of the vector), as obtained from
  function
  [`get_shape_R_flat()`](https://mrc-ide.github.io/EpiEstim/reference/get_shape_R_flat.md).

- t_min:

  an integer \> 1 giving the minimum time step to consider in the
  estimation. Default value is 2 (as the estimation is conditional on
  observations at time step 1 and can therefore only start at time step
  2).

- t_max:

  an integer \> `t_min` and \<= `nrow(incid)` giving the maximum time
  step to consider in the estimation. Default value is `nrow(incid)`.

- seed:

  a numeric value used to fix the random seed

## Value

a matrix of the instantaneous reproduction number R for the reference
pathogen/strain/variant for each time step (row) and each location
(column) drawn from the marginal posterior distribution

## Examples

``` r
n_v <- 2
n_loc <- 3 # 3 locations
T <- 100 # 100 time steps
priors <- default_priors()
# constant incidence 10 per day everywhere
incid <- array(10, dim = c(T, n_loc, n_v))
incid <- process_I_multivariant(incid)
# arbitrary serial interval, same for both variants
w_v <- c(0, 0.2, 0.5, 0.3)
si_distr <- cbind(w_v, w_v)
lambda <- compute_lambda(incid, si_distr)
# Epsilon = 1 i.e. no transmission advantage
epsilon <- 1
draw_R(epsilon, incid$local, lambda, priors, seed = 1, t_min = 2L)
#>             [,1]      [,2]      [,3]
#>   [1,]        NA        NA        NA
#>   [2,] 4.1754778 6.4177517 5.2923993
#>   [3,] 1.8419050 1.3249988 1.3698939
#>   [4,] 1.2759194 0.7448514 1.2165933
#>   [5,] 1.0686560 1.1385853 0.9931790
#>   [6,] 0.6649529 0.9588906 0.9569318
#>   [7,] 1.0855305 0.7220602 1.1054765
#>   [8,] 1.1447094 0.9048344 1.0594580
#>   [9,] 1.1061911 0.8729437 0.8835008
#>  [10,] 0.9088510 1.0958110 1.1002343
#>  [11,] 0.8067693 0.6688852 0.7408918
#>  [12,] 0.8428318 0.7756092 1.1733477
#>  [13,] 0.9101662 0.9097909 1.3912541
#>  [14,] 0.9651637 0.7926315 1.0355282
#>  [15,] 0.7882348 0.9625400 0.8841484
#>  [16,] 1.1646074 0.9902207 0.6335370
#>  [17,] 1.1104524 0.9049448 0.7787620
#>  [18,] 1.1882926 0.8754851 1.3867636
#>  [19,] 1.1552044 0.8130454 0.8797009
#>  [20,] 0.9915667 1.4897800 1.6469520
#>  [21,] 0.5856106 0.9788908 1.0097665
#>  [22,] 1.0892993 0.7119596 1.2402854
#>  [23,] 1.2292548 1.3263248 0.5354871
#>  [24,] 0.8285957 0.7072327 0.9102313
#>  [25,] 0.8724321 0.9061527 1.1884938
#>  [26,] 0.9237510 1.1813990 1.0648465
#>  [27,] 0.9525089 0.9202169 0.8872294
#>  [28,] 0.8370521 1.2081112 1.0916053
#>  [29,] 0.6949556 0.8428264 1.0582477
#>  [30,] 1.2061801 0.8941177 1.3923786
#>  [31,] 0.8826097 1.0710726 0.9301088
#>  [32,] 1.2327881 0.9231198 1.1946829
#>  [33,] 1.1506565 1.3809056 1.0730621
#>  [34,] 0.9390971 1.1545916 1.2093712
#>  [35,] 0.8777327 1.1387221 0.8908964
#>  [36,] 1.1017039 0.8589781 0.6169645
#>  [37,] 0.8290427 0.6861502 0.6858241
#>  [38,] 0.6606188 0.7653708 1.1000737
#>  [39,] 1.1519406 1.5505936 0.8274025
#>  [40,] 0.9504261 0.9985273 1.1971140
#>  [41,] 1.0346379 1.0784597 0.9639281
#>  [42,] 0.8937274 0.9581058 0.9477953
#>  [43,] 1.0517454 1.1116471 0.7577208
#>  [44,] 0.7418473 0.6679318 1.4543461
#>  [45,] 1.0310064 1.4865314 0.9166080
#>  [46,] 1.1161898 1.2148392 1.3371404
#>  [47,] 0.9373445 1.2596910 0.9934363
#>  [48,] 1.1047674 0.7223597 1.1041808
#>  [49,] 0.9454872 0.8487927 0.7621509
#>  [50,] 0.9726150 0.6782616 0.6581815
#>  [51,] 1.1331269 0.9551992 0.9970270
#>  [52,] 0.9812364 1.3248666 0.8774496
#>  [53,] 0.8179912 0.8133895 1.0525876
#>  [54,] 0.6860385 1.2333217 1.3442366
#>  [55,] 1.3251154 0.9713591 0.9037205
#>  [56,] 1.0091474 1.0657415 0.9485403
#>  [57,] 1.5131681 0.8203213 0.6695261
#>  [58,] 1.0827581 1.1984967 1.2848857
#>  [59,] 0.8247382 1.4099866 1.3332242
#>  [60,] 0.9489409 0.7544445 1.1630366
#>  [61,] 0.7181302 0.7638487 0.6062275
#>  [62,] 1.0070208 0.8697879 1.3527410
#>  [63,] 0.9492489 1.0674464 1.2099130
#>  [64,] 0.9915170 1.3831627 0.9157278
#>  [65,] 0.8493494 1.3564213 0.7937749
#>  [66,] 1.2712944 0.9034245 1.2361158
#>  [67,] 1.1528061 1.3488304 1.1685096
#>  [68,] 0.9272305 0.8961155 1.0456598
#>  [69,] 1.1104629 1.0981111 1.0246826
#>  [70,] 1.0498748 0.9720965 0.7978442
#>  [71,] 1.2236467 1.1900186 1.1433775
#>  [72,] 0.9091075 0.9832146 1.1950250
#>  [73,] 1.1229097 0.7540497 0.9760204
#>  [74,] 1.2827117 0.6962321 0.8988835
#>  [75,] 0.7919397 1.2507439 0.8486861
#>  [76,] 1.2478087 0.6800184 0.7545819
#>  [77,] 1.1356190 0.7181952 0.9782834
#>  [78,] 1.3564850 0.6496650 0.9343853
#>  [79,] 1.1021314 0.6302136 0.6890480
#>  [80,] 0.7137906 1.0374863 1.0307920
#>  [81,] 1.0262212 0.9327329 1.0082427
#>  [82,] 1.2072557 0.7366215 0.7741553
#>  [83,] 1.1551378 0.9434132 1.3272573
#>  [84,] 0.8112364 0.6530636 1.1301635
#>  [85,] 0.7844711 1.5287749 1.0606628
#>  [86,] 0.7416506 0.6695974 0.9329863
#>  [87,] 1.4038401 1.3822027 1.1024619
#>  [88,] 1.1395487 1.0876609 0.9585883
#>  [89,] 1.1861500 1.3852849 0.9410206
#>  [90,] 1.0616345 0.9183850 0.7569440
#>  [91,] 1.3814038 1.2319826 0.4870707
#>  [92,] 0.8398616 0.9738715 1.0986567
#>  [93,] 0.9806569 0.7642993 1.2182067
#>  [94,] 0.9811198 1.0249681 1.0191023
#>  [95,] 0.6396538 0.7908576 0.6487266
#>  [96,] 0.8903295 0.7461245 1.1304180
#>  [97,] 1.0504853 0.8585829 1.1043524
#>  [98,] 1.0872457 0.9527889 0.8528483
#>  [99,] 1.0057392 1.2131339 0.6324443
#> [100,] 0.9490229 1.0053291 1.0373155
```
