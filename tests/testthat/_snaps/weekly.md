# you can't run estimate_R_agg unless using parametric/non-parametric SI methods

    Code
      estimate_R_agg(incid = weekly_inc, dt = 7L, dt_out = 7L, iter = 10L, config = config,
        method = method, grid = list(precision = 0.001, min = -1, max = 1))
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "non_parametric_si", "parametric_si"

