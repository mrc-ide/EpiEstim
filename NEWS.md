# EpiEstim 2.2-0

## MISC

* Added a `NEWS.md` file to track changes to the package. (#74, @zkamvar)
* Added tests for plotting with vdiffr. (#74, @zkamvar)
* Remove un-used dependencies that were added during hackout3: plyr, grid, 
  and plotly (#74, @zkamvar)
* Remove compare package from suggests and use test_that version. (#74, @zkamvar)
* Bump minimum required version of coarseDataTools to 0.6-4 (#71, @zkamvar)
* Incidence objects are now handled appropriately with accessors (#65, @zkamvar)

# EpiEstim 2.1-0

* TBD

# EpiEstim 2.0-0

* Changed function names to snake_case (only exception is that R remains capital letter to avoid confusion between the reproduction number R and the growth rate 
r) and to be more explicit; so `EstimateR` becomes `estimate_R`, `OverallInfectivity` becomes `oberall_infectivity`, `WT` becomes `wallinga_teunis`, and `DiscrSI` becomes `discr_si`
* Additional methods available for estimate_R: in addition to `non_parametric_si`, `parametric_si` and `uncertain_si`, which were already available in EpiEstim 1.0-0, two new methods have been added: `si_from_data` or `si_from_sample`. These allow feeding function `estimate_R` data on observed serial intervals (method `si_from_data`) or posterior samples of serial interval distributions obtained from such data (method `si_from_sample`). 
