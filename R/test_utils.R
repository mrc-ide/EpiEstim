#' Expect doppelgangers skipped on CI if missing snapshot
#' 
#' See `tests/testthat/test-plots.R`.
#'
#' @param name Character. Snapshot name.
#' @param code Code to create image to snapshot.
#' @param variant Character. Variant for the snapshot
#' 
#' @returns snapshot or skip
#' 
#' @noRd

epi_expect_doppelganger <- function(name, code, variant) {

  testthat::skip_if_not_installed("vdiffr")

  on_ci <- as.logical(Sys.getenv("CI", "false"))
  no_snap <- !file.exists(
    testthat::test_path(
      "_snaps",
      variant,
      paste0(name, ".svg")
    )
  )
  if (on_ci && no_snap) {
    testthat::skip(
      paste0("On CI and no snapshot for this variant (", variant, ")")
    )
  }

  vdiffr::expect_doppelganger(
    name,
    code,
    variant = variant
  )
}

