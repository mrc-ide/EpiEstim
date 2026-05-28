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

#' Custom snapshot utility function
#' 
#' Creates mini snapshots by sampling `x` down to `n`. If `x` is a list, the
#' function recurses down the list sampling applies to all sub-items (which may
#' or may not be appropriate, depending on the situation). Also rounds to 8
#' digits to avoid mismatches in re-loading values with long trailing decimals.
#'
#' @param x Vector or list of vectors to take a snapshot of.
#' @param style Snapshot style (see `?testthat::snapshot_value()`)
#' @param n Numeric. Number of samples to keep
#'
#' @noRd
epi_snapshot_value <- function(x, style = "json2", n = 5) {
  set.seed(1) # TODO: could use withr::with_seed()

  sample_x <- function(xx) {
    # Recurse if this is still a list
    if(is.list(xx)) return(lapply(xx, sample_x))
    
    # Take random sample of points or max number
    if(!is.null(n)) {
      nn <- if(length(xx) < n) length(xx) else n
      xx <- sample(xx, size = nn)
    }

    # Round to avoid mismatches from truncated records
    if(is.numeric(xx)) xx <- round(xx, digits = 8) 
    xx
  }

  x <- sample_x(x)
  
  testthat::expect_snapshot_value(x, style = style)
}
