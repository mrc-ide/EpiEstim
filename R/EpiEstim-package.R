#' EpiEstim: Estimate Time Varying Reproduction Numbers from Epidemic Curves
#'
#' The EpiEstim package provides tools to Tools to quantify transmissibility throughout
#' an epidemic from the analysis of time series of incidence.
#' 
#' @section Bibliography:
#' A BibTeX file of  papers that describe the methodology underlying EpiEstim is available via
#' \code{system.file("epiestimpapers.bib", package = "EpiEstim")}.
#'
#' @importFrom ggplot2 last_plot ggplot aes geom_step ggtitle
#' @importFrom ggplot2 geom_ribbon geom_line xlab ylab xlim geom_hline ylim geom_histogram
#' @importFrom ggplot2 scale_colour_manual scale_fill_manual scale_linetype_manual lims theme
#' @importFrom ggplot2 margin element_rect theme_light %+replace% element_blank element_line
#' @importFrom ggplot2 element_text scale_y_continuous
#' 
#' @importFrom rlang .data
#' 
#' @keywords internal
"_PACKAGE"
