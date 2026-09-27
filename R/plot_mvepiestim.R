plot_r_trace <- function(rmat, options_r) {

  if (!is.array(rmat) || length(dim(rmat)) != 3L || !is.numeric(rmat) ||
    any(dim(rmat) == 0L) || any(is.infinite(rmat))) {
    stop("x$R must be a non-empty numeric array: time, location, MCMC sample.")
  }
  n_locations <- dim(rmat)[2L]
  n_samples <- dim(rmat)[3L]
  location_names <- dimnames(rmat)[[2L]]
  if (is.null(location_names)) {
    location_names <- paste("Location", seq_len(n_locations))
  }

  # Matrices are flattened column by column: locations vary within each sample.
  r_long <- data.frame(
    location = factor(rep(location_names, times = n_samples),
      levels = location_names
    ),
    sample = rep(seq_len(n_samples), each = n_locations),
    R = as.vector(rmat)
  )

  if (is.null(options_r$col)) options_r$col <- grDevices::palette()
  if (is.null(options_r$transp)) options_r$transp <- 0.2
  p <- ggplot(
    r_long, aes(sample, R, colour = location)
  ) +
    geom_line(alpha = options_r$transp) +
    scale_colour_manual(values = rep_len(options_r$col, n_locations)) +
    ggplot2::coord_cartesian(
      xlim = options_r$xlim,
      ylim = options_r$ylim
    ) +
    labs(x = "Retained MCMC sample", y = "R", colour = "Location")

  p
}

plot_epsilon_trace <- function(epsilon, options_epsilon) {
  
  if (!is.matrix(epsilon) || !is.numeric(epsilon) ||
      any(dim(epsilon) == 0L) || any(!is.finite(epsilon))) {
    stop("x$epsilon must be a non-empty numeric matrix with finite values.")
  }

  n_variants <- nrow(epsilon)
  n_samples <- ncol(epsilon)
  variant_names <- rownames(epsilon)
  if (is.null(variant_names)) {
    variant_names <- paste("Variant", seq_len(n_variants))
  }

  # Matrices are flattened column by column: variants vary within each sample.
  epsilon_long <- data.frame(
    variant = factor(rep(variant_names, times = n_samples),
                     levels = variant_names),
    sample = rep(seq_len(n_samples), each = n_variants),
    epsilon = as.vector(epsilon)
  )

  if (is.null(options_epsilon$col)) {
    options_epsilon$col <- grDevices::palette()
  }
  if (is.null(options_epsilon$transp)) options_epsilon$transp <- 0.7
  p <- ggplot(
    epsilon_long, aes(sample, epsilon, colour = variant)
  ) +
    geom_line(alpha = options_epsilon$transp) +
    scale_colour_manual(values = rep_len(options_epsilon$col, n_variants)) +
    coord_cartesian(xlim = options_epsilon$xlim,
                            ylim = options_epsilon$ylim) +
    labs(x = "Retained MCMC sample", y = "epsilon", colour = "Variant")

  p
}

#' Plot epsilon traces and posterior R means with 95% credible intervals.
#' Returns a named list of ggplot objects for the requested components.
#' @export
#' @method plot estimate_advantage
plot.estimate_advantage <- function(x,
                                    what = c("all", "epsilon", "R"),
                                    r_trace_at = NULL,
                                    plot_theme = "v2",
                                    options_epsilon = list(
                                      col = grDevices::palette(),
                                      transp = 0.7,
                                      xlim = NULL,
                                      ylim = NULL
                                    ),
                                    options_R = list(
                                      col = grDevices::palette(),
                                      transp = 0.2,
                                      xlim = NULL,
                                      ylim = NULL
                                    ), ...) {
  what <- match.arg(what)

  theme_epiestim <- if (plot_theme == "v2") theme_epiestimv2 else theme
  plots <- list()

  if (what %in% c("all", "epsilon")) {
    pepsilon <- plot_epsilon_trace(x$epsilon, options_epsilon)
    plots$epsilon <- pepsilon + theme_epiestim()
  }

  if (what %in% c("all", "R")) {
    if (is.null(r_trace_at)) {
      warning("r_trace_at is NULL, plotting trace at the last time point.", call. = FALSE)
      r_trace_at <- dim(x$R)[1L]
    }
    if (!is.null(r_trace_at)) {
      prtrace <- plot_r_trace(x$R[r_trace_at, , , drop = FALSE], options_R)
      plots$R <- prtrace + theme_epiestim()
    }
  }

  plots
}
