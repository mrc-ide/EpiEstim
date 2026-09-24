plot.estimate_advantage <- function(x,
                                    what = c("all", "epsilon", "R", "diagnostics"),
                                    plot_theme = "v2",
                                    options_epsilon = list(
                                      col = grDevices::palette(),
                                      transp = 0.7,
                                      xlim = NULL,
                                      ylim = NULL
                                    ), ...) {
  what <- match.arg(what)
  if (!what %in% c("all", "epsilon")) {
    stop("Only epsilon plots are currently implemented.")
  }

  epsilon <- x$epsilon
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
    variant = rep(variant_names, times = n_samples),
    sample = rep(seq_len(n_samples), each = n_variants),
    epsilon = as.vector(epsilon)
  )

  if (is.null(options_epsilon$col)) {
    options_epsilon$col <- grDevices::palette()
  }
  if (is.null(options_epsilon$transp)) options_epsilon$transp <- 0.7
  if (is.null(options_epsilon$xlim)) {
    options_epsilon$xlim <- range(epsilon_long$sample)
  }
  if (is.null(options_epsilon$ylim)) {
    options_epsilon$ylim <- range(epsilon_long$epsilon)
  }
  colours <- grDevices::adjustcolor(
    rep_len(options_epsilon$col, length.out = n_variants),
    alpha.f = options_epsilon$transp
  )

  theme_epiestim <- ifelse(plot_theme == "v2", theme_epiestimv2, theme)
  pmcmc_samples <- ggplot() +
    geom_line(aes(sample, epsilon, color = variant)) +
    theme_epiestim()
  
  list(pmcmc_samples)
}
