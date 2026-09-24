# Plot epsilon traces and posterior R means with 95% credible intervals.
# Returns a named list of ggplot objects for the requested components.
plot.estimate_advantage <- function(x,
                                    what = c("all", "epsilon", "R", "diagnostics"),
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
  if (what == "diagnostics") {
    stop("Diagnostic plots are not yet implemented.")
  }
  theme_epiestim <- if (plot_theme == "v2") theme_epiestimv2 else theme
  plots <- list()

  if (what %in% c("all", "epsilon")) {
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
      variant = factor(rep(variant_names, times = n_samples),
                       levels = variant_names),
      sample = rep(seq_len(n_samples), each = n_variants),
      epsilon = as.vector(epsilon)
    )

    if (is.null(options_epsilon$col)) {
      options_epsilon$col <- grDevices::palette()
    }
    if (is.null(options_epsilon$transp)) options_epsilon$transp <- 0.7
    plots$epsilon <- ggplot(
      epsilon_long, aes(sample, epsilon, colour = variant)
    ) +
      geom_line(alpha = options_epsilon$transp) +
      scale_colour_manual(values = rep_len(options_epsilon$col, n_variants)) +
      coord_cartesian(xlim = options_epsilon$xlim,
                              ylim = options_epsilon$ylim) +
      labs(x = "Retained MCMC sample", y = "epsilon", colour = "Variant") +
      theme_epiestim()
  }

  repro_num <- x$R
  if (!is.array(repro_num) || length(dim(repro_num)) != 3L || !is.numeric(repro_num) ||
      any(dim(repro_num) == 0L) || any(is.infinite(repro_num))) {
      stop("x$R must be a non-empty numeric array: time, location, MCMC sample.")
  }
  n_times <- dim(repro_num)[1L]
  n_locations <- dim(repro_num)[2L]
  location_names <- dimnames(repro_num)[[2L]]
  if (is.null(location_names)) {
    location_names <- paste("Location", seq_len(n_locations))
  }

  # Summarise over the third dimension, retaining time and location.
  # Times outside the estimation window can contain only missing draws.
  summaries <- apply(repro_num, c(1L, 2L), function(draws) {
    draws <- draws[!is.na(draws)]
    if (!length(draws)) return(c(mean = NA_real_, lower = NA_real_, upper = NA_real_))
    c(mean = mean(draws),
      lower = unname(stats::quantile(draws, 0.025)),
      upper = unname(stats::quantile(draws, 0.975)))
  })
  repro_num_summary <- data.frame(
    time = rep(seq_len(n_times), times = n_locations),
    location = factor(rep(location_names, each = n_times),
                      levels = location_names),
    mean = as.vector(summaries[1L, , ]),
    lower = as.vector(summaries[2L, , ]),
    upper = as.vector(summaries[3L, , ])
  )

  if (is.null(options_repro_num$col)) options_repro_num$col <- grDevices::palette()
  if (is.null(options_repro_num$transp)) options_R$transp <- 0.2
  colours <- rep_len(options_R$col, n_locations)
  plots$R <- ggplot(
    repro_num_summary, aes(time, mean, colour = location, fill = location)
  ) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                         alpha = options_R$transp, colour = NA, na.rm = TRUE) +
    geom_line(na.rm = TRUE) +
    geom_hline(yintercept = 1, linetype = "dotted") +
    scale_colour_manual(values = colours) +
    scale_fill_manual(values = colours) +
    coord_cartesian(xlim = options_R$xlim, ylim = options_R$ylim) +
    labs(x = "Time", y = "R", colour = "Location", fill = "Location") +
    theme_epiestim()
  

  plots
}
