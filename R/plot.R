#' Plot outputs of estimate_r
#'
#' The plot method of \code{estimate_r} objects can be used to visualise three
#' types of information. The first one shows the epidemic curve. The second one
#' shows the posterior mean and 95\% credible interval of the reproduction
#' number. The estimate for a time window is plotted at the end of the time
#' window. The third plot shows the discrete distribution(s) of the serial
#' interval.
#'
#'
#' @param x The output of function \code{\link{estimate_R}} or function
#'   \code{\link{wallinga_teunis}}. To plot simultaneous outputs on the same
#'   plot use \code{\link{estimate_R_plots}} function
#'
#' @param what A string specifying what to plot, namely the incidence time
#'   series (\code{what='incid'}), the estimated reproduction number
#'   (\code{what='R'}), the serial interval distribution (\code{what='SI'}, or
#'   all three (\code{what='all'})).
#'
#' @param plot_theme A string specifying whether to use the original plot theme
#' (plot_theme = "original") or an alternative plot theme (plot_theme = "v2").
#' The plot_theme is "v2" by default.
#'
#' @param add_imported_cases A boolean to specify whether, on the incidence time
#'   series plot, to add the incidence of imported cases.
#'
#' @param options_I For what = "incid" or "all". A list of graphical options:
#'   \describe{ \item{col}{A color or vector of colors used for plotting
#'   incid. By default uses the default R colors.}  \item{transp}{A numeric
#'   value between 0 and 1 used to monitor transparency of the bars
#'   plotted. Defaults to 0.7.}  \item{xlim}{A parameter similar to that in
#'   \code{par}, to monitor the limits of the horizontal axis} \item{ylim}{A
#'   parameter similar to that in \code{par}, to monitor the limits of the
#'   vertical axis} \item{interval}{An integer or character indicating the
#'   (fixed) size of the time interval used for plotting the incidence;
#'   defaults to 1 day.} \item{xlab, ylab}{Labels for the axes of the
#'   incidence plot}}
#'
#' @param options_R For what = "R" or "all". A list of graphical options:
#'   \describe{ \item{col}{A color or vector of colors used for plotting R. By
#'   default uses the default R colors.}  \item{transp}{A numeric value between
#'   0 and 1 used to monitor transparency of the 95\%CrI. Defaults to 0.2.}
#'   \item{xlim}{A parameter similar to that in \code{par}, to monitor the
#'   limits of the horizontal axis} \item{ylim}{A parameter similar to that in
#'   \code{par}, to monitor the limits of the vertical axis}
#'   \item{xlab, ylab}{Labels for the axes of the R plot}}
#'
#' @param options_SI For what = "SI" or "all". A list of graphical options:
#'   \describe{ \item{prob_min}{A numeric value between 0 and 1. The SI
#'   distributions explored are only shown from time 0 up to the time t so that
#'   each distribution explored has probability < \code{prob_min} to be on any
#'   time step after t. Defaults to 0.001.}  \item{col}{A color or vector of
#'   colors used for plotting the SI. Defaults to black.}  \item{transp}{A
#'   numeric value between 0 and 1 used to monitor transparency of the
#'   lines. Defaults to 0.25} \item{xlim}{A parameter similar to that in
#'   \code{par}, to monitor the limits of the horizontal axis} \item{ylim}{A
#'   parameter similar to that in \code{par}, to monitor the limits of the
#'   vertical axis} \item{xlab, ylab}{Labels for the axes of the serial interval
#'    distribution plot}}
#'
#' @param legend A boolean (TRUE by default) governing the presence / absence of
#'   legends on the plots
#'
#' @param ... further arguments passed to other methods (currently unused).
#'
#' @return a plot (if \code{what = "incid"}, \code{"R"}, or \code{"SI"}) or a
#'   \code{\link[grid]{grob}} object (if \code{what = "all"}).
#'
#' @seealso \code{\link{estimate_R}},
#'   \code{\link{wallinga_teunis}} and
#'   \code{\link{estimate_R_plots}}
#'
#' @author Rolina van Gaalen \email{rolina.van.gaalen@rivm.nl} and Anne Cori
#'   \email{a.cori@imperial.ac.uk}; S3 method by Thibaut Jombart; v2 theme by
#'   Rebecca Nash
#'
#' @import reshape2 grid gridExtra
#'
#' @importFrom scales alpha
#'
#' @importFrom grDevices palette
#'
#' @importFrom ggplot2 last_plot ggplot aes aes_string geom_step ggtitle
#'   geom_ribbon geom_line xlab ylab xlim geom_hline ylim geom_histogram
#'   scale_colour_manual scale_fill_manual scale_linetype_manual lims theme
#'   margin element_rect theme_light %+replace% element_blank element_line
#'   element_text scale_y_continuous 
#'
#' @importFrom graphics plot
#'
#' @importFrom incidence as.incidence
#'
#' @importFrom patchwork plot_layout
#'
#'
#' @export
#'
#' @examples
#' ## load data on pandemic flu in a school in 2009
#' data("Flu2009")
#'
#' ## estimate the instantaneous reproduction number
#' ## (method "non_parametric_si")
#' R_i <- estimate_R(Flu2009$incidence,
#'   method = "non_parametric_si",
#'   config = list(
#'     t_start = seq(2, 26),
#'     t_end = seq(8, 32),
#'     si_distr = Flu2009$si_distr
#'   )
#' )
#'
#' ## visualise results
#' plot(R_i, legend = FALSE)
#'
#' ## estimate the instantaneous reproduction number
#' ## (method "non_parametric_si")
#' R_c <- wallinga_teunis(Flu2009$incidence,
#'   method = "non_parametric_si",
#'   config = list(
#'     t_start = seq(2, 26),
#'     t_end = seq(8, 32),
#'     si_distr = Flu2009$si_distr,
#'     n_sim = 10
#'   )
#' )
#'
#' ## produce plot of the incidence
#' ## (with, on top of total incidence, the incidence of imported cases),
#' ## estimated instantaneous and case reproduction numbers
#' ## and serial interval distribution used
#' p_I <- plot(R_i, "incid", add_imported_cases=TRUE) # plots the incidence
#' p_SI <- plot(R_i, "SI") # plots the serial interval distribution
#' p_Ri <- plot(R_i, "R",
#'              options_R = list(ylim = c(0, 4)))
#'         # plots the estimated instantaneous reproduction number
#' p_Rc <- plot(R_c, "R",
#'              options_R = list(ylim = c(0, 4)))
#'         # plots the estimated case reproduction number
#' gridExtra::grid.arrange(p_I, p_SI, p_Ri, p_Rc, ncol = 2)
plot.estimate_R <- function(x, what = c("all", "incid", "R", "SI"), plot_theme = "v2",
                            add_imported_cases = FALSE,
                            options_I = list(
                              col = palette(), transp = 0.7,
                              xlim = NULL, ylim = NULL,
                              interval = 1L,
                              xlab = "Time",
                              ylab = "Incidence"
                            ),
                            options_R = list(
                              col = palette(), transp = 0.2,
                              xlim = NULL, ylim = NULL,
                              xlab = "Time",
                              ylab = "R"
                            ),
                            options_SI = list(
                              prob_min = 0.001,
                              col = "black", transp = 0.25,
                              xlim = NULL, ylim = NULL,
                              xlab = "Time",
                              ylab = "Frequency"
                            ),
                            legend = TRUE, ...) {
  ## dealing with the fact that some options may be left to default but others
  ## may have been specified by user
  if (is.null(options_I$col)) options_I$col <- palette()
  if (is.null(options_I$transp)) options_I$transp <- 0.7
  if (is.null(options_I$xlab)) options_I$xlab <- "Time"
  if (is.null(options_I$ylab)) options_I$ylab <- "Incidence"
  if (is.null(options_I$interval)) options_I$interval <- 1L

  if (is.null(options_R$col)) options_R$col <- palette()
  if (is.null(options_R$transp)) options_R$transp <- 0.2
  if (is.null(options_R$xlab)) options_R$xlab <- "Time"
  if (is.null(options_R$ylab)) options_R$ylab <- "R"


  if (is.null(options_SI$prob_min)) options_SI$prob_min <- 0.001
  if (is.null(options_SI$col)) options_SI$col <- "black"
  if (is.null(options_SI$transp)) options_SI$transp <- 0.25
  if (is.null(options_SI$xlab)) options_SI$xlab <- "Time"
  if (is.null(options_SI$ylab)) options_SI$ylab <- "Frequency"

  ## New theme

  if (plot_theme == "v2") {
    theme_epiestim <- function() {
      theme_light() %+replace%

        theme(
          panel.border = element_blank(),
          axis.line = element_line(
            colour = "black",
            size = 0.2
          ),
          plot.title = element_text(
            size = 12,
            hjust = 0,
            vjust = 4,
            margin = margin(5, b = 5, t = 10)
          ),
          axis.title = element_text(
            size = 11
          ),
          axis.text = element_text(
            size = 9
          ),
          axis.text.x = element_text(
            margin = margin(5, b = 10)
          ),
          axis.text.y = element_text(
            margin = margin(5, l = 10, r = 4)
          )
        )
    }

    options_I$col <- "#5983AB"
  } else {
    if (plot_theme == "original") {
      theme_epiestim <- function() {
        theme()
      }
    }
  }


  # check if x is a single output of EpiEstim or a list of such outputs
  if (is.data.frame(x[[1]])) # x is a single output of EpiEstim
    {
      multiple_input <- FALSE
      if (plot_theme == "v2") {
        options_R$col <- "#5983AB"
      } else {
        options_R$col <- options_R$col[1]
      }
    } else {
    multiple_input <- TRUE
    if (length(unique(vapply(x, function(e) nrow(e$R), integer(1)))) > 1) {
      stop("R estimates cannot be plotted simulatneously because
           they are of different sizes, i.e. they were obtained using
           t_start or t_end of different lengths")
    }
    x_list <- x
    x <- x_list[[1]]
    if (length(x_list) > length(col)) {
      warnings("color vector too short, recycling colors.")
      options_R$col <- rep(
        options_R$col,
        ceiling(length(x_list) / length(options_R$col))
      )
      options_R$col <- options_R$col[seq_len(length(x_list))]
    } else {
      options_R$col <- options_R$col[seq_len(length(x_list))]
    }
  }

  t_start <- x$R$t_start
  t_end <- x$R$t_end
  mean_posterior <- x$R[, "Mean(R)"]
  quantile_0.025_posterior <- x$R[, "Quantile.0.025(R)"]
  quantile_0.975_posterior <- x$R[, "Quantile.0.975(R)"]
  method <- x$method
  si_distr <- x$si_distr
  incid <- data.frame(local = x$I_local, imported = x$I_imported)
  T <- nrow(incid)
  if (!is.null(x$dates)) {
    dates <- x$dates
  } else {
    dates <- seq_len(T)
  }

  ########################################################################
  ### these few lines are to make CRAN checks happy with ggplot2... ###
  value <- NULL
  meanR <- NULL
  group <- NULL
  lower <- NULL
  upper <- NULL
  end <- NULL
  ########################################################################

  ## temp fix for estimate_R_agg to be able to plot disaggreated incidence <1
  if (any(rowSums(incid) < 1 & rowSums(incid) > 0)) {
    idx_round <- which(rowSums(incid) < 1 & rowSums(incid) > 0)
    incid[idx_round, ] <- ceiling(incid[idx_round, ])
  }
  ## TODO: change plot to allow for non-integer incidence

  what <- match.arg(what)
  if (what %in% c("incid", "all")) {
    if (add_imported_cases) {
      p1 <- plot(
        as.incidence(incid,
          dates = x$dates,
          interval = options_I$interval
        ),
        ylab = options_I$ylab, xlab = options_I$xlab,
        color = options_I$col, alpha = options_I$transp
      ) +
        theme_epiestim() +
        ggtitle("Epidemic curve")
    } else {
      p1 <- plot(
        as.incidence(rowSums(incid),
          dates = x$dates,
          interval = options_I$interval
        ),
        ylab = options_I$ylab, xlab = options_I$xlab,
        color = options_I$col, alpha = options_I$transp
      ) +
        theme_epiestim() +
        ggtitle("Epidemic curve")
    }

    if (!is.null(options_I$xlim)) {
      p1 <- p1 +
        theme_epiestim() +
        lims(x = options_I$xlim)
    }

    if (!is.null(options_I$ylim)) {
      p1 <- p1 +
        theme_epiestim() +
        lims(y = options_I$ylim)
    }
  }
  if (what %in% c("R", "all")) {
    time.points <- apply(x$R[, c("t_start", "t_end")], 1, function(x) {
      seq(x[1], x[2] - 1)
    })
    if (length(time.points) == length(unique(matrix(time.points, ncol = 1)))) {
      # non sliding windows
      if (!multiple_input) {
        if (is.null(options_R$ylim)) {
          options_R$ylim <- c(0, max(quantile_0.975_posterior, na.rm = TRUE))
        }

        if (is.null(options_R$xlim)) {
          options_R$xlim <- c(min(dates), max(dates) + 1)
        }

        df <- melt(data.frame(
          start = dates[t_start] - 0.5, end = dates[t_end] + 0.5, meanR = mean_posterior,
          lower = quantile_0.025_posterior,
          upper = quantile_0.975_posterior
        ), id.vars = c("meanR", "lower", "upper"))
        df$group <- as.factor(rep(
          seq_len(length(t_start)),
          dim(df)[1] / length(t_start)
        ))

        p2 <- ggplot(df, aes(
          x = value, y = as.numeric(meanR),
          group = as.factor(group)
        )) +
          theme_epiestim() +
          geom_ribbon(aes(ymin = lower, ymax = upper, fill = "95%CrI")) +
          geom_line(aes(y = meanR, colour = "Mean")) +
          xlab(options_R$xlab) +
          ylab(options_R$ylab) +
          xlim(options_R$xlim) +
          ylim(options_R$ylim) +
          scale_colour_manual("", values = options_R$col) +
          scale_fill_manual("", values = alpha(
            options_R$col,
            options_R$transp
          )) +
          ggtitle("Estimated R")
      } else {
        df_tmp <- data.frame(
          start = dates[t_start], end = dates[t_end], meanR = mean_posterior,
          lower = quantile_0.025_posterior,
          upper = quantile_0.975_posterior
        )
        df <- df_tmp
        id_tmp <- c("meanR", "lower", "upper")
        id <- id_tmp

        for (i in seq(2, length(x_list)))
        {
          x2 <- x_list[[i]]
          t_start2 <- x2$R$t_start
          if (!is.null(x2$dates)) {
            dates2 <- x2$dates
          } else {
            dates2 <- seq_len(T)
          }
          mean_posterior2 <- x2$R[, "Mean(R)"]
          quantile_0.025_posterior2 <- x2$R[, "Quantile.0.025(R)"]
          quantile_0.975_posterior2 <- x2$R[, "Quantile.0.975(R)"]
          df_tmp2 <- data.frame(
            start2 = dates2[t_start2], end2 = dates2[t_end],
            meanR2 = mean_posterior2, lower2 = quantile_0.025_posterior2,
            upper2 = quantile_0.975_posterior2
          )
          names(df_tmp2) <- paste0(names(df_tmp), i)
          df <- cbind(df, df_tmp2)
          id_tmp2 <- paste0(id, i)
          id <- c(id, id_tmp2)
        }

        if (is.null(options_R$ylim)) {
          options_R$ylim <- c(0, max(df[, grep("upper", names(df), fixed = TRUE)],
            na.rm = TRUE
          ))
        }

        if (is.null(options_R$xlim)) {
          options_R$xlim <- c(min(dates), max(dates) + 1)
        }

        df <- melt(df, id = id)
        df$group <- as.factor(rep(
          seq_len(length(t_start)),
          dim(df)[1] / length(t_start)
        ))

        p2 <- ggplot(df, aes(
          x = value, y = as.numeric(meanR),
          group = as.factor(group)
        )) +
          theme_epiestim() +
          geom_ribbon(aes(ymin = lower, ymax = upper, fill = "95%CrI")) +
          geom_line(aes(y = meanR, colour = "Mean"))

        for (i in seq(2, length(x_list)))
        {
          p2 <- p2 +
            theme_epiestim() +
            geom_ribbon(aes_string(
              ymin = paste0("lower", i),
              ymax = paste0("upper", i),
              fill = shQuote(paste0("95%CrI", i))
            )) +
            geom_line(aes_string(
              y = paste0("meanR", i),
              colour = shQuote(paste0("Mean", i))
            ))
        }

        p2 <- p2 +
          theme_epiestim() +
          xlab(options_R$xlab) +
          ylab(options_R$ylab) +
          xlim(options_R$xlim) +
          ylim(options_R$ylim) +
          scale_colour_manual("", values = options_R$col) +
          scale_fill_manual("",
            values = alpha(options_R$col, options_R$transp)
          ) +
          ggtitle("Estimated R")
      }
    } else {
      if (!multiple_input) {
        if (is.null(options_R$ylim)) {
          options_R$ylim <- c(0, max(quantile_0.975_posterior, na.rm = TRUE))
        }

        if (is.null(options_R$xlim)) {
          options_R$xlim <- c(min(dates), max(dates) + 1)
        }

        p2 <- ggplot(data.frame(
          start = dates[t_start], end = dates[t_end], meanR = mean_posterior,
          lower = quantile_0.025_posterior,
          upper = quantile_0.975_posterior
        ), aes(end, meanR)) +
          theme_epiestim() +
          geom_ribbon(aes(ymin = lower, ymax = upper, fill = "95%CrI")) +
          geom_line(aes(colour = "Mean")) +
          geom_hline(yintercept = 1, linetype = "dotted") +
          xlab(options_R$xlab) +
          ylab(options_R$ylab) +
          xlim(options_R$xlim) +
          ylim(options_R$ylim) +
          ggtitle("Estimated R") +
          scale_colour_manual("", values = options_R$col) +
          scale_fill_manual("", values = alpha(options_R$col, options_R$transp))
      } else {
        df_tmp <- data.frame(
          start = dates[t_start], end = dates[t_end],
          meanR = mean_posterior, lower = quantile_0.025_posterior,
          upper = quantile_0.975_posterior
        )
        df <- df_tmp

        for (i in seq(2, length(x_list)))
        {
          x2 <- x_list[[i]]
          t_start2 <- x2$R$t_start
          if (!is.null(x2$dates)) {
            dates2 <- x2$dates
          } else {
            dates2 <- seq_len(T)
          }
          mean_posterior2 <- x2$R[, "Mean(R)"]
          quantile_0.025_posterior2 <- x2$R[, "Quantile.0.025(R)"]
          quantile_0.975_posterior2 <- x2$R[, "Quantile.0.975(R)"]
          df_tmp2 <- data.frame(
            start2 = dates2[t_start2], end2 = dates2[t_end],
            meanR2 = mean_posterior2, lower2 = quantile_0.025_posterior2,
            upper2 = quantile_0.975_posterior2
          )
          names(df_tmp2) <- paste0(names(df_tmp), i)
          df <- cbind(df, df_tmp2)
        }

        if (is.null(options_R$ylim)) {
          options_R$ylim <- c(0, max(df[, grep("upper", names(df), fixed = TRUE)],
            na.rm = TRUE
          ))
        }

        if (is.null(options_R$xlim)) {
          options_R$xlim <- c(min(dates), max(dates) + 1)
        }

        p2 <- ggplot(df, aes(end, meanR)) +
          theme_epiestim() +
          geom_ribbon(aes(ymin = lower, ymax = upper, fill = "95%CrI")) +
          geom_line(aes(y = meanR, colour = "Mean"))

        for (i in seq(2, length(x_list)))
        {
          p2 <- p2 +
            theme_epiestim() +
            geom_ribbon(aes_string(
              ymin = paste0("lower", i),
              ymax = paste0("upper", i),
              fill = shQuote(paste0("95%CrI", i))
            )) +
            geom_line(aes_string(
              y = paste0("meanR", i),
              colour = shQuote(paste0("Mean", i))
            ))
        }

        p2 <- p2 +
          theme_epiestim() +
          geom_hline(yintercept = 1, linetype = "dotted") +
          xlab(options_R$xlab) +
          ylab(options_R$ylab) +
          xlim(options_R$xlim) +
          ylim(options_R$ylim) +
          ggtitle("Estimated R") +
          scale_colour_manual("", values = options_R$col) +
          scale_fill_manual("", values = alpha(options_R$col, options_R$transp))
      }
    }
  }
  if (what %in% c("SI", "all")) {
    if (method %in% c("uncertain_si", "si_from_data", "si_from_sample")) {
      tmp <- cumsum(apply(si_distr, 2, max) >= options_SI$prob_min)
      stop_at <- min(which(tmp == tmp[length(tmp)]))

      si_distr_for_plot <- si_distr[, seq_len(stop_at)]

      dataL <- melt(t(si_distr_for_plot))
      dataL$Var1 <- seq(0, (ncol(si_distr_for_plot) - 1))
      p3 <- ggplot(dataL, aes_string(
        x = "Var1", y = "value",
        group = "Var2"
      )) +
        theme_epiestim() +
        geom_line(col = options_SI$col, alpha = options_SI$transp) +
        ggtitle("Explored SI distributions") +
        xlab(options_SI$xlab) +
        ylab(options_SI$ylab)

      if (!is.null(options_SI$xlim)) {
        p3 <- p3 +
          theme_epiestim() +
          lims(x = options_SI$xlim)
      }

      if (!is.null(options_SI$ylim)) {
        p3 <- p3 +
          theme_epiestim() +
          lims(y = options_SI$ylim)
      }
    } else {
      tmp <- cumsum(si_distr >= options_SI$prob_min)
      stop_at <- min(which(tmp == tmp[length(tmp)]))

      si_distr_for_plot <- si_distr[seq_len(stop_at)]

      dataL <- data.frame(
        Times = seq(0, length(si_distr_for_plot) - 1),
        SIDistr = si_distr_for_plot
      )
      p3 <- ggplot(dataL, aes_string(x = "Times", y = "SIDistr")) +
        theme_epiestim() +
        geom_line(col = options_SI$col, alpha = options_SI$transp) +
        ggtitle("Explored SI distribution") +
        xlab(options_SI$xlab) +
        ylab(options_SI$ylab)

      if (!is.null(options_SI$xlim)) {
        p3 <- p3 +
          theme_epiestim() +
          lims(x = options_SI$xlim)
      }

      if (!is.null(options_SI$ylim)) {
        p3 <- p3 +
          theme_epiestim() +
          lims(y = options_SI$ylim)
      }
    }
  }

  if (what == "incid") {
    p1 <- p1 +
      theme_epiestim() +
      theme(legend.position = "none") +
      scale_y_continuous(expand = c(0, 0))

    return(p1)
  }

  if (what == "R") {
    if (!legend) {
      p2 <- p2 +
        theme_epiestim() +
        theme(legend.position = "none")
    } else {
      p2 <- p2 +
        theme_epiestim() +
        theme(
          legend.background = element_blank(),
          legend.margin = margin(-0.8, 0, 0, 0, unit = "cm"),
          legend.key = element_rect(fill = NA)
        )
    }
    return(p2)
  }

  if (what == "SI") {
    p3 <- p3 +
      theme_epiestim() +
      theme(legend.position = "none") +
      scale_y_continuous(expand = c(0, 0))

    return(p3)
  }

  if (what == "all") {
    p1 <- p1 +
      theme_epiestim() +
      theme(legend.position = "none") +
      scale_y_continuous(expand = c(0, 0))

    p3 <- p3 +
      theme_epiestim() +
      theme(legend.position = "none") +
      scale_y_continuous(expand = c(0, 0))

    if (!legend) {
      p2 <- p2 +
        theme_epiestim() +
        theme(legend.position = "none")
    } else {
      p2 <- p2 +
        theme_epiestim() +
        theme(
          legend.background = element_blank(),
          legend.margin = margin(-0.8, 0, 0, 0, unit = "cm"),
          legend.key = element_rect(fill = NA)
        )
    }

    plot <- p1 + p2 + p3 + plot_layout(ncol = 1)

    return(plot)
  }
}
