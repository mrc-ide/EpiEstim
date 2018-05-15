EstimateRAmended2 <- function (I, T.Start, T.End, method = c("NonParametricSI", "ParametricSI",
                                                            "UncertainSI","NonParametricUncertainSI"), n1 = NULL, n2 = NULL, Mean.SI = NULL, Std.SI = NULL,
                              Std.Mean.SI = NULL, Min.Mean.SI = NULL, Max.Mean.SI = NULL,
                              Std.Std.SI = NULL, Min.Std.SI = NULL, Max.Std.SI = NULL,
                              SI.Distr = NULL, SI.Dist.Matrix = NULL, Mean.Prior = 5, Std.Prior = 5, CV.Posterior = 0.3,
                              plot = FALSE, leg.pos = "topright")
{
  CalculIncidencePerTimeStep <- function(I, T.Start, T.End) {
    NbTimePeriods <- length(T.Start)
    IncidencePerTimeStep <- sapply(1:NbTimePeriods, function(i) sum(I[T.Start[i]:T.End[i]]))
    return(IncidencePerTimeStep)
  }
  PosteriorFromSIDistr <- function(I, SI.Distr, a.Prior, b.Prior,
                                   T.Start, T.End) {
    NbTimePeriods <- length(T.Start)
    lambda <- OverallInfectivity(I, SI.Distr)
    FinalMean.SI <- sum(SI.Distr * (0:(length(SI.Distr) -
                                         1)))
    a.Posterior <- vector()
    b.Posterior <- vector()
    a.Posterior <- lapply(1:(NbTimePeriods), function(t) if (T.End[t] >
                                                             FinalMean.SI) {
      a.Prior + sum(I[T.Start[t]:T.End[t]])
    }
    else {
      NA
    })
    b.Posterior <- lapply(1:(NbTimePeriods), function(t) if (T.End[t] >
                                                             FinalMean.SI) {
      1/(1/b.Prior + sum(lambda[T.Start[t]:T.End[t]]))
    }
    else {
      NA
    })
    return(list(a.Posterior, b.Posterior))
  }

  SampleFromPosterior <- function(SampleSize, I, Mean.SI, Std.SI,
                                  a.Prior, b.Prior, T.Start, T.End) {
    NbTimePeriods <- length(T.Start)
    SI.Distr <- sapply(1:T, function(t) DiscrSI(t - 1, Mean.SI,
                                                Std.SI))
    FinalMean.SI <- sum(SI.Distr * (0:(length(SI.Distr) -
                                         1)))
    lambda <- OverallInfectivity(I, SI.Distr)
    a.Posterior <- vector()
    b.Posterior <- vector()
    a.Posterior <- sapply(1:(NbTimePeriods), function(t) if (T.End[t] >
                                                             FinalMean.SI) {
      a.Prior + sum(I[T.Start[t]:T.End[t]])
    }
    else {
      NA
    })
    b.Posterior <- sapply(1:(NbTimePeriods), function(t) if (T.End[t] >
                                                             FinalMean.SI) {
      1/(1/b.Prior + sum(lambda[T.Start[t]:T.End[t]], na.rm = TRUE))
    }
    else {
      NA
    })
    SampleR.Posterior <- sapply(1:(NbTimePeriods), function(t) if (!is.na(a.Posterior[t])) {
      rgamma(SampleSize, shape = unlist(a.Posterior[t]),
             scale = unlist(b.Posterior[t]))
    }
    else {
      rep(NA, SampleSize)
    })
    return(list(SampleR.Posterior, SI.Distr))
  }
   SampleFromPosterior2 <- function(SampleSize, I,  SI.Distr,
                                  a.Prior, b.Prior, T.Start, T.End) {
    NbTimePeriods <- length(T.Start)

    FinalMean.SI <- sum(SI.Distr * (0:(length(SI.Distr) -
                                         1)))
    lambda <- OverallInfectivity(I, SI.Distr)
    a.Posterior <- vector()
    b.Posterior <- vector()
    a.Posterior <- sapply(1:(NbTimePeriods), function(t) if (T.End[t] >
                                                             FinalMean.SI) {
      a.Prior + sum(I[T.Start[t]:T.End[t]])
    }
    else {
      NA
    })
    b.Posterior <- sapply(1:(NbTimePeriods), function(t) if (T.End[t] >
                                                             FinalMean.SI) {
      1/(1/b.Prior + sum(lambda[T.Start[t]:T.End[t]], na.rm = TRUE))
    }
    else {
      NA
    })
    SampleR.Posterior <- sapply(1:(NbTimePeriods), function(t) if (!is.na(a.Posterior[t])) {
      rgamma(SampleSize, shape = unlist(a.Posterior[t]),
             scale = unlist(b.Posterior[t]))
    }
    else {
      rep(NA, SampleSize)
    })
    return(list(SampleR.Posterior, SI.Distr))
  }
  method <- match.arg(method)
  if (is.vector(I) == FALSE) {
    stop("I must be a vector.")
  }
  T <- length(I)
  for (i in 1:T) {
    if (I[i] < 0) {
      stop("I must be a positive vector.")
    }
  }
  if (Mean.Prior <= 0) {
    stop("Mean.Prior must be >0.")
  }
  if (Std.Prior <= 0) {
    stop("Std.Prior must be >0.")
  }
  a.Prior <- (Mean.Prior/Std.Prior)^2
  b.Prior <- Std.Prior^2/Mean.Prior
  if (is.vector(T.Start) == FALSE) {
    stop("T.Start must be a vector.")
  }
  if (is.vector(T.End) == FALSE) {
    stop("T.End must be a vector.")
  }
  if (length(T.Start) != length(T.End)) {
    stop("T.Start and T.End must have the same length.")
  }
  NbTimePeriods <- length(T.Start)
  for (i in 1:NbTimePeriods) {
    if (T.Start[i] > T.End[i]) {
      stop("T.Start[i] must be <= T.End[i] for all i.")
    }
    if (T.Start[i] < 1 || T.Start[i]%%1 != 0) {
      stop("T.Start must be a vector of >0 integers.")
    }
    if (T.End[i] < 1 || T.End[i]%%1 != 0) {
      stop("T.End must be a vector of >0 integers.")
    }
  }
  if (method == "NonParametricSI") {
    if (is.null(SI.Distr) == TRUE) {
      stop("method NonParametricSI requires to specify the SI.Distr argument.")
    }
    if (is.vector(SI.Distr) == FALSE) {
      stop("method NonParametricSI requires that SI.Distr must be a vector.")
    }
    if (SI.Distr[1] != 0) {
      stop("method NonParametricSI requires that SI.Distr[1] = 0.")
    }
    if (length(SI.Distr) > 1) {
      for (i in 2:length(SI.Distr)) {
        if (SI.Distr[i] < 0) {
          stop("method NonParametricSI requires that SI.Distr must be a positive vector.")
        }
      }
    }
    if (abs(sum(SI.Distr) - 1) > 0.01) {
      stop("method NonParametricSI requires that SI.Distr must sum to 1.")
    }
  }
  if (method == "ParametricSI") {
    if (is.null(Mean.SI) == TRUE) {
      stop("method NonParametricSI requires to specify the Mean.SI argument.")
    }
    if (is.null(Std.SI) == TRUE) {
      stop("method NonParametricSI requires to specify the Std.SI argument.")
    }
    if (Mean.SI < 1) {
      stop("method ParametricSI requires a value >1 for Mean.SI.")
    }
    if (Std.SI < 0) {
      stop("method ParametricSI requires a >0 value for Std.SI.")
    }
  }
  if (method == "UncertainSI") {
    if (is.null(Mean.SI) == TRUE) {
      stop("method UncertainSI requires to specify the Mean.SI argument.")
    }
    if (is.null(Std.SI) == TRUE) {
      stop("method UncertainSI requires to specify the Std.SI argument.")
    }
    if (is.null(n1) == TRUE) {
      stop("method UncertainSI requires to specify the n1 argument.")
    }
    if (is.null(n2) == TRUE) {
      stop("method UncertainSI requires to specify the n2 argument.")
    }
    if (is.null(Std.Mean.SI) == TRUE) {
      stop("method UncertainSI requires to specify the Std.Mean.SI argument.")
    }
    if (is.null(Min.Mean.SI) == TRUE) {
      stop("method UncertainSI requires to specify the Min.Mean.SI argument.")
    }
    if (is.null(Max.Mean.SI) == TRUE) {
      stop("method UncertainSI requires to specify the Max.Mean.SI argument.")
    }
    if (is.null(Std.Std.SI) == TRUE) {
      stop("method UncertainSI requires to specify the Std.Std.SI argument.")
    }
    if (is.null(Min.Std.SI) == TRUE) {
      stop("method UncertainSI requires to specify the Min.Std.SI argument.")
    }
    if (is.null(Max.Std.SI) == TRUE) {
      stop("method UncertainSI requires to specify the Max.Std.SI argument.")
    }
    if (Mean.SI < 0) {
      stop("method UncertainSI requires a >0 value for Mean.SI.")
    }
    if (Std.SI < 0) {
      stop("method UncertainSI requires a >0 value for Std.SI.")
    }
    if (n2 <= 0 || n2%%1 != 0) {
      stop("method UncertainSI requires a >0 integer value for n2.")
    }
    if (n1 <= 0 || n1%%1 != 0) {
      stop("method UncertainSI requires a >0 integer value for n1.")
    }
    if (Std.Mean.SI < 0) {
      stop("method UncertainSI requires a >0 value for Std.Mean.SI.")
    }
    if (Min.Mean.SI < 1) {
      stop("method UncertainSI requires a value >=1 for Min.Mean.SI.")
    }
    if (Max.Mean.SI < Mean.SI) {
      stop("method UncertainSI requires that Max.Mean.SI >= Mean.SI.")
    }
    if (Mean.SI < Min.Mean.SI) {
      stop("method UncertainSI requires that Mean.SI >= Min.Mean.SI.")
    }
    if (signif(Max.Mean.SI - Mean.SI, 3) != signif(Mean.SI -
                                                   Min.Mean.SI, 3)) {
      warning("The distribution you chose for the mean SI is not centered around the mean.")
    }
    if (Std.Std.SI < 0) {
      stop("method UncertainSI requires a >0 value for Std.Std.SI.")
    }
    if (Min.Std.SI < 0) {
      stop("method UncertainSI requires a >0 value for Min.Std.SI.")
    }
    if (Max.Std.SI < Std.SI) {
      stop("method UncertainSI requires that Max.Std.SI >= Std.SI.")
    }
    if (Std.SI <= Min.Std.SI) {
      stop("method UncertainSI requires that Std.SI >= Min.Std.SI.")
    }
    if (signif(Max.Std.SI - Std.SI, 3) != signif(Std.SI -
                                                 Min.Std.SI, 3)) {
      warning("The distribution you chose for the std of the SI is not centered around the mean.")
    }
  }
  if (CV.Posterior < 0) {
    stop("CV.Posterior must be >0.")
  }
  MinNbCasesPerTimePeriod <- ceiling(1/CV.Posterior^2 - a.Prior)
  IncidencePerTimeStep <- CalculIncidencePerTimeStep(I, T.Start,
                                                     T.End)
  if (IncidencePerTimeStep[1] < MinNbCasesPerTimePeriod) {
    warning("You're estimating R too early in the epidemic to get the desired posterior CV.")
  }
  if (plot != TRUE && plot != FALSE) {
    stop("plot must be TRUE or FALSE.")
  }
  if (method == "NonParametricSI") {
    SIUncertainty <- "N"
    ParametricSI <- "N"
  }
  if (method == "ParametricSI") {
    SIUncertainty <- "N"
    ParametricSI <- "Y"
  }
  if (method == "UncertainSI") {
    SIUncertainty <- "Y"
  }
  if (method == "NonParametricUncertainSI") {
    ParametricSI <- "N"
    SIUncertainty <- "Y"
  }
  if (SIUncertainty == "Y") {
    if  (ParametricSI == "Y") {
    Mean.SI.sample <- rep(-1, n1)
    Std.SI.sample <- rep(-1, n1)
    for (k in 1:n1) {
      while (Mean.SI.sample[k] < Min.Mean.SI || Mean.SI.sample[k] >
             Max.Mean.SI) {
        Mean.SI.sample[k] <- rnorm(1, mean = Mean.SI,
                                   sd = Std.Mean.SI)
      }
      while (Std.SI.sample[k] < Min.Std.SI || Std.SI.sample[k] >
             Max.Std.SI || Std.SI.sample[k] > Mean.SI.sample[k]) {
        Std.SI.sample[k] <- rnorm(1, mean = Std.SI, sd = Std.Std.SI)
      }
    }
    temp <- lapply(1:n1, function(k) SampleFromPosterior(n2,
                                                         I, Mean.SI.sample[k], Std.SI.sample[k], a.Prior,
                                                         b.Prior, T.Start, T.End))
    SI.Distr <- cbind(t(sapply(1:n1, function(k) (temp[[k]])[[2]])),
                      rep(0, n1))
    Rsample <- matrix(NA, n2 * n1, NbTimePeriods)
    for (k in 1:n1) {
      Rsample[((k - 1) * n2 + 1):(k * n2), which(T.End >
                                                   Mean.SI.sample[k])] <- (temp[[k]])[[1]][, which(T.End >
                                                                                                     Mean.SI.sample[k])]
    }
    Mean.Posterior <- apply(Rsample, 2, mean, na.rm = TRUE)
    Std.Posterior <- apply(Rsample, 2, sd, na.rm = TRUE)
    Quantile.0.025.Posterior <- apply(Rsample, 2, quantile,
                                      0.025, na.rm = TRUE)
    Quantile.0.05.Posterior <- apply(Rsample, 2, quantile,
                                     0.05, na.rm = TRUE)
    Quantile.0.25.Posterior <- apply(Rsample, 2, quantile,
                                     0.25, na.rm = TRUE)
    Median.Posterior <- apply(Rsample, 2, median, na.rm = TRUE)
    Quantile.0.75.Posterior <- apply(Rsample, 2, quantile,
                                     0.75, na.rm = TRUE)
    Quantile.0.95.Posterior <- apply(Rsample, 2, quantile,
                                     0.95, na.rm = TRUE)
    Quantile.0.975.Posterior <- apply(Rsample, 2, quantile,
                                      0.975, na.rm = TRUE)
  }
  else {
    n1<-dim(SI.Dist.Matrix)[2]
    Mean.SI.sample <- rep(-1, n1)
    Std.SI.sample <- rep(-1, n1)
    for (k in 1:n1) {
        Mean.SI.sample[k] <- sum((1:dim(SI.Dist.Matrix)[1]-1)*SI.Dist.Matrix[,k])
        Std.SI.sample[k] <- sqrt(sum(SI.Dist.Matrix[,k]*((1:dim(SI.Dist.Matrix)[1]-1) - Mean.SI.sample[k])^2))
    }
    temp <- lapply(1:n1, function(k) SampleFromPosterior2(n2,
                                                         I, SI.Dist.Matrix[,k], a.Prior,
                                                         b.Prior, T.Start, T.End))
    SI.Distr <- cbind(t(sapply(1:n1, function(k) (temp[[k]])[[2]])),
                      rep(0, n1))
    Rsample <- matrix(NA, n2 * n1, NbTimePeriods)
    for (k in 1:n1) {
      Rsample[((k - 1) * n2 + 1):(k * n2), which(T.End >
                                                   Mean.SI.sample[k])] <- (temp[[k]])[[1]][, which(T.End >
                                                                                                     Mean.SI.sample[k])]
    }
    Mean.Posterior <- apply(Rsample, 2, mean, na.rm = TRUE)
    Std.Posterior <- apply(Rsample, 2, sd, na.rm = TRUE)
    Quantile.0.025.Posterior <- apply(Rsample, 2, quantile,
                                      0.025, na.rm = TRUE)
    Quantile.0.05.Posterior <- apply(Rsample, 2, quantile,
                                     0.05, na.rm = TRUE)
    Quantile.0.25.Posterior <- apply(Rsample, 2, quantile,
                                     0.25, na.rm = TRUE)
    Median.Posterior <- apply(Rsample, 2, median, na.rm = TRUE)
    Quantile.0.75.Posterior <- apply(Rsample, 2, quantile,
                                     0.75, na.rm = TRUE)
    Quantile.0.95.Posterior <- apply(Rsample, 2, quantile,
                                     0.95, na.rm = TRUE)
    Quantile.0.975.Posterior <- apply(Rsample, 2, quantile,
                                      0.975, na.rm = TRUE)
    }
  }


  else{
    # CertainSI
    if (ParametricSI == "Y") {
      SI.Distr <- sapply(1:T, function(t) DiscrSI(t - 1,
                                                  Mean.SI, Std.SI))
    }
    if (length(SI.Distr) < T + 1) {
      SI.Distr[(length(SI.Distr) + 1):(T + 1)] <- 0
    }
    FinalMean.SI <- sum(SI.Distr * (0:(length(SI.Distr) -
                                         1)))
    FinalStd.SI <- sqrt(sum(SI.Distr * (0:(length(SI.Distr) -
                                             1))^2) - FinalMean.SI^2)
    post <- PosteriorFromSIDistr(I, SI.Distr, a.Prior, b.Prior,
                                 T.Start, T.End)
    a.Posterior <- unlist(post[[1]])
    b.Posterior <- unlist(post[[2]])
    Mean.Posterior <- a.Posterior * b.Posterior
    Std.Posterior <- sqrt(a.Posterior) * b.Posterior
    Quantile.0.025.Posterior <- qgamma(0.025, shape = a.Posterior,
                                       scale = b.Posterior, lower.tail = TRUE, log.p = FALSE)
    Quantile.0.05.Posterior <- qgamma(0.05, shape = a.Posterior,
                                      scale = b.Posterior, lower.tail = TRUE, log.p = FALSE)
    Quantile.0.25.Posterior <- qgamma(0.25, shape = a.Posterior,
                                      scale = b.Posterior, lower.tail = TRUE, log.p = FALSE)
    Median.Posterior <- qgamma(0.5, shape = a.Posterior,
                               scale = b.Posterior, lower.tail = TRUE, log.p = FALSE)
    Quantile.0.75.Posterior <- qgamma(0.75, shape = a.Posterior,
                                      scale = b.Posterior, lower.tail = TRUE, log.p = FALSE)
    Quantile.0.95.Posterior <- qgamma(0.95, shape = a.Posterior,
                                      scale = b.Posterior, lower.tail = TRUE, log.p = FALSE)
    Quantile.0.975.Posterior <- qgamma(0.975, shape = a.Posterior,
                                       scale = b.Posterior, lower.tail = TRUE, log.p = FALSE)
  }
  results <- list()
  results$R <- as.data.frame(cbind(T.Start, T.End, Mean.Posterior,
                                   Std.Posterior, Quantile.0.025.Posterior, Quantile.0.05.Posterior,
                                   Quantile.0.25.Posterior, Median.Posterior, Quantile.0.75.Posterior,
                                   Quantile.0.95.Posterior, Quantile.0.975.Posterior))
  names(results$R) <- c("T.Start", "T.End", "Mean(R)", "Std(R)",
                        "Quantile.0.025(R)", "Quantile.0.05(R)", "Quantile.0.25(R)",
                        "Median(R)", "Quantile.0.75(R)", "Quantile.0.95(R)",
                        "Quantile.0.975(R)")
  if (SIUncertainty == "Y") {
    results$SIDistr <- as.data.frame(cbind(Mean.SI.sample,
                                           Std.SI.sample))
  }
  else {
    if (ParametricSI == "Y") {
      if (length(which(abs(cumsum(SI.Distr) - 1) < 0.01)) ==
          0) {
        warning("The serial interval distribution you have chosen is very wide compared to the duration of the epidemic.\nEstimation will be performed anyway but restults should be interpreted with care.")
        MaxT <- length(cumsum(SI.Distr))
      }
      else {
        MaxT <- min(which(abs(cumsum(SI.Distr) - 1) <
                            0.01))
      }
      results$SIDistr <- as.data.frame(cbind(0:(MaxT -
                                                  1), SI.Distr[1:MaxT]))
      names(results$SIDistr) <- c("k", "w[k]")
    }
    else {
      results$SIDistr <- as.data.frame(cbind(FinalMean.SI,
                                             FinalStd.SI))
      names(results$SIDistr) <- c("Mean Discrete SI", "Std Discrete SI")
    }
  }
  if (plot == TRUE) {
    grey <- "#999999"
    if (SIUncertainty == "Y") {
      par(mfcol = c(2, 2), las = 1, cex.main = 1.5, cex.lab = 1.2,
          cex.axis = 1, mar = c(4.8, 4.8, 2.4, 0.8), mgp = c(4,
                                                             1, 0))
      plot(I, type = "s", bty = "n", xlab = "", ylab = "",
           main = "Epidemic curve")
      title(xlab = "Time", ylab = "Incidence", line = 3)
      plot(T.End, Median.Posterior, type = "l", bty = "n",
           xlab = "", ylab = "", main = "Estimated R", ylim = c(0,
                                                                max(Quantile.0.975.Posterior, na.rm = TRUE)),
           xlim = c(1, T))
      title(xlab = "Time", ylab = "R", line = 3)
      polygon(c(T.End, rev(T.End)), c(Quantile.0.025.Posterior,
                                      rev(Quantile.0.975.Posterior)), col = grey, border = FALSE)
      lines(T.End, Median.Posterior)
      lines(0:T, rep(1, T + 1), lty = 2)
      legend(leg.pos, c("Median", "95%CrI"), col = c("Black",
                                                     grey), lwd = c(1, 10), bty = "n", cex = 1.2)
      hist(Mean.SI.sample, xlab = "", ylab = "", main = "Explored \n mean serial intervals",
           freq = FALSE)
      title(xlab = "Mean serial interval", ylab = "Density",
            line = 3)
      hist(Std.SI.sample, xlab = "", ylab = "", main = "Explored \n std serial intervals",
           freq = FALSE)
      title(xlab = "Std serial interval", ylab = "Density",
            line = 3)
    }
    else {
      par(mfrow = c(3, 1), las = 1, cex.main = 1.8, cex.lab = 1.5,
          cex.axis = 1.2, mar = c(6, 6, 3, 1), mgp = c(4,
                                                       1, 0))
      plot(I, type = "s", bty = "n", xlab = "Time", ylab = "Incidence",
           main = "Epidemic curve")
      plot(T.End, Median.Posterior, type = "l", bty = "n",
           xlab = "Time", ylab = "R", main = "Estimated R",
           ylim = c(0, max(Quantile.0.975.Posterior, na.rm = TRUE)),
           xlim = c(1, T))
      polygon(c(T.End, rev(T.End)), c(Quantile.0.025.Posterior,
                                      rev(Quantile.0.975.Posterior)), col = grey, border = FALSE)
      lines(T.End, Median.Posterior)
      lines(0:T, rep(1, T + 1), lty = 2)
      legend(leg.pos, c("Median", "95%CrI"), col = c("Black",
                                                     grey), lwd = c(1, 10), bty = "n", cex = 1.2)
      plot(0:(length(SI.Distr) - 1), SI.Distr, type = "h",
           lwd = 10, lend = 1, bty = "n", xlab = "Time",
           ylab = "Frequency", main = "Serial interval distribution",
           xlim = c(0, FinalMean.SI + 6 * FinalStd.SI))
    }
  }
  return(results)
}
