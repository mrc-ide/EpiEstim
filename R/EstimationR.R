#############################################
#############################################
# Functions                                 #
#############################################
#############################################

#########################################################
# Discretized serial interval (assuming a shifted gamma #
# distribution (with shift 1)                           #
#########################################################

DiscrSI<-function(k,mu,sigma)
{
	if(sigma<0)
	{
		stop("sigma must be >=0.")
	}
	a=((mu-1)/sigma)^2
	b=sigma^2/(mu-1)
	CDFGamma<-function(k,a,b)
	{
		return(pgamma(k,shape=a,scale=b))
	}
	res<-k*CDFGamma(k,a,b)+(k-2)*CDFGamma(k-2,a,b)-2*(k-1)*CDFGamma(k-1,a,b)
	res<-res+a*b*(2*CDFGamma(k-1,a+1,b)-CDFGamma(k-2,a+1,b)-CDFGamma(k,a+1,b))
	res<-max(0,res)
	return(res)
}

#########################################################
# Calculates Lambda_t = Sum_1^t I_{t-s} * w_s           #
# with I incidence and w discrete SI distribution       #
#########################################################

OverallInfectivity <-function (I,SI.Distr)
{
	if(is.vector(I)==FALSE)
	{
		stop("Incidence must be a vector.")
	}
	T<-length(I)
	for(i in 1:T)
	{
		if(I[i]<0)
		{
			stop("Incidence must be a positive vector.")
		}
	}
	if(is.vector(SI.Distr)==FALSE)
	{
		stop("SI.Distr must be a vector.")
	}
	if(SI.Distr[1]!=0)
	{
		stop("SI.Distr[1] needs to be 0.")
	}
	if(length(SI.Distr)>1)
	{
		for(i in 2:length(SI.Distr))
		{
			if(SI.Distr[i]<0)
			{
				stop("SI.Distr must be a positive vector.")
			}
		}
	}
	if(abs(sum(SI.Distr)-1)>0.01)
	{
		stop("SI.Distr must sum to 1.")
	}
	lambda <- vector()
	lambda[1]<-NA
	for (t in 2:length(I))
	{
		lambda[t] <- sum(SI.Distr[1:t]*I[t:1],na.rm=TRUE)
	}
	return(lambda)
}

#########################################################
# EstimateR_func: Doing the heavy work in EstimateR     #
#########################################################

EstimateR_func <- function (I, T.Start, T.End, method = c("NonParametricSI", "ParametricSI",
                                                         "UncertainSI","NonParametricUncertainSI"), n1 = NULL, n2 = NULL, Mean.SI = NULL, Std.SI = NULL,
                           Std.Mean.SI = NULL, Min.Mean.SI = NULL, Max.Mean.SI = NULL,
                           Std.Std.SI = NULL, Min.Std.SI = NULL, Max.Std.SI = NULL,
                           SI.Distr = NULL, SI.Dist.Matrix = NULL, Mean.Prior = 5, Std.Prior = 5, CV.Posterior = 0.3,
                           plot = FALSE, leg.pos = "topright")
{
 
  #########################################################
  # Calculates the cumulative incidence over time steps   #
  #########################################################
  
   CalculIncidencePerTimeStep <- function(I, T.Start, T.End) {
    NbTimePeriods <- length(T.Start)
    IncidencePerTimeStep <- sapply(1:NbTimePeriods, function(i) sum(I[T.Start[i]:T.End[i]]))
    return(IncidencePerTimeStep)
   }
   
   #########################################################
   # Calculates the parameters of the Gamma posterior      #
   # distribution from the discrete SI distribution        #
   #########################################################
   
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

  #########################################################
  # Samples from the Gamma posterior distribution for a   #
  # given mean SI and std SI                              #
  #########################################################
  
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
    ParametricSI <- "Y"
  }
  if (method == "NonParametricUncertainSI") {
    SIUncertainty <- "Y"
    ParametricSI <- "N"
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

#######################################################################################################################
# coarse2estim integates CoarseDataTools with EpiEstim using the amended version of EstimateR called EstimateR_func #
#######################################################################################################################

coarse2estim <- function(object, n_samples=1000){

  samples0 <- as.matrix(object@samples)
  index <- sample(1:nrow(samples0), size= n_samples)
  samples <- samples0[index, ]
  dist <- object@dist

  ##  Probability matrix that will be used in EpiEstim based on which distribution is specified by the user
  if (dist == "G"){
    max_interval <- c(10^(-4), 1:ceiling(qgamma(0.999, shape=object@ests[1,1], scale=object@ests[2,1])))
    prob_matrix <- apply(samples, 1, function(x) dgamma(max_interval, shape=x[1], scale=x[2]))

  } else if (dist == "off1G"){
    # offset gamma distribution with shifted min and max value of max serial interval
    max_interval <- c(10^(-4), 1:ceiling(qgamma(0.999, shape=object@ests[1,1], scale=object@ests[2,1])))
    prob_matrix <- apply(samples, 1, function(x) dgamma(max_interval, shape=x[1], scale=x[2]))
  }
  else if (dist == "W"){
    max_interval <- c(10^(-4), 1:ceiling(qgamma(0.999, shape=object@ests[1,1], scale=object@ests[2,1])))
    prob_matrix <- apply(samples, 1, function(x) dweibull(max_interval, shape=x[1], scale=x[2]))

  } else if (dist == "L"){
    max_interval <- c(10^(-4), 1:ceiling(qgamma(0.999, shape=object@ests[1,1], scale=object@ests[2,1])))
    prob_matrix <- apply(samples, 1, function(x) dlnorm(max_interval, meanlog=x[1], sdlog=x[2]))
  } else {
    stop(sprintf("Distribtion (%s) not supported",dist))
  }
  prob_matrix <- apply(prob_matrix, 2, function(x) x/sum(x))

  prob_matrix <- rbind(rep(0, ncol(prob_matrix)), prob_matrix)
  out <- list(prob_matrix = prob_matrix, dist = dist)

  return(out)
}


#########################################################################################################################
# EstimareR is a wraper which replace the old EstimateR with EstimateR_func and accepts an object of class "cd.fit.mcmc"#
#########################################################################################################################

EstimateR <- function(I, T.Start, T.End, method = c("NonParametricSI", "ParametricSI",
                                                    "UncertainSI","NonParametricUncertainSI"), 
                      n1 = NULL, n2 = NULL, Mean.SI = NULL, Std.SI = NULL,
                      Std.Mean.SI = NULL, Min.Mean.SI = NULL, Max.Mean.SI = NULL,
                      Std.Std.SI = NULL, Min.Std.SI = NULL, Max.Std.SI = NULL,
                      SI.Distr = NULL, Mean.Prior = 5, Std.Prior = 5, CV.Posterior = 0.3,
                      plot = FALSE, leg.pos = "topright", CDT = NULL) {
  
  ### Need to add warnings if method="NonParametricUncertainSI" and CDT is not null 
  
  if (!is.null(CDT)) {
    # Warning if the CDT object is not of the S4 class "cd.fit.mcmc"
    if (class(CDT)[1]!="cd.fit.mcmc")
      warning("CDT needs to be defined as an object of the S4 class 'cd.fit.mcmc??")
    
    c2e <- coarse2estim(CDT)
    
    EstimateR_func(I=I, T.Start=T.Start, T.End=T.End, method = "NonParametricUncertainSI", n1=n1 , n2=n2 , Mean.SI=Mean.SI , Std.SI=Std.SI ,
                   Std.Mean.SI=Std.Mean.SI , Min.Mean.SI=Min.Mean.SI , Max.Mean.SI=Max.Mean.SI ,
                   Std.Std.SI=Std.Std.SI , Min.Std.SI=Min.Std.SI , Max.Std.SI=Max.Std.SI ,
                   SI.Distr=SI.Distr , SI.Dist.Matrix= c2e$prob_matrix , Mean.Prior=Mean.Prior , Std.Prior=Std.Prior, CV.Posterior=CV.Posterior ,
                   plot=plot , leg.pos=leg.pos)
    
  } else {
    EstimateR_func(I=I, T.Start=T.Start, T.End=T.End, method = method, n1=n1 , n2=n2 , Mean.SI=Mean.SI , Std.SI=Std.SI ,
                   Std.Mean.SI=Std.Mean.SI , Min.Mean.SI=Min.Mean.SI , Max.Mean.SI=Max.Mean.SI ,
                   Std.Std.SI=Std.Std.SI , Min.Std.SI=Min.Std.SI , Max.Std.SI=Max.Std.SI ,
                   SI.Distr=SI.Distr , SI.Dist.Matrix= NULL, Mean.Prior=Mean.Prior , Std.Prior=Std.Prior, CV.Posterior=CV.Posterior ,
                   plot=plot , leg.pos=leg.pos)
  }
}

#########################################################################################################################
# WT function to estimate Rc the case reproduction number #
#########################################################################################################################

WT <- function(I,T.Start,T.End,method=c("NonParametricSI","ParametricSI"),Mean.SI=NULL,Std.SI=NULL,SI.Distr=NULL,nSim=10,plot=FALSE,leg.pos="topright")
{
    
    ### Functions ###
    
    #########################################################
    # Draws a possile transmission tree                     #
    #########################################################
    
    DrawOneSetOfAncestries <- function()
    {
        res <- vector()
        for(t in T.Start[1]:T.End[length(T.End)])
        {
            if(length(which(Onset==t))>0)
            {
                if(length(possibleAncesTime[[t]])>0)
                {
                	prob <- SI.Distr[t-possibleAncesTime[[t]]+1]*I[possibleAncesTime[[t]]]
                	if(any(prob>0))
                	{
                    		res[which(Onset==t)] <- possibleAncesTime[[t]][which(rmultinom(length(which(Onset==t)),size=1,prob=prob)==TRUE,arr.ind=TRUE)[,1]]
                    }else
                    {
                    		res[which(Onset==t)] <- NA
                    }
                }else
                res[which(Onset==t)] <- NA
            }
        }
        return(res)
    }
    
    ### Adjusting T.Start and T.End so that at least an incident case has been observed before T.Start[1] ###
    
    i <- 1
    while(sum(I[1:(T.Start[i]-1)])==0)
    {
        i <- i+1
    }
    temp <- which(T.Start<i)
    if(length(temp>0))
    {
        T.Start <- T.Start[-temp]
        T.End <- T.End[-temp]
    }
    
    ### Error messages ###
    
    method <- match.arg(method)
    
    if(is.vector(I)==FALSE)
    {
        stop("I must be a vector.")
    }
    T<-length(I)
    for(i in 1:T)
    {
        if(I[i]<0)
        {
            stop("I must be a positive vector.")
        }
    }
    I[which(is.na(I))] <- 0
    
    if(is.vector(T.Start)==FALSE)
    {
        stop("T.Start must be a vector.")
    }
    if(is.vector(T.End)==FALSE)
    {
        stop("T.End must be a vector.")
    }
    if(length(T.Start)!=length(T.End))
    {
        stop("T.Start and T.End must have the same length.")
    }
    NbTimePeriods<-length(T.Start)
    for(i in 1:NbTimePeriods)
    {
        if(T.Start[i]>T.End[i])
        {
            stop("T.Start[i] must be <= T.End[i] for all i.")
        }
        if(T.Start[i]<1 || T.Start[i]%%1!=0)
        {
            stop("T.Start must be a vector of >0 integers.")
        }
        if(T.End[i]<1 || T.End[i]%%1!=0)
        {
            stop("T.End must be a vector of >0 integers.")
        }
    }
    
    if(method=="NonParametricSI")
    {
        if(is.null(SI.Distr)==TRUE)
        {
            stop("method NonParametricSI requires to specify the SI.Distr argument.")
        }
        if(is.vector(SI.Distr)==FALSE)
        {
            stop("method NonParametricSI requires that SI.Distr must be a vector.")
        }
        if(SI.Distr[1]!=0)
        {
            stop("method NonParametricSI requires that SI.Distr[1] = 0.")
        }
        if(length(SI.Distr)>1)
        {
            for(i in 2:length(SI.Distr))
            {
                if(SI.Distr[i]<0)
                {
                    stop("method NonParametricSI requires that SI.Distr must be a positive vector.")
                }
            }
        }
        if(abs(sum(SI.Distr)-1)>0.01)
        {
            stop("method NonParametricSI requires that SI.Distr must sum to 1.")
        }
        SI.Distr <- c(SI.Distr,0)
    }
    
    if(method=="ParametricSI")
    {
        if(is.null(Mean.SI)==TRUE)
        {
            stop("method NonParametricSI requires to specify the Mean.SI argument.")
        }
        if(is.null(Std.SI)==TRUE)
        {
            stop("method NonParametricSI requires to specify the Std.SI argument.")
        }
        if(Mean.SI<1)
        {
            stop("method ParametricSI requires a value >1 for Mean.SI.")
        }
        if(Std.SI<0)
        {
            stop("method ParametricSI requires a >0 value for Std.SI.")
        }
    }
    
    if(plot!=TRUE && plot!=FALSE)
    {
        stop("plot must be TRUE or FALSE.")
    }
    
    if(is.numeric(nSim)==FALSE)
    {
        stop("nSim must be a positive integer.")
    }
    if(nSim<=0)
    {
        stop("nSim must be a positive integer.")
    }
    
    ### What does each method do ###
    
    if(method=="NonParametricSI")
    {
        ParametricSI<-"N"
    }
    if(method=="ParametricSI")
    {
        ParametricSI<-"Y"
    }
    
    if(ParametricSI=="Y")
    {
        SI.Distr <- sapply(1:T, function(t) DiscrSI(t-1,Mean.SI,Std.SI))
    }
    if(length(SI.Distr)<T+1){SI.Distr[(length(SI.Distr)+1):(T+1)]<-0}
    FinalMean.SI<-sum(SI.Distr*(0:(length(SI.Distr)-1)))
    FinalStd.SI<-sqrt(sum(SI.Distr*(0:(length(SI.Distr)-1))^2)-FinalMean.SI^2)
    
    TimePeriodsWithNoIncidence <- vector()
    for(i in 1:NbTimePeriods)
    {
        if(sum(I[T.Start[i]:T.End[i]])==0)
        {
            TimePeriodsWithNoIncidence <- c(TimePeriodsWithNoIncidence,i)
        }
    }
    if(length(TimePeriodsWithNoIncidence)>0)
    {
        T.Start <- T.Start[-TimePeriodsWithNoIncidence]
        T.End <- T.End[-TimePeriodsWithNoIncidence]
        NbTimePeriods <- length(T.Start)
    }
    
    Onset <- vector()
    for (t in 1:T) {
        Onset <- c(Onset, rep(t, I[t]))
    }
    NbCases <- length(Onset)
    
    delay <- outer (1:T, 1:T, "-")
    SIdelay <- apply(delay, 2, function(x) SI.Distr[pmin(pmax(x + 1, 1), length(SI.Distr))])
    SumOnColSIDelay_tmp <- sapply(1:nrow(SIdelay), function(i) sum(SIdelay [i,]* I, na.rm=TRUE) )
    SumOnColSIDelay <- vector()
    for (t in 1:T) {
        SumOnColSIDelay <- c(SumOnColSIDelay, rep(SumOnColSIDelay_tmp[t], I[t]))
    }
    MatSumOnColSIDelay <- matrix(rep(SumOnColSIDelay_tmp, T), nrow = T, ncol = T)
    p <- SIdelay/(MatSumOnColSIDelay)
    p[which(is.na(p))] <- 0
    p[which(is.infinite(p))] <- 0
    MeanRperIndexCaseDate <- sapply(1:ncol(p), function(j) sum(p[,j]*I, na.rm=TRUE))
    MeanRperDate.WT <- sapply(1:NbTimePeriods, function(i) mean(rep(MeanRperIndexCaseDate[which((1:T >= T.Start[i]) * (1:T <= T.End[i]) == 1)], I[which((1:T >= T.Start[i]) * (1:T <= T.End[i]) == 1)]) ) )
    
    possibleAncesTime <- sapply(1:T,function(t) (t-(which(SI.Distr!=0))+1)[which(t-(which(SI.Distr!=0))+1>0)])
    ancestriesTime <- t(sapply(1:nSim , function(i) DrawOneSetOfAncestries()))
    
    Rsim <- sapply(1:NbTimePeriods,function(i) rowSums((ancestriesTime[,]>=T.Start[i]) * (ancestriesTime[,]<=T.End[i]),na.rm=TRUE)/sum(I[T.Start[i]:T.End[i]]))
    
    R025.WT <- apply(Rsim, 2, quantile,0.025,na.rm=TRUE)
    R025.WT <- R025.WT[which(!is.na(R025.WT))]
    R975.WT <- apply(Rsim, 2, quantile,0.975,na.rm=TRUE)
    R975.WT <- R975.WT[which(!is.na(R975.WT))]
    std.WT <- apply(Rsim, 2, sd,na.rm=TRUE)
    std.WT <- std.WT[which(!is.na(std.WT))]
    
    results<-list()
    
    results$R<-as.data.frame(cbind(T.Start,T.End,MeanRperDate.WT,std.WT,R025.WT,R975.WT))
    names(results$R)<-c("T.Start","T.End","Mean(R)","Std(R)","Quantile.0.025(R)","Quantile.0.975(R)")
    
    if(ParametricSI=="Y") # method "ParametricSI"
    {
        if(length(which(abs(cumsum(SI.Distr)-1)<0.01))==0)
        {
            warning("The serial interval distribution you have chosen is very wide compared to the duration of the epidemic.\nEstimation will be performed anyway but restults should be interpreted with care.")
            MaxT <- length(cumsum(SI.Distr))
        }else
        {
            MaxT<-min(which(abs(cumsum(SI.Distr)-1)<0.01))
        }
        results$SIDistr<-as.data.frame(cbind(0:(MaxT-1),SI.Distr[1:MaxT]))
        names(results$SIDistr)<-c("k","w[k]")
    }else # method "NonParametricSI"
    {
        results$SIDistr<-as.data.frame(cbind(FinalMean.SI,FinalStd.SI))
        names(results$SIDistr)<-c("Mean Discrete SI","Std Discrete SI")
    }
    
    if(plot==TRUE)
    {
        grey <- "#999999"
        
        par(mfrow=c(3,1),las=1,cex.main=1.8,cex.lab=1.5,cex.axis=1.2,mar=c(6,6,3,1),mgp=c(4,1,0))
        plot(I,type="s",bty="n",xlab="Time",ylab="Incidence",main="Epidemic curve")		
        plot(T.End,MeanRperDate.WT,type="p",bty="n",xlab="Time",ylab=expression(R^c),main=expression(paste("Estimated ",R^c,sep="")),ylim=c(0,max(R975.WT,na.rm=TRUE)),xlim=c(1,T),pch=20)
        for(i in 1:length(T.End))
        {
            segments(T.End[i],R025.WT[i],T.End[i],R975.WT[i],col=grey)
        }
        points(T.End,MeanRperDate.WT,pch=20)
        lines(0:T,rep(1,T+1),lty=2)
        legend(leg.pos,c("Mean","95% CI"),col=c("Black",grey),lwd=c(0,1.5),pch=c(19,19),pt.cex=c(1,0),bty="n",cex=1.2)
        
        plot(0:(length(SI.Distr)-1),SI.Distr,type="h",lwd=10,lend=1,bty="n",xlab="Time",ylab="Frequency",main="Serial interval distribution",xlim=c(0,FinalMean.SI+6*FinalStd.SI))
        
    }
    
    return(results)
    
}
