#' Check MCMC chain convergence using the Gelman-Rubin algorithm
#'
#' This function is deprecated as [estimate_R()] estimates the serial interval
#' by maximum likelihood. See [primary2estim()] for using primarycensored fits.
#'
#' Splits an MCMC chain in two halves and uses the Gelman-Rubin algorithm to
#' assess convergence of the chain by comparing its two halves.
#'
#' @param cdt_samples the `@samples` slot of the output of
#'   `coarseDataTools::dic.fit.mcmc()`
#' @return TRUE if the Gelman Rubin test for convergence was successful, FALSE
#' otherwise
#'
#' @seealso [primary2estim()]
#' @author Anne Cori
#' @export
check_cdt_samples_convergence <- function(cdt_samples) {
  .Deprecated(msg = paste(
    "check_cdt_samples_convergence() is deprecated as estimate_R()",
    "estimates the serial interval by maximum likelihood. See",
    "primary2estim() for using primarycensored fits."
  ))
  ## checking convergence of the MCMC by using the Gelman-Rubin algorithm 
  ## between the first and second half of the MCMC sample
  spl1 <- cdt_samples[seq_len(floor(nrow(cdt_samples) / 2)), ]
  spl2 <- cdt_samples[seq(ceiling(nrow(cdt_samples) / 2) + 1, nrow(cdt_samples)), ]
  GRD <- coda::gelman.diag(coda::as.mcmc.list(list(coda::as.mcmc(spl1), coda::as.mcmc(spl2))))
  # Is any of the potential scale reduction factors >1.1 
  # (looking at the upper CI)?
  # If so this would suggest that the MCMC has not converged well.
  if (any(GRD$psrf[, "Upper C.I."] > 1.1)) {
    warning("The Gelman-Rubin algorithm suggests the MCMC may not have converged
within the number of iterations (MCMC.burnin + n1) specified.
            You can visualise the full MCMC chain using: \n
            > par(mfrow=c(2,1))
            > plot(res$SI.Moments[,'Mean'], type='l', xlab='Iterations', 
ylab='Mean SI')
            > plot(res$SI.Moments[,'Std'], type='l', xlab='Iterations', 
ylab='Std SI'),
            where res is the output of estimate_R
            and decide whether to rerun for longer.", call. = FALSE)
    return(FALSE)
  } else {
    cat("\nGelman-Rubin MCMC convergence diagnostic was successful.")
    return(TRUE)
  }
}
