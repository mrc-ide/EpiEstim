#' Check MCMC chain convergence using the Gelman-Rubin algorithm
#' 
#' This function splits an MCMC chain in two halves and uses the Gelman-Rubin 
#' algorithm to assess convergence of the chain by comparing its two halves.
#'
#' @param cdt_samples the `@sample` slot of a `cd.fit.mcmc` S4 object 
#' (see package `coarseDataTools`)
#' @return TRUE if the Gelman Rubin test for convergence was successful, FALSE 
#' otherwise
#' 
#' @seealso `estimate_R()`
#' @author Anne Cori
#' @export
#' @examples
#' \dontrun{
#' ## Note the following examples use an MCMC routine
#' ## to estimate the serial interval distribution from data,
#' ## so they may take a few minutes to run
#'
#' ## load data on rotavirus
#' data("MockRotavirus")
#'
#' ## estimate the serial interval from data
#' SI_fit <- coarseDataTools::dic.fit.mcmc(dat = MockRotavirus$si_data,
#'                      dist = "G",
#'                      init_pars = init_mcmc_params(MockRotavirus$si_data, "G"),
#'                      burnin = 1000,
#'                      n.samples = 5000)
#'
#' ## use check_cdt_samples_convergence to check convergence
#' converg_diag <- check_cdt_samples_convergence(SI_fit@samples)
#' converg_diag
#'
#' }

check_cdt_samples_convergence <- function(cdt_samples) {
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
            and decide whether to rerun for longer.")
    return(FALSE)
  } else {
    cat("\nGelman-Rubin MCMC convergence diagnostic was successful.")
    return(TRUE)
  }
}
