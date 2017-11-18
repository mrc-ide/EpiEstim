#######################################################################################################################
# check_cdtsamples_convergence runs a Gelman Rubin test to check convergence of the MCMC chain in coarseDataTools #
#######################################################################################################################

#' check_cdtsamples_convergence TITLE
#' 
#' \code{check_cdtsamples_convergence} DESCRIPTION TO COME.
#' 
#' @param cdtsamples XXXXXXX.
#' @return TRUE if the Gelman Rubin test for convergence was successful, FALSE otherwise
#' @details{
#' XXXXXXX. 
#' }
#' @seealso XXXXXXX.
#' @author XXXXXXX.
#' @references XXXXXXX.
#' @importFrom coda gelman.diag
#' @export
#' @examples
#' ## XXXXXXX.
check_cdtsamples_convergence <- function(cdtsamples)
{
  
  #############################################################################################################################
  # checking convergence of the MCMC by using the Gelman-Rubin algorithm between the first and second half of the MCMC sample
  spl1 <- cdtsamples[1:floor(nrow(cdtsamples)/2),]
  spl2 <- cdtsamples[(ceiling(nrow(cdtsamples)/2)+1):nrow(cdtsamples),]
  GRD <- gelman.diag(as.mcmc.list(list(as.mcmc(spl1), as.mcmc(spl2))))
  # Is any of the potential scale reduction factors >1.1 (looking at the upper CI)? 
  # If so this would suggest that the MCMC has not converged well. 
  if(any(GRD$psrf[,"Upper C.I."]>1.1))
  {
    warning("The Gelman-Rubin algorithm suggests the MCMC may not have converged within the number of iterations (MCMC.burnin + n1) specified. 
            You can visualise the full MCMC chain using: \n
            > par(mfrow=c(2,1))
            > plot(res$SI.Moments[,'Mean'], type='l', xlab='Iterations', ylab='Mean SI') 
            > plot(res$SI.Moments[,'Std'], type='l', xlab='Iterations', ylab='Std SI'),
            where res is the output of EstimateR
            and decide whether to rerun for longer.")
    return(FALSE)
  }else
  {
    cat("\nGelman-Rubin MCMC convergence diagnostic was successful.")
    return(TRUE)
  }
  
}
