#########################################################
# Discretized serial interval (assuming a shifted gamma #
# distribution (with shift 1)                           #
#########################################################

#' Discretized Generation Time Distribution Assuming A Shifted Gamma 
#' Distribution
#'
#' \code{discr_si} computes the discrete distribution of the serial interval, 
#' assuming that the serial interval is shifted Gamma distributed, with shift 1.
#'
#' @param k Positive integer, or vector of positive integers for which the 
#' discrete distribution is desired.
#' @param mu A positive real giving the mean of the Gamma distribution.
#' @param sigma A non-negative real giving the standard deviation of the Gamma 
#' distribution.
#' @return Gives the discrete probability \eqn{w_k} that the serial interval is 
#' equal to \eqn{k}.
#' @details{
#' Assuming that the serial interval is shifted Gamma distributed with mean 
#' \eqn{\mu}, standard deviation \eqn{\sigma} and shift \eqn{1},
#' the discrete probability \eqn{w_k} that the serial interval is equal to 
#' \eqn{k} is:
#'
#' \deqn{w_k = kF_{\{\mu-1,\sigma\}}(k)+(k-2)F_{\{\mu-1,\sigma\}}
#' (k-2)-2(k-1)F_{\{\mu-1,\sigma\}}(k-1)\\
#' +(\mu-1)(2F_{\{\mu-1+\frac{\sigma^2}{\mu-1},
#' \sigma\sqrt{1+\frac{\sigma^2}{\mu-1}}\}}(k-1)-
#' F_{\{\mu-1+\frac{\sigma^2}{\mu-1},
#' \sigma\sqrt{1+\frac{\sigma^2}{\mu-1}}\}}(k-2)-
#' F_{\{\mu-1+\frac{\sigma^2}{\mu-1},
#' \sigma\sqrt{1+\frac{\sigma^2}{\mu-1}}\}}(k))}
#'
#' where \eqn{F_{\{\mu,\sigma\}}} is the cumulative density function of a Gamma 
#' distribution with mean \eqn{\mu} and standard deviation \eqn{\sigma}.
#' }
#' @seealso \code{\link{overall_infectivity}}, \code{\link{estimate_R}}
#' @author Anne Cori \email{a.cori@imperial.ac.uk}
#' @references Cori, A. et al. A new framework and software to estimate 
#' time-varying reproduction numbers during epidemics (AJE 2013).
#' @export
#' @examples
#' ## Computing the discrete serial interval of influenza
#' mean_flu_si <- 2.6
#' sd_flu_si <- 1.5
#' dicrete_si_distr <- discr_si(seq(0, 20), mean_flu_si, sd_flu_si)
#' plot(seq(0, 20), dicrete_si_distr, type = "h",
#'           lwd = 10, lend = 1, xlab = "time (days)", ylab = "frequency")
#' title(main = "Discrete distribution of the serial interval of influenza")
discr_si <- function(k, mu, sigma) 
{
  if (sigma < 0) {
    stop("sigma must be >=0.")
  }
  if (mu <= 1) {
    stop("mu must be >1")
  }
  if (any(k < 0)) {
    stop("all values in k must be >=0.")
  }

  a <- ((mu - 1) / sigma)^2
  b <- sigma^2 / (mu - 1)

  cdf_gamma <- function(k, a, b) stats::pgamma(k, shape = a, scale = b)

  res <- k * cdf_gamma(k, a, b) + 
    (k - 2) * cdf_gamma(k - 2, a, b) - 2 * (k - 1) * cdf_gamma(k - 1, a, b)
  res <- res + a * b * (2 * cdf_gamma(k - 1, a + 1, b) - 
                          cdf_gamma(k - 2, a + 1, b) - cdf_gamma(k, a + 1, b))
  res <- vnapply(res, function(e) max(0, e))

  return(res)
}
