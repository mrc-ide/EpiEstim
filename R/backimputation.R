#' Impute unobserved generations of infection
#' 
#' This function imputes incidence prior the first date of reported cases to
#' tackle bias in R_0 estimates. 
#' It does so by fitting a simple linear model on the logged-incidence cases,
#' based on observation in an initial window.
#' No cases are assumed to be imported.
#'
#' @param incid the raw, reported incidence cases.
#' @param window length of the observation window to fit the exponential growth
#' model for back-imputation
#'
#' @return an incidence data.frame, combining back-imputed cases for 100 time 
#' points(with rows indexed by a negative integer rowname) and cases (with rows
#' indexed by a non-negative integer) 
#' @export
#'
#' @examples
#' incid_all <- ceiling( exp(.3 * 0:20) )
#' incid_trunc <- tail(incid, 10)
#' x <- backimpute_I(incid=incid_trunc, window=6)
#' idx <- as.integer(rownames(x)) > -10
#' x[idx, ]$local 
#' incid_all
#'

backimpute_I <- function(incid, window) {

    stopifnot("Backimputation window needs to have integer length"=! is.integer(window))

    # process observed incidence
    incid_processed <- process_I(incid)
    incid_processed[1, ] <- c(sum(incid_processed[1, ]), 0)

    # some cases may be 0, implying -infinite logs
    safe_shift <- .5

    # backimpute unobserved, previous cases based on first window of observations
    log_incid_start <- data.frame(
        t = seq(window),
        li = log(incid[1:window] + safe_shift )
    )
    imputed_t <- seq(from = -100, to = 0)
    fit_backimpute <- lm(li ~ t, data = log_incid_start)

    predict_backimpute_log <- predict.lm(fit_backimpute, newdata = list(t = imputed_t))
    predict_backimpute <- data.frame(
        local = exp(predict_backimpute_log) - safe_shift,
        imported = rep(0, length(imputed_t))
    )
    rownames(predict_backimpute) <- imputed_t

    # exclude negative cases arising from shift before logs.
    idx_nonnegative <- predict_backimpute$local >= 0
    predict_backimpute <- predict_backimpute[ idx_nonnegative ]

    # first observation must be imported, otherwise later warning
    predict_backimpute[1, ] <- c(0, sum(predict_backimpute[1, ]))

    # bind imputed with observed
    incid_with_imputed <- rbind(predict_backimpute, incid_processed)
    return(incid_with_imputed)
}
