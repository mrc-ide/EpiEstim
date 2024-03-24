#' Impute unobserved generations of infection
#' 
#' This function imputes incidence prior the first date of reported cases to
#' address early bias in R_t estimates. 
#' A simple linear model is fitted on shifted, logged-incidence cases, based on 
#' an initial observation window.
#' Currently, no cases are assumed to be imported.
#'
#' @param incid the raw, reported incidence cases.
#' @param window_b length of the observation window to fit the exponential 
#' growth model for back-imputation
#'
#' @return an incidence data.frame, combining back-imputed cases for 100 time 
#' points(with rows indexed by a negative integer rowname) and cases (with rows
#' indexed by a non-negative integer) 
#' @export
#'
#' @examples
#' incid_all <- ceiling( exp(.3 * 0:20) )
#' incid_trunc <- tail(incid_all, 10)
#' x <- backimpute_I(incid=incid_trunc, window_b=6)
#' idx <- as.integer(rownames(x)) > -10
#' x[idx, ]$local 
#' incid_all
#'

backimpute_I <- function(incid, window_b) {

    if( inherits(incid, 'incidence') ){
      msg <- "incidence objects are currently not supported by backimpute_I()."
      stop(msg)
    }
  
    stopifnot("Backimputation window needs to contain at least 2 timepoints"=
        window_b >= 2)
    stopifnot("Backimputation window needs to have integer length"= 
        window_b %% 1 == 0 )
    stopifnot("Backimputation window should be shorter than observed incidence" =
        length(incid) >= window_b )
    
    
    if( window_b <= 5 ){
        msg <- "The backimputation window is short and may lead to an inaccurate estimate of the growth rate."
        warning(msg)
    }
    
    # process observed incidence
    incid_processed <- process_I(incid)
    incid_processed[1, ] <- c(sum(incid_processed[1, ]), 0)

    # some cases may be 0, implying -infinite logs
    safe_shift <- .5

    # backimpute unobserved, previous cases based on first window_b of observations
    log_incid_start <- data.frame(
        t = seq(window_b),
        li = log(incid[1:window_b] + safe_shift )
    )
    imputed_t <- seq(from = -100, to = 0)
    fit_backimpute <- lm(li ~ t, data = log_incid_start)

    if(fit_backimpute$coefficients[2] < 0) {
        msg <- "Estimate of the growth rate is negative, consider removing backimputation, or extending the backimputation window"
        warning(msg)
    }

    predict_backimpute_log <- predict.lm(fit_backimpute, newdata = list(t = imputed_t))
    predict_backimpute <- data.frame(
        local = exp(predict_backimpute_log) - safe_shift,
        imported = rep(0, length(imputed_t))
    )
    rownames(predict_backimpute) <- imputed_t

    # exclude negative cases arising from shift before logs.
    idx_nonnegative <- which(predict_backimpute$local >= 0)
    predict_backimpute <- predict_backimpute[ idx_nonnegative,  ]

    # first observation must be imported, otherwise later warning
    predict_backimpute[1, ] <- c(0, sum(predict_backimpute[1, ]))

    # bind imputed with observed
    incid_with_imputed <- rbind(predict_backimpute, incid_processed)
    return(incid_with_imputed)
}