process_si_data <- function(si_data) {
  # NULL entries
  if (is.null(si_data)) {
    stop("Method si_from_data requires non NULL argument si_data")
  }

  # wrong number of columns
  si_data <- as.data.frame(si_data)
  num_cols <- dim(si_data)[2]
  if (num_cols < 4 || num_cols > 5) {
    stop("si_data should have 4 or 5 columns")
  }

  # entries with incorrect column names
  if (!all(c("EL", "ER", "SL", "SR") %in% names(si_data))) {
    names <- c("EL", "ER", "SL", "SR", "type")
    names(si_data) <- names[seq_len(num_cols)]
    warning("column names for si_data were not as expected; they were 
            automatically interpreted as 'EL', 'ER', 'SL', 'SR', and 'type' 
            (the last one only if si_data had five columns). ")
  }

  # non integer entries in date columns
  if (!all(vlapply(seq_len(4), function(e) class(si_data[, e]) == "integer"))) {
    stop("si_data has entries for which EL, ER, SL or SR are non integers.")
  }

  # entries with wrong order in lower and upper bounds of dates
  if (any(si_data$ER - si_data$EL < 0)) {
    stop("si_data has entries for which ER<EL.")
  }
  if (any(si_data$SR - si_data$SL < 0)) {
    stop("si_data has entries for which SR<SL.")
  }

  # entries with negative serial interval
  if (any(si_data$SR - si_data$EL <= 0)) {
    stop("You cannot fit any of the supported distributions to this SI dataset, 
         because for some data points the maximum serial interval is <=0.")
  }

  ## check that the types [0: double censored, 1; single censored, 
  ## 2: exact observation] are correctly specified, and if not present 
  ## put them in.
  tmp_type <- 2 - rowSums(cbind(si_data$ER - si_data$EL != 0, 
                                si_data$SR - si_data$SL != 0))
  if (!("type" %in% names(si_data))) {
    warning("si_data contains no 'type' column. This is inferred automatically 
            from the other columns.")
    si_data$type <- tmp_type
  } else if (any(is.na(si_data$type)) | !all(si_data$type == tmp_type)) {
    warning("si_data contains unexpected entries in the 'type' column. This is 
            inferred automatically from the other columns.")
    si_data$type <- tmp_type
  }

  return(si_data)
}


process_I <- function(incid) {
  # If the input is an incidence object, we want to convert it to a data frame
  # that EpiEstim understands, which contains a single column for the I counts.
  if (inherits(incid, "incidence")) {
    I_inc   <- incid
    incid   <- as.data.frame(I_inc)
    incid$I <- rowSums(incidence::get_counts(I_inc))
  }
  vector_I        <- FALSE
  single_col_df_I <- FALSE
  if (is.vector(incid)) {
    vector_I <- TRUE
  } else if (is.data.frame(incid)) {
    if (ncol(incid) == 1) {
      single_col_df_I <- TRUE
    }
  }
  if (vector_I | single_col_df_I) {
    if (single_col_df_I) {
      I_tmp <- incid[[1]]
    } else {
      I_tmp <- incid
    }
    incid      <- data.frame(local = I_tmp, imported = rep(0, length(I_tmp)))
    I_init     <- sum(incid[1, ])
    incid[1, ] <- c(0, I_init)
  } else {
    if (!is.data.frame(incid) | 
        (!("I" %in% names(incid)) &
         !all(c("local", "imported") %in% names(incid)))) {
      stop("incid must be a vector or a dataframe with either i) a column 
           called 'I', or ii) 2 columns called 'local' and 'imported'.")
    }
    if (("I" %in% names(incid)) & 
        !all(c("local", "imported") %in% names(incid))) {
      incid$local    <- incid$I
      incid$local[1] <- 0
      incid$imported <- c(incid$I[1], rep(0, nrow(incid) - 1))
    }
    if (incid$local[1] > 0) {
      warning("incid$local[1] is >0 but must be 0, as all cases on the first 
              time step are assumed imported. This is corrected automatically 
              by cases being transferred to incid$imported.")
      I_init <- sum(incid[1, c("local", "imported")])
      incid[1, c("local", "imported")] <- c(0, I_init)
    }
  }

  incid[which(is.na(incid))] <- 0
  date_col <- names(incid) == "dates"
  if (any(date_col)) {
    if (any(incid[, !date_col] < 0)) {
      stop("incid must contain only non negative integer values.")
    }
  } else {
    if (any(incid < 0)) {
      stop("incid must contain only non negative integer values.")
    }
  }

  return(incid)
}

process_I_vector <- function(incid) {
  # here, the incident counts are being forced into a vector.
  if (inherits(incid, "incidence")) {
    incid <- rowSums(incidence::get_counts(incid))
  }
  if (!is.vector(incid)) {
    if (is.data.frame(incid)) {
      if (ncol(incid) == 1) {
        incid <- as.vector(incid[, 1])
      } else if ("I" %in% names(incid)) {
        incid <- as.vector(incid$I)
      } else if (!all(c("local", "imported") %in% names(incid))) {
        stop("incid must be a vector or a dataframe with at least a column named
             'I' or two columns named 'local' and 'imported'.")
      }
    } else {
      stop("incid must be a vector or a dataframe with at least a column named 
           'I' or two columns named 'local' and 'imported'.")
    }
  }
  incid[which(is.na(incid))] <- 0
  date_col <- names(incid) == "dates"
  if (any(date_col)) {
    if (any(incid[, !date_col] < 0)) {
      stop("incid must contain only non negative integer values.")
    }
  } else {
    if (any(incid < 0)) {
      stop("incid must contain only non negative integer values.")
    }
  }

  return(incid)
}

process_si_sample <- function(si_sample) {
  if (is.null(si_sample)) {
    stop("method si_from_sample requires to specify the si_sample argument.")
  }

  si_sample <- as.matrix(si_sample)

  if (any(si_sample[1, ] != 0)) {
    stop("method si_from_sample requires that si_sample[1,] contains only 0.")
  }
  if (any(si_sample < 0)) {
    stop("method si_from_sample requires that si_sample must contain only non 
         negtaive values.")
  }
  if (any(abs(colSums(si_sample) - 1) > 0.01)) {
    stop("method si_from_sample requires the sum of each column in si_sample to 
         be 1.")
  }

  return(si_sample)
}

check_times <- function(t_start, t_end, T) 
 ## this only produces warnings and errors, does not return anything
{
  if (!is.vector(t_start)) {
    stop("t_start must be a vector.")
  }
  if (!is.vector(t_end)) {
    stop("t_end must be a vector.")
  }
  if (length(t_start) != length(t_end)) {
    stop("t_start and t_end must have the same length.")
  }
  if (any(t_start > t_end)) {
    stop("t_start[i] must be <= t_end[i] for all i.")
  }
  if (any(t_start < 2 | t_start > T | t_start %% 1 != 0)) {
    stop("t_start must be a vector of integers between 2 and the number of 
         timesteps in incid.")
  }
  if (any(t_end < 2 | t_end > T | t_end %% 1 != 0)) {
    stop("t_end must be a vector of integers between 2 and the number of 
         timesteps in incid.")
  }
}

check_si_distr <- function(si_distr, sumToOne = c("error", "warning"), 
                           method = "non_parametric_si") 
  ## this only produces warnings and errors, does not return anything
{
  sumToOne <- match.arg(sumToOne)
  if (is.null(si_distr)) {
    stop(paste0("si_distr argument is missing but is required for method ", 
                method, "."))
  }
  if (!is.vector(si_distr)) {
    stop("si_distr must be a vector.")
  }
  if (si_distr[1] != 0) {
    stop("si_distr should be so that si_distr[1] = 0.")
  }
  if (any(si_distr < 0)) {
    stop("si_distr must be a positive vector.")
  }
  if (abs(sum(si_distr) - 1) > 0.01) {
    if (sumToOne == "error") {
      stop("si_distr must sum to 1.")
    }
    else if (sumToOne == "warning") {
      warning("si_distr does not sum to 1.")
    }
  }
}

check_dates <- function(incid) {
  dates <- incid$dates
  if (!inherits(dates, "Date") && !is.numeric(dates)) {
    stop("incid$dates must be an object of class date or numeric.")
  } else {
    if (unique(diff(dates)) != 1) {
      stop("incid$dates must contain dates which are all in a row.")
    } else {
      return(dates)
    }
  }
}

process_config <- function(config) {
  if (!("mean_prior" %in% names(config))) {
    config$mean_prior <- 5
  }

  if (!("std_prior" %in% names(config))) {
    config$std_prior <- 5
  }

  if (config$mean_prior <= 0) {
    stop("config$mean_prior must be >0.")
  }
  if (config$std_prior <= 0) {
    stop("config$std_prior must be >0.")
  }

  if (!("cv_posterior" %in% names(config))) {
    config$cv_posterior <- 0.3
  }

  if (!("mcmc_control" %in% names(config))) {
    config$mcmc_control <- make_mcmc_control()
  }

  return(config)
}

process_config_si_from_data <- function(config, si_data) {
  config$si_parametric_distr <- match.arg(
    config$si_parametric_distr,
    c("G", "W", "L", "off1G", "off1W", "off1L")
  )
  if (is.null(config$n1)) {
    stop("method si_from_data requires to specify the config$n1 argument.")
  }
  if (is.null(config$n2)) {
    stop("method si_from_data requires to specify the config$n2 argument.")
  }
  if (config$n2 <= 0 || config$n2 %% 1 != 0) {
    stop("method si_from_data requires a >0 integer value for config$n2.")
  }
  if (config$n1 <= 0 || config$n1 %% 1 != 0) {
    stop("method si_from_data requires a >0 integer value for config$n1.")
  }
  if (is.null(config$mcmc_control$init_pars)) {
    config$mcmc_control$init_pars <-
      init_mcmc_params(si_data, config$si_parametric_distr)
  }
  if ((config$si_parametric_distr == "off1G" |
    config$si_parametric_distr == "off1W" |
    config$si_parametric_distr == "off1L") &
    any(si_data$SR - si_data$EL <= 1)) {
    stop(paste(
      "You cannot fit a distribution with offset 1 to this SI",
      "dataset, because for some data points the maximum serial",
      "interval is <=1.\nChoose a different distribution"
    ))
  }
  return(config)
}

check_config <- function(config, method) {
  if (method == "non_parametric_si") {
    check_si_distr(config$si_distr, method = method)
  }
  if (method == "parametric_si") {
    if (is.null(config$mean_si)) {
      stop("method parametric_si requires to specify the config$mean_si 
           argument.")
    }
    if (is.null(config$std_si)) {
      stop("method parametric_si requires to specify the config$std_si 
           argument.")
    }
    if (config$mean_si <= 1) {
      stop("method parametric_si requires a value >1 for config$mean_si.")
    }
    if (config$std_si <= 0) {
      stop("method parametric_si requires a >0 value for config$std_si.")
    }
  }
  if (method == "uncertain_si") {
    if (is.null(config$mean_si)) {
      stop("method uncertain_si requires to specify the config$mean_si 
           argument.")
    }
    if (is.null(config$std_si)) {
      stop("method uncertain_si requires to specify the config$std_si 
           argument.")
    }
    if (is.null(config$n1)) {
      stop("method uncertain_si requires to specify the config$n1 argument.")
    }
    if (is.null(config$n2)) {
      stop("method uncertain_si requires to specify the config$n2 argument.")
    }
    if (is.null(config$std_mean_si)) {
      stop("method uncertain_si requires to specify the config$std_mean_si 
           argument.")
    }
    if (is.null(config$min_mean_si)) {
      stop("method uncertain_si requires to specify the config$min_mean_si 
           argument.")
    }
    if (is.null(config$max_mean_si)) {
      stop("method uncertain_si requires to specify the config$max_mean_si 
           argument.")
    }
    if (is.null(config$std_std_si)) {
      stop("method uncertain_si requires to specify the config$std_std_si 
           argument.")
    }
    if (is.null(config$min_std_si)) {
      stop("method uncertain_si requires to specify the config$min_std_si 
           argument.")
    }
    if (is.null(config$max_std_si)) {
      stop("method uncertain_si requires to specify the config$max_std_si 
           argument.")
    }
    if (config$mean_si <= 0) {
      stop("method uncertain_si requires a >0 value for config$mean_si.")
    }
    if (config$std_si <= 0) {
      stop("method uncertain_si requires a >0 value for config$std_si.")
    }
    if (config$n2 <= 0 || config$n2 %% 1 != 0) {
      stop("method uncertain_si requires a >0 integer value for config$n2.")
    }
    if (config$n1 <= 0 || config$n1 %% 1 != 0) {
      stop("method uncertain_si requires a >0 integer value for config$n1.")
    }
    if (config$std_mean_si <= 0) {
      stop("method uncertain_si requires a >0 value for config$std_mean_si.")
    }
    if (config$min_mean_si < 1) {
      stop("method uncertain_si requires a value >=1 for config$min_mean_si.")
    }
    if (config$max_mean_si < config$mean_si) {
      stop("method uncertain_si requires that config$max_mean_si >= 
           config$mean_si.")
    }
    if (config$mean_si < config$min_mean_si) {
      stop("method uncertain_si requires that config$mean_si >= 
           config$min_mean_si.")
    }
    if (signif(config$max_mean_si - config$mean_si, 3) != signif(config$mean_si -
      config$min_mean_si, 3)) {
      warning("The distribution you chose for the mean SI is not centered around
              the mean.")
    }
    if (config$std_std_si <= 0) {
      stop("method uncertain_si requires a >0 value for config$std_std_si.")
    }
    if (config$min_std_si <= 0) {
      stop("method uncertain_si requires a >0 value for config$min_std_si.")
    }
    if (config$max_std_si < config$std_si) {
      stop("method uncertain_si requires that config$max_std_si >= 
           config$std_si.")
    }
    if (config$std_si < config$min_std_si) {
      stop("method uncertain_si requires that config$std_si >= 
           config$min_std_si.")
    }
    if (signif(config$max_std_si - config$std_si, 3) != signif(config$std_si -
      config$min_std_si, 3)) {
      warning("The distribution you chose for the std of the SI is not centered 
              around the mean.")
    }
  }
  if (config$cv_posterior < 0) {
    stop("config$cv_posterior must be >0.")
  }
}

viapply <- function(X, FUN, ...) {
  vapply(X, FUN, integer(1), ...)
}

vlapply <- function(X, FUN, ...) {
  vapply(X, FUN, logical(1), ...)
}

vnapply <- function(X, FUN, ...) {
  vapply(X, FUN, numeric(1), ...)
}

vcapply <- function(X, FUN, ...) {
  vapply(X, FUN, character(1), ...)
}

## This function was contributed by Rich Fitzjohn. It modifies default arguments
## using user-provided values. The argument 'strict' triggers and error
## behaviour: if strict==TRUE: all new values need to be part of the defaults.

modify_defaults <- function(defaults, x, strict = TRUE) {
  extra <- setdiff(names(x), names(defaults))
  if (strict && (length(extra) > 0L)) {
    stop("Additional invalid options: ", paste(extra, collapse=", "))
  }
  utils::modifyList(defaults, x, keep.null = TRUE) # keep.null is needed here
}
