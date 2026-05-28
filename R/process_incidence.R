process_I <- function(incid) {
  # If the input is an incidence object, we want to convert it to a data frame
  # that EpiEstim understands, which contains a single column for the I counts.
  ## TODO: remove this bit once estimate_R and wallinga_teunis have methods
  if (inherits(incid, "incidence")) {
    I_inc   <- incid
    incid   <- as.data.frame(I_inc)
    incid$I <- rowSums(incidence::get_counts(I_inc))
  }
  vector_I        <- FALSE
  single_col_df_I <- FALSE
  if (is.vector(incid)) {
    vector_I <- TRUE
  } else if (is.data.frame(incid) && ncol(incid) == 1) {
    single_col_df_I <- TRUE
  }
  if (vector_I || single_col_df_I) {
    if (single_col_df_I) {
      I_tmp <- incid[[1]]
    } else {
      I_tmp <- incid
    }
    incid      <- data.frame(local = I_tmp, imported = rep(0, length(I_tmp)))
    I_init     <- sum(incid[1, ])
    incid[1, ] <- c(0, I_init)
  } else {
    if (!is.data.frame(incid) || 
        (!("I" %in% names(incid)) &&
         !all(c("local", "imported") %in% names(incid)))) {
      stop("incid must be a vector or a dataframe with either i) a column 
           called 'I', or ii) 2 columns called 'local' and 'imported'.")
    }
    if (("I" %in% names(incid)) && 
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

  incid[is.na(incid)] <- 0
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

  incid
}

process_I_vector <- function(incid) {
  # incidence2 object
  if (inherits(incid, "incidence2")) {
    I_col <- incidence2::get_count_value_name(incid)
    date_col <- incidence2::get_date_index_name(incid)
    dates <- unique(incid[[date_col]])
    incid <- tapply(incid[[I_col]], incid[[date_col]], sum)[as.character(dates)]
    return(as.vector(incid))
  }

  incid_processed <- process_I(incid)
  as.vector(rowSums(incid_processed[, c("local", "imported")]))
}
