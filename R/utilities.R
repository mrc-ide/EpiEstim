process_SI.Data <- function(SI.Data)
{
  # NULL entries
  if(is.null(SI.Data))
  {
    stop("Method SIFromData requires non NULL argument SI.Data") 
  }
  
  # wrong number of columns
  SI.Data <- as.data.frame(SI.Data)
  num_cols = dim(SI.Data)[2]
  if (num_cols < 4 || num_cols > 5) {
    stop("SI.Data should have 4 or 5 columns")
  }
  
  # entries with incorrect column names
  if(!all(c("EL", "ER", "SL", "SR") %in% names(SI.Data)))
  {
    names <- c("EL", "ER", "SL", "SR", "type")
    names(SI.Data) <- names[1:num_cols]
    warning("column names for SI.Data were not as expected; they were automatically interpreted as 'EL', 'ER', 'SL', 'SR', and 'type' (the last one only if SI.Data had five columns). ") 
  }
  
  # non integer entries in date columns
  if(!all(sapply(1:4, function(e) class(SI.Data[,e])=="integer" )))
  {
    stop("SI.Data has entries for which EL, ER, SL or SR are non integers.") 
  }
  
  # entries with wrong order in lower and upper bounds of dates
  if(any(SI.Data$ER-SI.Data$EL<0))
  {
    stop("SI.Data has entries for which ER<EL.")
  }
  if(any(SI.Data$SR-SI.Data$SL<0))
  {
    stop("SI.Data has entries for which SR<SL.")
  }
  
  # entries with negative serial interval
  if(any(SI.Data$SR-SI.Data$EL<=0))
  {
    stop("You cannot fit any of the supported distributions to this SI dataset, because for some data points the maximum serial interval is <=0.")
  }
  
  # check that the types [0: double censored, 1; single censored, 2: exact observation] are correctly specified, and if not present put them in.
  tmp_type <- 2 - rowSums(cbind(SI.Data$ER-SI.Data$EL!=0, SI.Data$SR-SI.Data$SL!=0))
  if(!("type" %in% names(SI.Data)))
  {
    warning("SI.Data contains no 'type' column. This is inferred automatically from the other columns.")
    SI.Data$type <- tmp_type
  }else if(any(is.na(SI.Data$type)) | !all(SI.Data$type == tmp_type))
  {
    warning("SI.Data contains unexpected entries in the 'type' column. This is inferred automatically from the other columns.")
    SI.Data$type <- tmp_type
  }
  
  return(SI.Data)
}


process_I <- function(I)
{
  if(class(I)=="incidence")
  {
    I_inc <- I
    I <- as.data.frame(I_inc)
    I$I <- rowSums(I_inc$counts)
  }
  vector_I <- FALSE
  single_col_df_I <- FALSE
  if(is.vector(I)) 
  {
    vector_I <- TRUE
  }else if(is.data.frame(I))
  {
    if(ncol(I)==1)
    {
      single_col_df_I <- TRUE
    }
  }
  if(vector_I | single_col_df_I)
  {
    if(single_col_df_I)
    {
      I_tmp <- as.vector(I[,1])
    }else
    {
      I_tmp <- I
    }
    I <- data.frame(local=I_tmp, imported=rep(0, length(I_tmp)))
    I_init <- sum(I[1,])
    I[1,] <- c(0, I_init)
  }else
  {
    if(!is.data.frame(I) | (!("I" %in% names(I)) & !all(c("local","imported") %in% names(I)) ) ) 
    {
      stop("I must be a vector or a dataframe with either i) a column called 'I', or ii) 2 columns called 'local' and 'imported'.")
    }
    if(("I" %in% names(I)) & !all(c("local","imported") %in% names(I)))
    {
      I$local <- I$I
      I$local[1] <- 0
      I$imported <- c(I$I[1], rep(0, nrow(I)-1))
    }
    if(I$local[1]>0)
    {
      warning("I$local[1] is >0 but must be 0, as all cases on the first time step are assumed imported. This is corrected automatically by cases being transferred to I$imported.")
      I_init <- sum(I[1,c('local','imported')])
      I[1,c('local','imported')] <- c(0, I_init)
    }
  }
  
  I[which(is.na(I))] <- 0
  date_col <- names(I)=='dates'
  if(any(date_col))
  {
    if(any(I[,!date_col]<0))
    {
      stop("I must contain only non negative integer values.")
    }
  }else
  {
    if(any(I<0))
    {
      stop("I must contain only non negative integer values.")
    }
  }
  
  return(I)
}

process_I_vector <- function(I)
{
  if(class(I)=="incidence")
  {
    I <- rowSums(I$counts)
  }
  if(!is.vector(I))
  {
    if(is.data.frame(I))
    {
      if(ncol(I)==1)
      {
        I <- as.vector(I[,1])
      }else if('I' %in% names(I))
      {
        I <- as.vector(I$I)
      }else if(!all(c('local', 'imported') %in% names(I)))
      {
        stop("I must be a vector or a dataframe with at least a column named 'I' or two columns named 'local' and 'imported'.")
      }
    }else
    {
      stop("I must be a vector or a dataframe with at least a column named 'I' or two columns named 'local' and 'imported'.")
    }
  }
  I[which(is.na(I))] <- 0
  date_col <- names(I)=='dates'
  if(any(date_col))
  {
    if(any(I[,!date_col]<0))
    {
      stop("I must contain only non negative integer values.")
    }
  }else
  {
    if(any(I<0))
    {
      stop("I must contain only non negative integer values.")
    }
  }
  
  return(I)
}

process_SI.Sample <- function(SI.Sample)
{
  if (is.null(SI.Sample)) {
    stop("method SIFromSample requires to specify the SI.Sample argument.")
  }
  
  SI.Sample <- as.matrix(SI.Sample)
  
  if (any(SI.Sample[1,] != 0)) {
    stop("method SIFromSample requires that SI.Sample[1,] contains only 0.")
  }
  if (any(SI.Sample < 0)) {
    stop("method SIFromSample requires that SI.Sample must contain only non negtaive values.")
  }
  if (any(abs(colSums(SI.Sample) - 1) > 0.01)) {
    stop("method SIFromSample requires the sum of each column in SI.Sample to be 1.")
  }
  
  return(SI.Sample)
}

check_times <- function(T.Start, T.End, T) # this only produces warnings and errors, does not return anything
{
  if (!is.vector(T.Start)) {
    stop("T.Start must be a vector.")
  }
  if (!is.vector(T.End)) {
    stop("T.End must be a vector.")
  }
  if (length(T.Start) != length(T.End)) {
    stop("T.Start and T.End must have the same length.")
  }
  if (any(T.Start > T.End)) {
    stop("T.Start[i] must be <= T.End[i] for all i.")
  }
  if (any(T.Start < 2 | T.Start > T | T.Start%%1 != 0 )) {
    stop("T.Start must be a vector of integers between 2 and the number of timesteps in I.")
  }
  if (any(T.End < 2 | T.End > T | T.End%%1 != 0)) {
    stop("T.End must be a vector of integers between 2 and the number of timesteps in I.")
  }
}

check_SI.Distr <- function(SI.Distr, sumToOne = c("error", "warning")) # this only produces warnings and errors, does not return anything
{
  sumToOne <- match.arg(sumToOne)
  if (is.null(SI.Distr)) {
    stop("SI.Distr argument missing.")
  }
  if (!is.vector(SI.Distr)) {
    stop("SI.Distr must be a vector.")
  }
  if (SI.Distr[1] != 0) {
    stop("SI.Distr should be so that SI.Distr[1] = 0.")
  }
  if (any(SI.Distr < 0)) {
    stop("SI.Distr must be a positive vector.")
  }
  if (abs(sum(SI.Distr) - 1) > 0.01) {
    if(sumToOne == "error") 
    {
      stop("SI.Distr must sum to 1.")
    }
    else if(sumToOne == "warning") 
    {
      warning("SI.Distr does not sum to 1.")
    }
  }
}

check_dates <- function(I)
{
  dates <- I$dates
  if(class(dates) != "Date" & class(dates) != "numeric")
  {
    stop("I$dates must be an object of class date or numeric.")
  }else
  {
    if(unique(diff(dates)) != 1)
    {
      stop("I$dates must contain dates which are all in a row.")
    }else
    {
      return(dates)
    }
  }
}