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