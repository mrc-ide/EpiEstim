test_that("process_I copes with NA inputs", {
  
  ## Incidence as vector without or with NAs
  incid <- c(25,19,8,10,13,0,0)
  processed_incid <- process_I(incid)
  expect_true(all(!is.na(processed_incid)))
  
  incid <- c(NA,19,8,10,13,0,0)
  processed_incid <- process_I(incid)
  expect_true(all(!is.na(processed_incid)))
  
  incid <- c(25,19,8,10,13,NA,NA)
  processed_incid <- process_I(incid)
  expect_true(all(!is.na(processed_incid)))
  
  incid2 <- c(25,19,8,NA,10,13,NA,NA)
  processed_incid2 <- process_I(incid2)
  expect_true(all(!is.na(processed_incid2)))
  
  ## Incidence as dataframe without or with NAs
  incid_df <- data.frame(local = c(0, incid), imported = c(1, rep(0, length(incid))))
  incid_df$imported[3] <- NA
  processed_incid_df <- process_I(incid_df)
  expect_true(all(!is.na(processed_incid_df)))
  
  incid_df2 <- data.frame(local = c(0, incid2), imported = c(1, rep(0, length(incid2))))
  incid_df2$imported[c(3, 7)] <- NA
  processed_incid_df2 <- process_I(incid_df2)
  expect_true(all(!is.na(processed_incid_df2)))
  
  incid_df3 <- data.frame(I = incid)
  incid_df3$I[c(3, 7)] <- NA
  processed_incid_df3 <- process_I(incid_df3)
  expect_true(all(!is.na(processed_incid_df3)))
  
  ## Incidence as incidence object
  data <- c(0,1,1,2,1,3,4,5,5,5,5,4,NA,26,6,7,9)
  location <- sample(c("local","imported"), length(data), replace=TRUE)
  location[1] <- "imported" # forcing the first case to be imported
  
  ## get incidence per group (location)
  incid_incidence <- incidence(data, groups = location)
  processed_incid_incidence <- process_I(incid_incidence)
  expect_true(all(!is.na(processed_incid_incidence)))
  
  incid_incidence2 <- incid_incidence
  incid_incidence2$counts[12,"imported"] <- NA
  processed_incid_incidence2 <- process_I(incid_incidence2)
  expect_true(all(!is.na(processed_incid_incidence2)))
  
})


