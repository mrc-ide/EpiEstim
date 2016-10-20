##  The simulated data we want is:
##  1.  Incidence data, which is found in the matrix casesPerDay
##  2.  Interval data of infection times of individual's parents and the individual themselves, which is found in the matrix dataMatrix  

Rzero <- 2

simulateTransmissionModel <- function(Rzero){
  
  N <- 1000
  E <- 1
  I <- 0
  R <- 0
  S <- N - E - I - R
  
  gamma <- 1/2
  mu <- 1
  beta <- mu*Rzero/N
  
  dataMatrix <- matrix(NA, nrow = N, ncol = 5)
  dataMatrix <- as.data.frame(dataMatrix)
  names(dataMatrix) <- c("Individual ID", "I1L", "I1R", "I2L", "I2R")
  
  currentStatus <- matrix(0, nrow = N, ncol = 2)
  currentStatus <- as.data.frame(currentStatus)
  names(currentStatus) <- c("Individual ID", "Current status")  # 0=susceptible, 1=exposed, 2=infected, 3=removed
  currentStatus$`Current status`[1] = 1;
  
  dataMatrix$I1L[1] = NA
  dataMatrix$I1R[1] = NA  # Don't know when the infector of individual 1 became infected  
  
  casesPerDay <- matrix(0, nrow = 100000, ncol = 2)
  casesPerDay <- as.data.frame(casesPerDay)
  names(casesPerDay) <- c("Time", "Cases")
  
  casesPerDay$Cases[1] <- 0
  
  casesPerDay$Time[1:100000] = 0:(100000 - 1)
  dataMatrix$`Individual ID`[1:N] = 1:N
  currentStatus$`Individual ID`[1:N] = 1:N
  
  timings <- 0
  recordNum <- 1
  
  numberNewCasesSinceYesterday <- 0
  
  while (((E+I) > 0) && (recordNum < 100000)) {
    
    proposedTime <- timings - (1/(beta*I*S + gamma*E + mu*I))*log(runif(1, 0, 1))
    
    if (proposedTime < casesPerDay$Time[recordNum + 1]) {
      
      timings <- proposedTime
      
      r <- runif(1, 0, 1);
      
      if (r < beta*I*S/(beta*I*S + gamma*E + mu*I)) {
        
        S <- S - 1
        E <- E + 1
        
        # Pick an S to change to an E.  Also decide who infected them, and assign lower and upper bounds of infection times
        madeChange = 0
        while (madeChange < 0.5) {
          
          whichOne = ceiling(N*runif(1, 0, 1))
          
          
          if (currentStatus$`Current status`[whichOne] == 0) {
            
            madeChange = 1
            
            madeChange2 = 0
            while (madeChange2 < 0.5) {
              
              whichInfector = ceiling(N*runif(1, 0, 1))
              if (currentStatus$`Current status`[whichInfector] == 2) {
                dataMatrix$I1L[whichOne] = dataMatrix$I2L[whichInfector]
                dataMatrix$I1R[whichOne] = dataMatrix$I2R[whichInfector]
                madeChange2 = 1
              }
              
            }
            
            currentStatus$`Current status`[whichOne] = 1
            
          }
          
        }
        
      } else {
        
        if ((r > beta*I*S/(beta*I*S + gamma*E + mu*I)) && (r < (beta*I*S + gamma*E)/(beta*I*S + gamma*E + mu*I))) {
          E <- E - 1
          I <- I + 1
          numberNewCasesSinceYesterday <- numberNewCasesSinceYesterday + 1
          
          # Pick an E to change to an I
          madeChange = 0
          while (madeChange < 0.5) {
            
            whichOne = ceiling(N*runif(1, 0, 1))
            
            if (currentStatus$`Current status`[whichOne] == 1) {
              currentStatus$`Current status`[whichOne] = 2
              
              dataMatrix$I2L[whichOne] = floor(timings)
              dataMatrix$I2R[whichOne] = ceiling(timings)
              
              madeChange = 1
            }
            
          }
          
          
          
        } else {
          I <- I - 1
          R <- R + 1
          
          # Pick an I to change to an R
          madeChange = 0
          while (madeChange < 0.5) {
            
            whichOne = ceiling(N*runif(1, 0, 1))
            if (currentStatus$`Current status`[whichOne] == 2) {
              currentStatus$`Current status`[whichOne] = 3
              madeChange = 1
            }
            
          }
          
        }
        
      } } else {
        
        casesPerDay$Cases[recordNum + 1] <- numberNewCasesSinceYesterday
        
        numberNewCasesSinceYesterday <- 0
        
        timings <- casesPerDay$Time[recordNum + 1];
        recordNum <- recordNum + 1;
        
      }
    
    print(timings)
    
  }
  
  #plot(casesPerDay$Time[1:(recordNum - 1)], casesPerDay$Cases[1:(recordNum - 1)], xlab="Time (days)", ylab="Number of new cases", type='l')
  
  out <- list(casesPerDay = casesPerDay, dataMatrix = dataMatrix)
  
  return(out)
  
}
