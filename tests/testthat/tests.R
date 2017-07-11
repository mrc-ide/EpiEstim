library(EpiEstim)
library(testthat)
library(compare)
library(incidence)

# The following compare_output function makes writing tests very easy.
# The first argument should be the output of a call to EpiEstim.
# The second should be a unique ID used to save the data to a file.

# When a new ID is used for the first time, the test will fail, but will
# also write the output from EpiEstim to a file with that ID. This means that
# next time you run the test it will pass. (It fails to begin with in case
# this was unintentional).

# If a test is failing, but you're sure the change is correct, simply delete
# the file with that ID and let it be re-written. (You'll have to run tests
# twice before they pass again).

# Ensure the correct saved files are always checked in to version control.

compare_output <- function(output, id) {
  filename <- paste("../expected_output/", id, ".RData", sep="")
  if(file.exists(filename)) {
    load(filename)
    expect_true(compare(saved, output)$result)
  } else {
    saved <- output
    save(saved, file=filename)
    stop(paste(sep="", "An existing file was not found matching '",
                  id,
                  "', so a new one has been created. In the future, tests will compare against this file"))
  }
}

test_that("Can import data from EpiEstim", {
  data("Flu2009")
})

# The following examples are from EpiEstim's documentation.
set.seed(1)

test_that("Example 1 matches saved output", {
  out <- EstimateR(Flu2009$Incidence, T.Start=2:26, T.End=8:32, method="NonParametricSI", SI.Distr=Flu2009$SI.Distr, plot=FALSE, seed=1)
  compare_output(out, "Example1")
})

test_that("Example 2 matches saved output", {
  data <- data <- c(0,1,1,2,1,3,4,5,5,5,5,4,4,26,6,7,9)
  location <- sample(c("local","imported"), length(data), replace=TRUE)
  location[1] <- "imported" # forcing the first case to be imported
  # get incidence per group (location)
  I <- incidence(data, groups = location)
  out <- EstimateR(I, T.Start = 2:21, T.End = 8:27, method = "ParametricSI", Mean.SI = 2.6, Std.SI = 1.5, plot = FALSE, seed=1)
  compare_output(out, "Example2")
})

test_that("Example 3 matches saved output", {
  ## estimate the reproduction number (method "ParametricSI")
  out <- EstimateR(Flu2009$Incidence, T.Start=2:26, T.End=8:32, method="ParametricSI", Mean.SI=2.6, Std.SI=1.5, plot=FALSE, seed=1)
  compare_output(out, "Example3")
})

test_that("Example 4 matches saved output", {
  ## estimate the reproduction number (method "UncertainSI")
  out <- EstimateR(Flu2009$Incidence, T.Start=2:26, T.End=8:32, method="UncertainSI",
           Mean.SI=2.6, Std.Mean.SI=1, Min.Mean.SI=1, Max.Mean.SI=4.2,
           Std.SI=1.5, Std.Std.SI=0.5, Min.Std.SI=0.5, Max.Std.SI=2.5,
           n1=100, n2=100, plot=FALSE, seed=1)
  compare_output(out, "Example4")
})



