## ---- echo = FALSE-------------------------------------------------------

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>", 
  fig.width = 6, 
  fig.height = 6, 
  fig.path = "figs-demo/"
)


## ----data----------------------------------------------------------------
library(EpiEstim)
library(ggplot2)

## load data
data(Flu2009)
## incidence:
head(Flu2009$incidence)
## serial interval (SI) distribution:
Flu2009$si_distr
## interval-ceonsored serial interval data:
## each line represents a transmission event, 
## EL/ER show the lower/upper bound of the symptoms onset date in the infector
## SL/SR show the same for the secondary case
## type has entries 0 corresponding to doubly interval-censored data
## (see Reich et al. Statist. Med. 2009).
head(Flu2009$si_data)

## ---- incidence----------------------------------------------------------
library(incidence)
plot(as.incidence(Flu2009$incidence$I, dates = Flu2009$incidence$dates))

## ----estimate_r_parametric_si--------------------------------------------
res_parametric_si <- estimate_R(Flu2009$incidence, 
                                method="parametric_si",
                                config = make_config(list(
                                  mean_si = 2.6, 
                                  std_si = 1.5))
)

head(res_parametric_si$R)

## ----plot_r_parametric_si------------------------------------------------
plot(res_parametric_si, legend = FALSE)
# use `type = "R"`, `type = "incid"` or `type = "SI"` 
# to generate only one of the 3 plots

## ----estimate_r_non_parametric_si----------------------------------------
res_non_parametric_si <- estimate_R(Flu2009$incidence, 
                                    method="non_parametric_si",
                                    config = make_config(list(
                                      si_distr = Flu2009$si_distr))
)
# si_distr gives the probability mass function of the serial interval for 
# time intervals 0, 1, 2, etc.
plot(res_non_parametric_si, "R")

## ----discr_si------------------------------------------------------------
discr_si(0:20, mu = 2.6, sigma = 1.5)

## ----estimate_r_uncertain_si---------------------------------------------
## we choose to draw:
## - the mean of the SI in a Normal(2.6, 1), truncated at 1 and 4.2
## - the sd of the SI in a Normal(1.5, 0.5), truncated at 0.5 and 2.5
config <- make_config(list(mean_si = 2.6, std_mean_si = 1,
                           min_mean_si = 1, max_mean_si = 4.2,
                           std_si = 1.5, std_std_si = 0.5,
                           min_std_si = 0.5, max_std_si = 2.5))
res_uncertain_si <- estimate_R(Flu2009$incidence,
                               method = "uncertain_si",
                               config = config)
plot(res_uncertain_si, legend = FALSE) 
## the third plot now shows all the SI distributions considered

## ----si_data-------------------------------------------------------------
## interval-ceonsored serial interval data:
## each line represents a transmission event, 
## EL/ER show the lower/upper bound of the symptoms onset date in the infector
## SL/SR show the same for the secondary case
## type has entries 0 corresponding to doubly interval-censored data
## (see Reich et al. Statist. Med. 2009).
head(Flu2009$si_data)

## ----estimate_r_si_from_data---------------------------------------------
## fixing the random seeds
MCMC_seed <- 1
overall_seed <- 2
mcmc_control <- make_mcmc_control(seed = MCMC_seed, 
                                  burnin = 1000)
dist <- "G" # fitting a Gamma dsitribution for the SI
config <- make_config(list(si_parametric_distr = dist,
                           mcmc_control = mcmc_control,
                           seed = overall_seed, 
                           n1 = 50, 
                           n2 = 50))
res_si_from_data <- estimate_R(Flu2009$incidence,
                               method = "si_from_data",
                               si_data = Flu2009$si_data,
                               config = config)

plot(res_si_from_data, legend = FALSE)
## the third plot now shows the posterior sample of SI distributions 
## that were integrated over

## ----estimate_r_si_from_sample-------------------------------------------
## using the same random seeds as before to be able to compare results

## first estimate the SI distribution using function dic.fit.mcmc fron 
## coarseDataTools package:
n_mcmc_samples <- config$n1*mcmc_control$thin
SI_fit <- coarseDataTools::dic.fit.mcmc(dat = Flu2009$si_data,
                  dist = dist,
                  init.pars = init_mcmc_params(Flu2009$si_data, dist),
                  burnin = mcmc_control$burnin,
                  n.samples = n_mcmc_samples,
                  seed = mcmc_control$seed)
## thinning the output of the MCMC and converting using coarse2estim function
si_sample <- coarse2estim(SI_fit, thin = mcmc_control$thin)$si_sample
res_si_from_sample <- estimate_R(Flu2009$incidence,
                                method = "si_from_sample",
                                si_sample = si_sample,
                                config = make_config(list(n2 = 50, 
                                seed = overall_seed)))

## check that res_si_from_sample is the same as res_si_from_data
## since they were generated using the same MCMC algorithm to generate the SI
## sample (either internally to EpiEstim or externally)
all(res_si_from_sample$R$`Mean(R)` == res_si_from_data$R$`Mean(R)`)

## ----estimate_r_weekly---------------------------------------------------
T <- nrow(Flu2009$incidence)
t_start <- seq(2, T-6) # starting at 2 as conditional on the past observations
t_end <- t_start + 6 # adding 6 to get 7-day windows as bounds included in window
res_weekly <- estimate_R(Flu2009$incidence, 
                         method="parametric_si",
                         config = make_config(list(
                           t_start = t_start,
                           t_end = t_end,
                           mean_si = 2.6, 
                           std_si = 1.5))
)
plot(res_weekly, "R") 

## ----estimate_r_biweekly-------------------------------------------------
t_start <- seq(2, T-13) # starting at 2 as conditional on the past observations
t_end <- t_start + 13 
res_biweekly <- estimate_R(Flu2009$incidence, 
                           method="parametric_si",
                           config = make_config(list(
                             t_start = t_start,
                             t_end = t_end,
                             mean_si = 2.6, 
                             std_si = 1.5))
)
plot(res_biweekly, "R") 

## ----estimate_r_before_during_after_closure------------------------------
t_start <- c(2, 18, 25) # starting at 2 as conditional on the past observations
t_end <- c(17, 24, 32)
res_before_during_after_closure <- estimate_R(Flu2009$incidence, 
                                              method="parametric_si",
                                              config = make_config(list(
                                                t_start = t_start,
                                                t_end = t_end,
                                                mean_si = 2.6, 
                                                std_si = 1.5))
)
plot(res_before_during_after_closure, "R") +
  geom_hline(aes(yintercept = 1), color = "red", lty = 2)

## ----incid_table---------------------------------------------------------
head(Flu2009$incidence)
config <- make_config(list(mean_si = 2.6, std_si = 1.5))
res_incid_table <- estimate_R(Flu2009$incidence, 
                              method="parametric_si",
                              config = config)
plot(res_incid_table, "R")

## ----incid_vector--------------------------------------------------------
res_incid_vector <- estimate_R(Flu2009$incidence$I, 
                               method="parametric_si",
                               config = config)
plot(res_incid_vector, "R")

## ----line_list-----------------------------------------------------------
dates_onset <- Flu2009$incidence$dates[unlist(lapply(1:nrow(Flu2009$incidence), function(i) 
  rep(i, Flu2009$incidence$I[i])))]

## ----incid_class---------------------------------------------------------
last_date <- Flu2009$incidence$date[T]
res_incid_class <- estimate_R(incidence(dates_onset, last_date = last_date), 
                              method="parametric_si",
                              config = config)
plot(res_incid_class, "R")

## ----imports-------------------------------------------------------------
# generating fake information on our cases:
location <- sample(c("local","imported"), length(dates_onset), replace=TRUE)
location[1] <- "imported" # forcing the first case to be imported

## get incidence per group (location)
incid <- incidence(dates_onset, groups = location)

plot(incid)

## Estimate R with assumptions on serial interval
res_with_imports <- estimate_R(incid, method = "parametric_si",
                   config = make_config(list(
                   mean_si = 2.6, std_si = 1.5)))
plot(res_with_imports, add_imported_cases=TRUE)


