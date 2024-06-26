# get imported
incid_raw <- as.Date("2001-01-01") + c(0, 1, 5, 5, 6, 7, 7, 7, 9, 10, 12, 13)
location <- sample(c("local", "imported"), length(incid_raw), replace = TRUE)
location[1] <- "imported"
incid_imported <- incidence::incidence(incid_raw, groups = location)

data("covid_deaths_2020_uk")

incid_covid <- covid_deaths_2020_uk$incidence$Incidence
config_covid <- make_config(list(si_distr = covid_deaths_2020_uk$si_distr))

incid_constant <- rep(10, 20)
config_constant <- make_config(list(mean_si = 5, std_si = 2))



test_that("warnings and errors are working as expected", {

    expect_error(
        estimate_R(incid = incid_covid, backimputation_window = 1),
        "Backimputation window needs to contain at least 2 timepoints"
    )

    expect_error(
        estimate_R(incid = incid_covid, backimputation_window = 10.4),
        "Backimputation window needs to have integer length"
    )

    expect_error(
        estimate_R(incid = incid_covid[1:10], backimputation_window = 100L),
        "Backimputation window should be shorter than observed incidence"
    )

    expect_error(
        estimate_R(incid = incid_covid[1:10], backimputation_window = 100L),
        "Backimputation window should be shorter than observed incidence"
    )

    expect_error(
        estimate_R(incid = incid_imported, backimputation_window = 7L),
        "incidence objects are currently not supported by backimpute_I()."
    )

  expect_warning(
        estimate_R(incid = incid_covid, backimputation_window = 5, config = config_covid),
        "Estimate of the growth rate is negative, consider removing backimputation, or extending the backimputation window"
    )

    expect_warning(
        estimate_R(incid = incid_covid, backimputation_window = 3, config = config_covid),
        "The backimputation window is short and may lead to an inaccurate estimate of the growth rate."
    )
})

test_that("outputs are working as expected", {

    # check default is 0.
    expect_equal(
        estimate_R(incid = incid_covid, backimputation_window = 0, config = config_covid),
        estimate_R(incid = incid_covid, config = config_covid)
    )

    # check adjusted estimates are lower than original...
    expect_lt({
        adjusted <- estimate_R(
            incid = incid_covid[1:10],
            backimputation_window = 10,
            config = config_covid)$R[1, "Mean(R)"]
        original <- estimate_R(
            incid = incid_covid[1:10],
            config = config_covid)$R[1, "Mean(R)"]
        adjusted-original
    },
        0
    )

    # ... even if growth rate is negative!
    suppressWarnings(
        expect_lt({
            adjusted <- estimate_R(
                incid = incid_covid,
                backimputation_window = 5,
                config = config_covid)$R[1, "Mean(R)"]
            original <- estimate_R(
                incid = incid_covid,
                config = config_covid)$R[1, "Mean(R)"]
            adjusted-original
        }, 0)
    )
    
    expect_lt({
      R_constant <- estimate_R(
        incid = incid_constant,
        backimputation_window = 7,
        config = config_constant,
        method = "parametric_si")
      diff(range(R_constant$R$`Mean(R)`)) 
    }, 1e-6)

})

