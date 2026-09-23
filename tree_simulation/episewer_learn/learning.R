# measurements: EpiSewer requires a time series of concentration measurements. 
#               Measurements should ideally be in gc/mL (gc = gene copies).
# flows: Data about the daily wastewater volume flowing through the sampling site (should be in mL/day). 
#        Flow data is highly recommended, but not strictly necessary (see documentation of flows_assume()).
# cases: Data about confirmed cases in the catchment area. 
#        This is optional but can help informing some of our prior assumptions.

data_zurich <- SARS_CoV_2_Zurich
names(data_zurich)

data_zurich$measurements

measurements_sparse <- data_zurich$measurements[,weekday := weekdays(data_zurich$measurements$date)][weekday %in% c("Monday","Thursday"),]
head(measurements_sparse, 10)

data_zurich$flows

data_zurich$cases

ww_data <- sewer_data(
    measurements = measurements_sparse,
    flows = data_zurich$flows,
    cases = data_zurich$cases # cases are optional
)

ww_assumptions <- sewer_assumptions(
    generation_dist = get_discrete_gamma_shifted(gamma_mean = 3, gamma_sd = 2.4),
    shedding_dist = get_discrete_gamma(gamma_shape = 0.929639, gamma_scale = 7.241397),
    shedding_reference = "symptom_onset", # shedding load distribution is relative to symptom onset
    incubation_dist = get_discrete_gamma(gamma_shape = 8.5, gamma_scale = 0.4)
)

options(mc.cores = 4) # allow stan to use 4 cores, i.e. one for each chain
ww_result <- EpiSewer(
    data = ww_data,
    assumptions = ww_assumptions,
    fit_opts = set_fit_opts(sampler = sampler_stan_mcmc(
        iter_warmup = 500, iter_sampling = 500, chains = 4, seed = 42
    ))
)

plot_concentration(ww_result, measurements = data_zurich$measurements)

plot_concentration(ww_result, measurements = data_zurich$measurements, flows = data_zurich$flows, normalized = TRUE)

plot_R(ww_result)

plot_growth_report(ww_result)

plot_growth_report(ww_result, date = "2022-02-19")

plot_growth_rate(ww_result)

plot_doubling_time(ww_result)

plot_load(ww_result)
plot_infections(ww_result) # Do not interpret as true absolute incidence or prevalence!!!
plot_infections(ww_result) + geom_step(data = data_zurich$cases, aes(x = date, y = cases), color = "#000080")

plot_prior_posterior(ww_result, "measurement_noise_cv")



# Helper functions
EpiSewer::component_functions("infection_noise")
EpiSewer::component_functions("generation_dist")
EpiSewer::component_functions("concentrations")
EpiSewer::component_functions("R")










