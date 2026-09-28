#!/usr/bin/env Rscript
#
# Usage (all arguments optional, "--key value" pairs):
#   Rscript run_episewer.R --prefix sim_ep --outdir results \
#       --shed_mean 8 --shed_sd 5 --gen_mean 3 --gen_sd 2.4 \
#       --iter_warmup 500 --iter_sampling 500 --chains 4 --cores 4 --seed 42
#
# --prefix        prefix of <prefix>_measurements.csv and <prefix>_flows.csv
# --outdir        folder for results (created if missing)
# --shed_mean/sd  mean and sd (days) of the shedding distribution; use the
#                 same values as in the simulation
# --gen_mean/sd   mean and sd (days) of the generation time distribution
#                 (placeholders: use the one matching your tree)
# --iter_warmup, --iter_sampling, --chains, --cores, --seed  Stan settings

suppressPackageStartupMessages({
    library(EpiSewer)
    library(data.table)
    library(ggplot2)
})

# ---- minimal argument parser (base R, no extra dependencies) --------------
defaults <- list(data_dir = "test",
                 iter_warmup = 500, iter_sampling = 500, chains = 4,
                 cores = 4, seed = 42)

parse_cli <- function(defaults) {
    args <- commandArgs(trailingOnly = TRUE)
    out <- defaults
    i <- 1
    while (i <= length(args)) {
        key <- gsub("-", "_", sub("^--", "", args[i]))
        if (!key %in% names(defaults) || i == length(args))
            stop("Unknown or incomplete argument: ", args[i])
        val <- args[i + 1]
        out[[key]] <- if (is.numeric(defaults[[key]])) as.numeric(val) else val
        i <- i + 2
    }
    out
}
opt <- parse_cli(defaults)

# ---- data -----------------------------------------------------------------
measurements <- fread(paste0(opt$data_dir, "/episewer_measurements.csv"))
flows        <- fread(paste0(opt$data_dir, "/episewer_flows.csv"))
measurements[, date := as.Date(date)]
measurements_sparse <- measurements[seq(3, nrow(measurements), by = 3), ]
flows[, date := as.Date(date)]

ww_data <- sewer_data(
    measurements = measurements,
    flows = flows,
    # cases = data_zurich$cases # cases are optional
)

ww_assumptions <- sewer_assumptions(
    generation_dist = get_discrete_gamma_shifted(gamma_mean = 3, gamma_sd = 2.4),
    shedding_dist = get_discrete_gamma(gamma_shape = 0.929639, gamma_scale = 7.241397),
    shedding_reference = "symptom_onset", # infection
    incubation_dist = get_discrete_gamma(gamma_shape = 8.5, gamma_scale = 0.4)
)

# ---- fit ------------------------------------------------------------------
options(mc.cores = opt$cores)
ww_result <- EpiSewer(
    data = ww_data,
    assumptions = ww_assumptions,
    fit_opts = set_fit_opts(sampler = sampler_stan_mcmc(
        iter_warmup = opt$iter_warmup, iter_sampling = opt$iter_sampling,
        chains = opt$chains, seed = opt$seed))
)

# ---- outputs ----------------------------------------------------------------
saveRDS(ww_result, file.path(opt$data_dir, "episewer_result.rds"))

ggsave(file.path(opt$data_dir, "R_estimate.png"), plot_R(ww_result),
       width = 8, height = 4, dpi = 150)
ggsave(file.path(opt$data_dir, "growth_rate.png"), plot_growth_rate(ww_result),
       width = 8, height = 4, dpi = 150)
ggsave(file.path(opt$data_dir, "concentration_fit.png"),
       plot_concentration(ww_result, measurements = measurements),
       width = 8, height = 4, dpi = 150)

# Rt summary as CSV. NOTE: the structure of ww_result$summary is not something
# I verified; if this fails, inspect names(ww_result$summary) or load the RDS.
tryCatch(
    fwrite(ww_result$summary$R, file.path(opt$data_dir, "R_estimate.csv")),
    error = function(e) message("Could not write R_estimate.csv: ", conditionMessage(e))
)

message("EpiSewer finished. Results in: ", opt$outdir)





# 
# 
# # measurements: EpiSewer requires a time series of concentration measurements. 
# #               Measurements should ideally be in gc/mL (gc = gene copies).
# # flows: Data about the daily wastewater volume flowing through the sampling site (should be in mL/day). 
# #        Flow data is highly recommended, but not strictly necessary (see documentation of flows_assume()).
# # cases: Data about confirmed cases in the catchment area. 
# #        This is optional but can help informing some of our prior assumptions.
# 
# # Basic workflow
# data_zurich <- SARS_CoV_2_Zurich
# names(data_zurich)
# 
# data_zurich$measurements
# 
# measurements_sparse <- data_zurich$measurements[,weekday := weekdays(data_zurich$measurements$date)][weekday %in% c("Monday","Thursday"),]
# head(measurements_sparse, 10)
# 
# data_zurich$flows
# 
# data_zurich$cases
# 
# ww_data <- sewer_data(
#     measurements = measurements_sparse,
#     flows = data_zurich$flows,
#     cases = data_zurich$cases # cases are optional
# )
# 
# ww_assumptions <- sewer_assumptions(
#     generation_dist = get_discrete_gamma_shifted(gamma_mean = 3, gamma_sd = 2.4),
#     shedding_dist = get_discrete_gamma(gamma_shape = 0.929639, gamma_scale = 7.241397),
#     shedding_reference = "symptom_onset", # shedding load distribution is relative to symptom onset
#     incubation_dist = get_discrete_gamma(gamma_shape = 8.5, gamma_scale = 0.4)
# )
# 
# options(mc.cores = 4) # allow stan to use 4 cores, i.e. one for each chain
# ww_result <- EpiSewer(
#     data = ww_data,
#     assumptions = ww_assumptions,
#     fit_opts = set_fit_opts(sampler = sampler_stan_mcmc(
#         iter_warmup = 500, iter_sampling = 500, chains = 4, seed = 42
#     ))
# )
# 
# plot_concentration(ww_result, measurements = data_zurich$measurements)
# 
# plot_concentration(ww_result, measurements = data_zurich$measurements, flows = data_zurich$flows, normalized = TRUE)
# 
# plot_R(ww_result)
# 
# plot_growth_report(ww_result)
# 
# plot_growth_report(ww_result, date = "2022-02-19")
# 
# plot_growth_rate(ww_result)
# 
# plot_doubling_time(ww_result)
# 
# plot_load(ww_result)
# plot_infections(ww_result) # Do not interpret as true absolute incidence or prevalence!!!
# plot_infections(ww_result) + geom_step(data = data_zurich$cases, aes(x = date, y = cases), color = "#000080")
# 
# plot_prior_posterior(ww_result, "measurement_noise_cv")
# 
# 
# 
# # Helper functions
# EpiSewer::component_functions("infection_noise")
# EpiSewer::component_functions("generation_dist")
# EpiSewer::component_functions("concentrations")
# EpiSewer::component_functions("R")
# 
# 
# # More precise workflow
# ww_measurements <- model_measurements(
#     concentrations = concentrations_observe(measurements = measurements_sparse),
#     noise = noise_estimate(),
#     LOD = LOD_none()
# )
# 
# ww_sampling <- model_sampling(
#     outliers = outliers_none(),
#     sample_effects = sample_effects_none()
# )
# 
# ww_sewage <- model_sewage(
#     flows = flows_observe(flows = data_zurich$flows),
#     residence_dist = residence_dist_assume(residence_dist = c(1))
# )
# 
# ww_shedding <- model_shedding(
#     shedding_dist = shedding_dist_assume(
#         get_discrete_gamma(gamma_shape = 0.929639, gamma_scale = 7.241397), shedding_reference = "symptom_onset"),
#     incubation_dist = incubation_dist_assume(get_discrete_gamma(gamma_shape = 8.5, gamma_scale = 0.4)),
#     load_per_case = load_per_case_calibrate(cases = data_zurich$cases),
#     load_variation = load_variation_estimate()
# )
# 
# ww_infections <- model_infections(
#     generation_dist = generation_dist_assume(get_discrete_gamma_shifted(gamma_mean = 3, gamma_sd = 2.4)),
#     R = R_estimate_gp(),
#     seeding = seeding_estimate_rw(),
#     infection_noise = infection_noise_estimate()
# )
# 
# ww_forecast <- model_forecast(
#     horizon = horizon_assume(horizon = 7),
#     damping = damping_assume(0.95)
# )
# 
# ww_fit_opts <- set_fit_opts(
#     model = model_stan_opts(package = "EpiSewer"),
#     sampler = sampler_stan_mcmc(
#         iter_warmup = 500,
#         iter_sampling = 500,
#         chains = 4,
#         parallel_chains = 4, # run all chains in parallel
#         seed = 42
#     )
# )
# 
# ww_results_opts <- set_results_opts(
#     fitted = TRUE,
#     summary_intervals = c(0.5, 0.95),
#     samples_ndraws = 50
# )
# 
# ww_result <- EpiSewer(
#     measurements = ww_measurements,
#     sampling = ww_sampling,
#     sewage = ww_sewage,
#     shedding = ww_shedding,
#     infections = ww_infections,
#     forecast = ww_forecast,
#     fit_opts = ww_fit_opts,
#     results_opts = ww_results_opts
# )
# 
# 
# ### TESTING
# 
# concentrations <- read_csv("episewer_measurements.csv")
# concentrations_sparse <- concentrations[seq(3, nrow(concentrations), by = 3), ]
# flows <- read_csv("episewer_flows.csv")
# 
# ww_data <- sewer_data(
#     measurements = concentrations_sparse,
#     flows = flows,
#     # cases = data_zurich$cases # cases are optional
# )
# 
# ww_assumptions <- sewer_assumptions(
#     generation_dist = get_discrete_gamma_shifted(gamma_mean = 3, gamma_sd = 2.4),
#     shedding_dist = get_discrete_gamma(gamma_shape = 0.929639, gamma_scale = 7.241397),
#     shedding_reference = "symptom_onset", # shedding load distribution is relative to symptom onset
#     incubation_dist = get_discrete_gamma(gamma_shape = 8.5, gamma_scale = 0.4)
# )
# 
# options(mc.cores = 4) # allow stan to use 4 cores, i.e. one for each chain
# ww_result <- EpiSewer(
#     data = ww_data,
#     assumptions = ww_assumptions,
#     fit_opts = set_fit_opts(sampler = sampler_stan_mcmc(
#         iter_warmup = 500, iter_sampling = 500, chains = 4, seed = 42
#     ))
# )
# 
# 
# plot_concentration(ww_result, measurements = concentrations)
# 
# plot_concentration(ww_result, measurements = concentrations, flows = flows, normalized = TRUE)
# 
# plot_R(ww_result)
# 
# plot_growth_report(ww_result)
# 
# plot_growth_report(ww_result, date = "2002-02-19")
# 
# plot_growth_rate(ww_result)
# 
# plot_doubling_time(ww_result)
# 
# plot_load(ww_result)
# plot_infections(ww_result) # Do not interpret as true absolute incidence or prevalence!!!
# # plot_infections(ww_result) + geom_step(data = data_zurich$cases, aes(x = date, y = cases), color = "#000080")
# 
# plot_prior_posterior(ww_result, "measurement_noise_cv")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
