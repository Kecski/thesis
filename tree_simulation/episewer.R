#!/usr/bin/env Rscript
# Run EpiSewer on the CSV files written by run_simulation.py.
#
# Usage:
#   Rscript episewer.R config.yml [--data_dir data/wastewater] [--tree_id 3] \
#       [--sequencing_delay 5]
#
# Most parameters come from the YAML config file (first argument). All config
# keys are optional and fall back to the defaults defined in this script;
# unknown keys cause an error. See config_example.yml.
#
# Config file keys:
#   shed_mean/sd    mean and sd (days) of the shedding distribution; use the
#                   same values as in the simulation
#   gen_mean/sd     mean and sd (days) of the generation time distribution
#   incub_mean/sd   mean and sd (days) of the incubation period distribution
#   load_per_case   particles shed by a single infected
#   iter_warmup, iter_sampling, chains, cores, seed   Stan settings
#   num_segments    number of equally long sections (used if breakpoints is empty)
#   breakpoints     list of dates (or comma-separated string) defining the sections
#
# Command-line arguments (these are NOT read from the config file):
#   --data_dir         folder containing the tree_<id> subfolders; results are
#                      written into <data_dir>/tree_<id> (default data/wastewater)
#   --tree_id          which tree_<id> folder to analyse (default 0)
#   --sequencing_delay days to skip at the start when building equal sections
#                      (default 0)

# ----- EpiSewer base parameters -----
# measurements:
#     concentrations - observe from data
#     noise - estimate from data
#     LOD - none
# sampling:
#     outliers - estimate from data
#     sample_effects - none
# sewage:
#     flows - observe from data
#     residence_dist - assumes that particles arrive at the sampling site within the day of shedding
# shedding:
#     incubation_dist - assumes discrete gamma distribution with mean=3.4, sd=1.17 (alpha=8.5, beta=0.4)
#     shedding_dist - assumes discrete gamma distribution with mean=8, sd=5 (alpha=2.56, beta=3.125)
#     load_per_case - assumes 1e12 particles shed by a single infected
#     load_variation - estimate from data
# infections:
#     generation_dist - assumes discrete gamma distribution with mean=3, sd=2.4 (alpha=1.6525, beta=1.92)
#     R - estimate R from data using Gaussian proccesses
#     seeding - estimate initial infections with geometric random walk
#     infection_noise - estimate from data

suppressPackageStartupMessages({
    library(EpiSewer)
    library(data.table)
    library(ggplot2)
    library(yaml)
})

# ---- config + command-line parsing ----------------------------------------
# breakpoints: dates defining n sections. n+1 dates give n half-open sections
# [breakpoints[i], breakpoints[i+1]). Either a YAML list of dates or a
# comma-separated string, e.g. "2020-03-01,2020-04-01,2020-05-15,2020-07-01"
config_defaults <- list(shed_mean = 8, shed_sd = 5, gen_mean = 3, gen_sd = 2.4,
                        incub_mean = 3.4, incub_sd = 1.17, load_per_case = 1e10,
                        iter_warmup = 500, iter_sampling = 500, chains = 4,
                        cores = 4, seed = 42, num_segments = 1,
                        breakpoints = "") # format for breakpoints: "YYYY-MM-DD,YYYY-MM-DD,..." 
cli_defaults <- list(data_dir = "data/wastewater", tree_id = 0,
                     sequencing_delay = 0)

read_config <- function(config_file, defaults, cli_keys) {
    if (!file.exists(config_file))
        stop("Config file not found: ", config_file)
    cfg <- yaml::read_yaml(config_file)
    if (is.null(cfg)) cfg <- list()  # empty file -> all defaults
    
    misplaced <- intersect(names(cfg), cli_keys)
    if (length(misplaced) > 0)
        stop("These keys must be given on the command line, not in the config ",
             "file: ", paste(misplaced, collapse = ", "))
    unknown <- setdiff(names(cfg), names(defaults))
    if (length(unknown) > 0)
        stop("Unknown key(s) in config file: ", paste(unknown, collapse = ", "),
             ". Allowed keys: ", paste(names(defaults), collapse = ", "))
    
    out <- defaults
    for (key in names(cfg)) {
        val <- cfg[[key]]
        if (key == "breakpoints") {
            # accept a YAML list of dates or a comma-separated string
            val <- paste(as.character(unlist(val)), collapse = ",")
        } else if (is.numeric(defaults[[key]])) {
            # as.numeric also handles values like 1e10 that YAML may read as text
            val <- suppressWarnings(as.numeric(val))
            if (length(val) != 1 || is.na(val))
                stop("Config key '", key, "' must be a single number.")
        }
        out[[key]] <- val
    }
    out
}

parse_cli <- function(args, defaults) {
    out <- defaults
    i <- 1
    while (i <= length(args)) {
        key <- gsub("-", "_", sub("^--", "", args[i]))
        if (!key %in% names(defaults) || i == length(args))
            stop("Unknown or incomplete argument: ", args[i],
                 ". Allowed arguments: ",
                 paste0("--", names(defaults), collapse = ", "))
        val <- args[i + 1]
        out[[key]] <- if (is.numeric(defaults[[key]])) as.numeric(val) else val
        i <- i + 2
    }
    out
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 || startsWith(args[1], "--"))
    stop("Usage: Rscript episewer.R <config.yml> [--tree_id N] ",
         "[--data_dir DIR] [--sequencing_delay N]")
opt <- c(read_config(args[1], config_defaults, names(cli_defaults)),
         parse_cli(args[-1], cli_defaults))
message("Read config from ", args[1])

# ---- data -----------------------------------------------------------------
tree_dir = paste0(opt$data_dir, "/tree_", opt$tree_id)
message(paste0("Working in ", tree_dir))

measurements <- fread(paste0(tree_dir, "/measurements.csv"))
flows        <- fread(paste0(tree_dir, "/flows.csv"))
measurements[, date := as.Date(date)]
measurements_sparse <- measurements[seq(3, nrow(measurements), by = 3), ]
flows[, date := as.Date(date)]

if(opt$breakpoints == ""){
    break_list <- seq(min(flows$date) + opt$sequencing_delay, max(flows$date), length.out = opt$num_segments+1)
    breakpoints <- paste(break_list, collapse = ",")
} else{
    breakpoints <- opt$breakpoints
}


# ---- write the parameters used for this run ---------------------------------
# Written before the model is fitted, so the file exists even if the fit fails.
# Note: the file sits in data_dir (not tree_dir), so it is overwritten when
# another tree_id is analysed with the same data_dir.
write_parameters <- function(file, opt, config_file, breakpoints) {
    fmt <- function(x) format(x, scientific = FALSE, trim = TRUE)
    lines <- c(
        "EpiSewer run parameters",
        paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
        paste0("Config file: ", config_file),
        paste0("EpiSewer version: ", as.character(packageVersion("EpiSewer"))),
        "",
        "[Command-line parameters]",
        paste0("data_dir = ", opt$data_dir),
        paste0("tree_dir = ", tree_dir),
        paste0("tree_id = ", fmt(opt$tree_id)),
        paste0("sequencing_delay = ", fmt(opt$sequencing_delay)),
        "",
        "[Config parameters]",
        paste0("shed_mean = ", fmt(opt$shed_mean)),
        paste0("shed_sd = ", fmt(opt$shed_sd)),
        paste0("gen_mean = ", fmt(opt$gen_mean)),
        paste0("gen_sd = ", fmt(opt$gen_sd)),
        paste0("incub_mean = ", fmt(opt$incub_mean)),
        paste0("incub_sd = ", fmt(opt$incub_sd)),
        paste0("load_per_case = ", format(opt$load_per_case, scientific = TRUE)),
        paste0("iter_warmup = ", fmt(opt$iter_warmup)),
        paste0("iter_sampling = ", fmt(opt$iter_sampling)),
        paste0("chains = ", fmt(opt$chains)),
        paste0("cores = ", fmt(opt$cores)),
        paste0("seed = ", fmt(opt$seed)),
        paste0("num_segments = ", fmt(opt$num_segments)),
        paste0("breakpoints (as given in config) = ",
               if (opt$breakpoints == "") "<empty>" else opt$breakpoints),
        "",
        "[Derived parameters]",
        paste0("breakpoints (used) = ", breakpoints),
        paste0("generation_dist: discrete shifted gamma, mean = ", fmt(opt$gen_mean),
               ", sd = ", fmt(opt$gen_sd)),
        paste0("shedding_dist: discrete gamma, shape = ",
               fmt((opt$shed_mean / opt$shed_sd)^2),
               ", scale = ", fmt(opt$shed_sd^2 / opt$shed_mean)),
        paste0("incubation_dist: discrete gamma, shape = ",
               fmt((opt$incub_mean / opt$incub_sd)^2),
               ", scale = ", fmt(opt$incub_sd^2 / opt$incub_mean)),
        "",
        "[Fixed settings in this script]",
        "shedding_reference = symptom_onset",
        "summary_intervals = 0.5, 0.95",
        paste0("samples_ndraws = ", fmt(opt$chains * opt$iter_sampling))
    )
    writeLines(lines, file)
}
 
params_file <- file.path(opt$data_dir, "parameters_episewer.txt")
write_parameters(params_file, opt, args[1], breakpoints)
message("Parameters written to ", params_file)

ww_data <- sewer_data(measurements = measurements, flows = flows)

# ---- assumptions (gamma shape = (mean/sd)^2, scale = sd^2/mean) ----------
ww_assumptions <- sewer_assumptions(
    generation_dist = get_discrete_gamma_shifted(gamma_mean = opt$gen_mean,
                                                 gamma_sd   = opt$gen_sd),
    shedding_dist   = get_discrete_gamma(gamma_shape = (opt$shed_mean / opt$shed_sd)^2,
                                         gamma_scale = opt$shed_sd^2 / opt$shed_mean),
    shedding_reference = "symptom_onset", #infection
    incubation_dist = get_discrete_gamma(gamma_shape = (opt$incub_mean / opt$incub_sd)^2, 
                                         gamma_scale = opt$incub_sd^2 / opt$incub_mean),
    load_per_case = opt$load_per_case
)

# ---- fit ------------------------------------------------------------------
options(mc.cores = opt$cores)
ww_result <- EpiSewer(
    data = ww_data,
    assumptions = ww_assumptions,
    shedding = model_shedding(load_per_case = load_per_case_assume()),
    fit_opts = set_fit_opts(sampler = sampler_stan_mcmc(
        iter_warmup = opt$iter_warmup, iter_sampling = opt$iter_sampling,
        chains = opt$chains, seed = opt$seed)),
    results_opts = set_results_opts(
        fitted = TRUE,
        summary_intervals = c(0.5, 0.95),
        samples_ndraws = opt$chains * opt$iter_sampling)
)

# ---- outputs ----------------------------------------------------------------
# saveRDS(ww_result, file.path(tree_dir, "episewer_result.rds"))

ggsave(file.path(tree_dir, "episewer_R_estimate.png"), plot_R(ww_result),
        width = 8, height = 4, dpi = 150)
ggsave(file.path(tree_dir, "episewer_growth_rate.png"), plot_growth_rate(ww_result),
        width = 8, height = 4, dpi = 150)
ggsave(file.path(tree_dir, "episewer_concentration_fit.png"),
        plot_concentration(ww_result, measurements = measurements),
        width = 8, height = 4, dpi = 150)

# Daily summary tables (mean/median/quantiles per date) as provided by
# EpiSewer directly - these are documented and safe to rely on.
fwrite(ww_result$summary$R, file.path(tree_dir, "episewer_R_daily_summary.csv"))
fwrite(ww_result$summary$growth_rate,
        file.path(tree_dir, "episewer_growth_rate_daily_summary.csv"))

# ---------------------------------------------------------------------------
# Per-section posterior distributions
# ---------------------------------------------------------------------------
# NOTE ON UNCERTAINTY: EpiSewer's result object exposes a `summary$samples`
# entry (see names(ww_result$summary)), which should hold posterior draws in
# long format, but its exact column names are not confirmed here against the
# package documentation. The block below inspects it once and stops with a
# clear message if the expected columns are not found, so you can adapt the
# column names below (search for "ADAPT HERE") instead of getting silently
# wrong results.
summarise_sections <- function(samples, variable_name, breakpoints,
                                section_labels = NULL, only_estimate = TRUE,
                                add_whole_period = TRUE) {
    if (!variable_name %in% names(samples))
        stop("Column '", variable_name, "' not found in ww_result$summary$samples. ",
                "Available columns: ", paste(names(samples), collapse = ", "))
    
    s <- samples[, c("date", ".draw", "type", variable_name), with = FALSE]
    setnames(s, variable_name, "value")
    s[, date := as.Date(date)]
    if (only_estimate && "type" %in% names(s)) s <- s[type == "estimate"]
    
    bp <- as.Date(breakpoints)
    n_sections <- length(bp) - 1
    if (is.null(section_labels))
        section_labels <- paste0("section_", seq_len(n_sections))
    
    # intervals to summarise: every section, plus (optionally) the whole period
    # from the first breakpoint to the last breakpoint
    starts <- bp[-length(bp)]
    ends   <- bp[-1]
    labels <- section_labels
    if (add_whole_period) {
        starts <- c(starts, bp[1])
        ends   <- c(ends,   bp[length(bp)])
        labels <- c(labels, "whole_period")
    }
    n_intervals <- length(labels)
    
    draws_list <- vector("list", n_intervals)
    summary_list <- vector("list", n_intervals)
    for (i in seq_len(n_intervals)) {
        sec <- s[date >= starts[i] & date < ends[i] & !is.na(value)]
        if (nrow(sec) == 0) {
            warning("No non-NA draws found for '", variable_name, "' in '",
                    labels[i], "' (", starts[i], " to ", ends[i],
                    "). This date range may fall entirely in the seeding period ",
                    "or outside the data. Skipping.")
            next
        }
        # one value per posterior draw: mean over the days in this interval
        per_draw <- sec[, .(value = mean(value)), by = .draw]
        per_draw[, `:=`(section = labels[i], variable = variable_name,
                        start_date = starts[i], end_date = ends[i])]
        draws_list[[i]] <- per_draw
        
        q <- quantile(per_draw$value, c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE)
        summary_list[[i]] <- data.table(
            section = labels[i], variable = variable_name,
            start_date = starts[i], end_date = ends[i],
            mean = mean(per_draw$value), sd = sd(per_draw$value),
            q2.5 = q[1], q25 = q[2], median = q[3], q75 = q[4], q97.5 = q[5],
            n_draws = nrow(per_draw)
        )
    }
    list(draws = rbindlist(draws_list, fill = TRUE),
            summary = rbindlist(summary_list, fill = TRUE))
}

breakpoints <- strsplit(breakpoints, ",")[[1]]
if (length(breakpoints) < 2)
    stop("'breakpoints' needs at least 2 dates to define 1 section.")

samples <- ww_result$summary$samples

results <- lapply(c("R", "growth_rate"), function(v)
    summarise_sections(samples, v, breakpoints, section_labels = NULL))
names(results) <- c("R", "growth_rate")

all_summary <- rbindlist(lapply(results, `[[`, "summary"), fill = TRUE)
all_draws   <- rbindlist(lapply(results, `[[`, "draws"), fill = TRUE)

fwrite(all_summary, file.path(tree_dir, "episewer_summary.csv"))
# fwrite(all_draws, file.path(tree_dir, "episewer_section_posterior_draws.csv"))
message("Section-level summary written to episewer_section_summary.csv")
# message("Full per-draw distributions written to episewer_section_posterior_draws.csv")


message("EpiSewer finished. Results in: ", tree_dir)