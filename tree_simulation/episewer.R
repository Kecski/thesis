#!/usr/bin/env Rscript
# Run EpiSewer on the CSV files written by run_simulation.py.
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
# breakpoints: comma-separated dates defining n sections. n+1 dates give n
# half-open sections [breakpoints[i], breakpoints[i+1]). Example for 3
# sections of unequal length: "2020-03-01,2020-04-01,2020-05-15,2020-07-01"
defaults <- list(data_dir = "data/wastewater", tree_id = 0,
                 shed_mean = 8, shed_sd = 5, gen_mean = 3, gen_sd = 2.4,
                 iter_warmup = 500, iter_sampling = 500, chains = 4,
                 cores = 4, seed = 42, num_segments = 1, sequencing_delay = 0,
                 breakpoints = "") # format for breakpoints: "YYYY-MM-DD,YYYY-MM-DD,..." 

parse_cli <- function(defaults) {
    args <- commandArgs(trailingOnly = TRUE)
    out <- defaults
    i <- 1
    while (i <= length(args)) {
        key <- gsub("-", "_", sub("^--", "", args[i]))
        if (!key %in% names(defaults) || i == length(args))
            stop("Unknown or incomplete argument: ", args[i])
        val <- args[i + 1]
        out[[key]] <- if (is.numeric(defaults[[key]]) && key != "breakpoints")
            as.numeric(val) else val
        i <- i + 2
    }
    out
}
opt <- parse_cli(defaults)

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


ww_data <- sewer_data(measurements = measurements, flows = flows)

# ---- assumptions (gamma shape = (mean/sd)^2, scale = sd^2/mean) ----------
ww_assumptions <- sewer_assumptions(
    generation_dist = get_discrete_gamma_shifted(gamma_mean = opt$gen_mean,
                                                    gamma_sd   = opt$gen_sd),
    shedding_dist   = get_discrete_gamma(gamma_shape = (opt$shed_mean / opt$shed_sd)^2,
                                            gamma_scale = opt$shed_sd^2 / opt$shed_mean),
    shedding_reference = "symptom_onset", #infection
    incubation_dist = get_discrete_gamma(gamma_shape = 8.5, gamma_scale = 0.4)
)

# ---- fit ------------------------------------------------------------------
options(mc.cores = opt$cores)
ww_result <- EpiSewer(
    data = ww_data,
    assumptions = ww_assumptions,
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
    stop("--breakpoints needs at least 2 dates to define 1 section.")

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


message("EpiSewer finished. Results in: ", opt$data_dir, "/tree_", opt$tree_id)


