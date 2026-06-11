#!/usr/bin/env Rscript

# Compare GAMs using MPD_Z and MNTD_Z as deep- and tip-level phylogenetic structure summaries.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(mgcv)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
config_file <- ifelse(length(args) >= 1, args[1], "config/analysis_paths_template.yml")
cfg <- yaml::read_yaml(config_file)
out_dir <- cfg$output$tables

input_file <- file.path(out_dir, "phylogenetic_diversity_sample_table.tsv")
df <- fread(input_file) |> as.data.frame()

required <- c("ARG_log", "VF_log", "MPD_Z", "MNTD_Z")
missing <- setdiff(required, names(df))
if (length(missing) > 0) {
  stop("Missing required columns: ", paste(missing, collapse = ", "))
}

fit_set <- function(response) {
  list(
    mpd = gam(as.formula(paste0(response, " ~ s(MPD_Z, k = 6, bs = 'cs')")), data = df, method = "REML"),
    mntd = gam(as.formula(paste0(response, " ~ s(MNTD_Z, k = 6, bs = 'cs')")), data = df, method = "REML"),
    additive = gam(as.formula(paste0(response, " ~ s(MPD_Z, k = 6, bs = 'cs') + s(MNTD_Z, k = 6, bs = 'cs')")), data = df, method = "REML"),
    tensor = gam(as.formula(paste0(response, " ~ te(MPD_Z, MNTD_Z, k = c(6, 6))")), data = df, method = "REML")
  )
}

summarise_model <- function(model, endpoint, model_name) {
  s <- summary(model)
  data.frame(
    endpoint = endpoint,
    model = model_name,
    AIC = AIC(model),
    adj_r2 = s$r.sq,
    deviance_explained = s$dev.expl
  )
}

arg_models <- fit_set("ARG_log")
vf_models <- fit_set("VF_log")

model_summary <- bind_rows(
  lapply(names(arg_models), function(n) summarise_model(arg_models[[n]], "ARG", n)),
  lapply(names(vf_models), function(n) summarise_model(vf_models[[n]], "Virulence", n))
)

fwrite(model_summary, file.path(out_dir, "MPD_MNTD_GAM_model_comparison.tsv"), sep = "\t")
message("MPD/MNTD GAM comparisons complete.")
