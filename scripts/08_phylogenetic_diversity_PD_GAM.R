#!/usr/bin/env Rscript

# Phylogenetic diversity, burden/prevalence GAMs and segmented-regression analyses.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ape)
  library(picante)
  library(mgcv)
  library(segmented)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
config_file <- ifelse(length(args) >= 1, args[1], "config/analysis_paths_template.yml")
cfg <- yaml::read_yaml(config_file)
`%||%` <- function(a, b) if (!is.null(a)) a else b
set.seed(cfg$analysis$random_seed %||% 2025)
out_dir <- cfg$output$tables
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

presence <- fread(cfg$input$mag_presence) |> as.data.frame()
rownames(presence) <- presence$sample_id
presence$sample_id <- NULL
presence[] <- lapply(presence, as.numeric)

sample_summary <- fread(file.path(out_dir, "sample_level_burden_summary.tsv"))
tree <- read.tree(cfg$input$mag_tree)

common_tips <- intersect(colnames(presence), tree$tip.label)
presence <- presence[, common_tips, drop = FALSE]
tree <- keep.tip(tree, common_tips)

pd_tbl <- pd(presence, tree, include.root = FALSE)
pd_tbl$sample_id <- rownames(pd_tbl)

analysis_df <- sample_summary |>
  left_join(pd_tbl, by = "sample_id") |>
  mutate(
    ARG_log = log1p(ARG_burden),
    VF_log = log1p(VF_burden)
  )

fit_gam <- function(response, data) {
  formula <- as.formula(paste0(response, " ~ s(PD, k = 6, bs = 'cs')"))
  gam(formula, data = data, method = "REML")
}

arg_gam <- fit_gam("ARG_log", analysis_df)
vf_gam <- fit_gam("VF_log", analysis_df)

summarise_gam <- function(model, endpoint) {
  s <- summary(model)
  data.frame(
    endpoint = endpoint,
    adj_r2 = s$r.sq,
    deviance_explained = s$dev.expl,
    edf = s$s.table[1, "edf"],
    F = s$s.table[1, "F"],
    p_value = s$s.table[1, "p-value"]
  )
}

gam_summary <- bind_rows(
  summarise_gam(arg_gam, "ARG"),
  summarise_gam(vf_gam, "Virulence")
)

fwrite(analysis_df, file.path(out_dir, "phylogenetic_diversity_sample_table.tsv"), sep = "\t")
fwrite(gam_summary, file.path(out_dir, "PD_GAM_summary.tsv"), sep = "\t")

# Segmented regression as sensitivity analysis.
segmented_fit <- function(response, data) {
  base <- lm(as.formula(paste0(response, " ~ PD")), data = data)
  tryCatch({
    segmented(base, seg.Z = ~ PD)
  }, error = function(e) NULL)
}

seg_arg <- segmented_fit("ARG_log", analysis_df)
seg_vf <- segmented_fit("VF_log", analysis_df)

extract_breakpoint <- function(model, endpoint) {
  if (is.null(model)) return(data.frame(endpoint = endpoint, breakpoint = NA_real_))
  data.frame(endpoint = endpoint, breakpoint = as.numeric(model$psi[1, "Est."]))
}

seg_summary <- bind_rows(
  extract_breakpoint(seg_arg, "ARG"),
  extract_breakpoint(seg_vf, "Virulence")
)

fwrite(seg_summary, file.path(out_dir, "PD_segmented_breakpoints.tsv"), sep = "\t")
message("PD and GAM analyses complete.")
