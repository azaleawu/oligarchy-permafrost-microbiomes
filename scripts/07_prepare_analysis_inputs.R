#!/usr/bin/env Rscript

# Prepare downstream analysis matrices from processed MAG, annotation and metadata tables.
# This script is intentionally generic; edit config/analysis_paths_template.yml before use.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
config_file <- ifelse(length(args) >= 1, args[1], "config/analysis_paths_template.yml")
cfg <- yaml::read_yaml(config_file)

out_dir <- cfg$output$tables
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

presence <- fread(cfg$input$mag_presence)
arg_matrix <- fread(cfg$input$arg_matrix)
vf_matrix <- fread(cfg$input$vf_matrix)
sample_meta <- fread(cfg$input$sample_metadata)
mag_meta <- fread(cfg$input$mag_metadata)

stopifnot("sample_id" %in% names(presence))
stopifnot("mag_id" %in% names(arg_matrix))
stopifnot("mag_id" %in% names(vf_matrix))

# Convert sample x MAG presence table to long format.
pres_long <- presence |>
  pivot_longer(-sample_id, names_to = "mag_id", values_to = "present") |>
  mutate(present = as.numeric(present))

arg_burden <- arg_matrix |>
  mutate(arg_total = rowSums(across(-mag_id), na.rm = TRUE)) |>
  select(mag_id, arg_total)

vf_burden <- vf_matrix |>
  mutate(vf_total = rowSums(across(-mag_id), na.rm = TRUE)) |>
  select(mag_id, vf_total)

sample_burden <- pres_long |>
  left_join(arg_burden, by = "mag_id") |>
  left_join(vf_burden, by = "mag_id") |>
  group_by(sample_id) |>
  summarise(
    n_mags = sum(present > 0, na.rm = TRUE),
    ARG_burden = sum(present * replace_na(arg_total, 0), na.rm = TRUE),
    VF_burden = sum(present * replace_na(vf_total, 0), na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(sample_meta, by = "sample_id")

fwrite(sample_burden, file.path(out_dir, "sample_level_burden_summary.tsv"), sep = "\t")
fwrite(mag_meta, file.path(out_dir, "mag_metadata_used.tsv"), sep = "\t")
fwrite(arg_matrix, file.path(out_dir, "mag_by_arg_class_used.tsv"), sep = "\t")
fwrite(vf_matrix, file.path(out_dir, "mag_by_vf_module_used.tsv"), sep = "\t")

message("Prepared analysis input summaries in: ", out_dir)
