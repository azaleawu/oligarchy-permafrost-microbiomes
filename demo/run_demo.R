#!/usr/bin/env Rscript

# Lightweight demo using simulated data.
# This script illustrates table structures and selected downstream calculations.
# It does not reproduce the manuscript results.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
})

base_dir <- "demo"
out_dir <- "results/demo_outputs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

presence <- fread(file.path(base_dir, "demo_mag_presence.tsv"))
arg <- fread(file.path(base_dir, "demo_arg_class_matrix.tsv"))
vf <- fread(file.path(base_dir, "demo_vf_module_matrix.tsv"))
env <- fread(file.path(base_dir, "demo_environmental_metadata.tsv"))

arg_total <- arg |>
  mutate(ARG_total = rowSums(across(-mag_id), na.rm = TRUE)) |>
  select(mag_id, ARG_total)

vf_total <- vf |>
  mutate(VF_total = rowSums(across(-mag_id), na.rm = TRUE)) |>
  select(mag_id, VF_total)

pres_long <- presence |>
  pivot_longer(-sample_id, names_to = "mag_id", values_to = "present") |>
  mutate(present = as.numeric(present))

sample_burdens <- pres_long |>
  left_join(arg_total, by = "mag_id") |>
  left_join(vf_total, by = "mag_id") |>
  group_by(sample_id) |>
  summarise(
    MAG_richness = sum(present > 0),
    ARG_burden = sum(present * ARG_total, na.rm = TRUE),
    VF_burden = sum(present * VF_total, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(env, by = "sample_id") |>
  mutate(ARG_log = log1p(ARG_burden), VF_log = log1p(VF_burden))

fwrite(sample_burdens, file.path(out_dir, "demo_sample_burdens.tsv"), sep = "\t")

model_summary <- list()
model_summary[["ARG_vs_GST"]] <- summary(lm(ARG_log ~ GST, data = sample_burdens))$coefficients
model_summary[["VF_vs_pH"]] <- summary(lm(VF_log ~ pH, data = sample_burdens))$coefficients

model_df <- rbind(
  data.frame(model = "ARG_log ~ GST", term = rownames(model_summary[["ARG_vs_GST"]]), model_summary[["ARG_vs_GST"]], row.names = NULL),
  data.frame(model = "VF_log ~ pH", term = rownames(model_summary[["VF_vs_pH"]]), model_summary[["VF_vs_pH"]], row.names = NULL)
)

fwrite(model_df, file.path(out_dir, "demo_model_summary.tsv"), sep = "\t")

writeLines(c(
  "Demo completed successfully.",
  "The simulated dataset illustrates input table structures and selected downstream calculations.",
  "It is not intended to reproduce the manuscript results."
), con = file.path(out_dir, "demo_readme.txt"))

message("Demo complete. Outputs written to: ", out_dir)
