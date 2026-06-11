#!/usr/bin/env Rscript

# Overlap and enrichment analyses for ARG and virulence-associated hotspot memberships.

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

arg_file <- file.path(out_dir, "ARG_hotspot_membership.tsv")
vf_file <- file.path(out_dir, "Virulence_hotspot_membership.tsv")
mag_meta_file <- cfg$input$mag_metadata

if (!file.exists(arg_file) || !file.exists(vf_file)) {
  message("Hotspot membership files not found. Expected ARG_hotspot_membership.tsv and Virulence_hotspot_membership.tsv in output tables directory.")
  quit(save = "no", status = 0)
}

arg <- fread(arg_file)
vf <- fread(vf_file)
mag_meta <- fread(mag_meta_file)

jaccard_pair <- function(a_set, b_set) {
  inter <- length(intersect(a_set, b_set))
  union <- length(union(a_set, b_set))
  if (union == 0) return(NA_real_)
  inter / union
}

pairs <- expand.grid(ARG_hotspot = unique(arg$hotspot), VF_hotspot = unique(vf$hotspot), stringsAsFactors = FALSE)
overlap <- pairs |>
  rowwise() |>
  mutate(
    jaccard = jaccard_pair(arg$mag_id[arg$hotspot == ARG_hotspot], vf$mag_id[vf$hotspot == VF_hotspot]),
    overlap_n = length(intersect(arg$mag_id[arg$hotspot == ARG_hotspot], vf$mag_id[vf$hotspot == VF_hotspot]))
  ) |>
  ungroup()

fwrite(overlap, file.path(out_dir, "ARG_Virulence_hotspot_overlap.tsv"), sep = "\t")

if ("phylum" %in% names(mag_meta)) {
  merged <- arg |> left_join(mag_meta, by = "mag_id")
  enrich <- merged |>
    count(hotspot, phylum, name = "n") |>
    group_by(hotspot) |>
    mutate(fraction = n / sum(n)) |>
    ungroup()
  fwrite(enrich, file.path(out_dir, "ARG_hotspot_taxonomic_composition.tsv"), sep = "\t")
}

message("Hotspot overlap and taxonomy summaries complete.")
