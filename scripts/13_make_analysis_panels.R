#!/usr/bin/env Rscript

# Export source-data tables and simple analytical panels for manuscript figure assembly.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
config_file <- ifelse(length(args) >= 1, args[1], "config/analysis_paths_template.yml")
cfg <- yaml::read_yaml(config_file)

table_dir <- cfg$output$tables
panel_dir <- cfg$output$panels
figure_dir <- cfg$output$figures

dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

pd_file <- file.path(table_dir, "phylogenetic_diversity_sample_table.tsv")
if (file.exists(pd_file)) {
  df <- fread(pd_file)
  fwrite(df, file.path(panel_dir, "source_data_phylogenetic_diversity.tsv"), sep = "\t")
  if (all(c("PD", "ARG_log") %in% names(df))) {
    p <- ggplot(df, aes(PD, ARG_log)) + geom_point() + geom_smooth(method = "gam", formula = y ~ s(x, k = 6))
    ggsave(file.path(figure_dir, "panel_PD_ARG.pdf"), p, width = 4, height = 3)
  }
  if (all(c("PD", "VF_log") %in% names(df))) {
    p <- ggplot(df, aes(PD, VF_log)) + geom_point() + geom_smooth(method = "gam", formula = y ~ s(x, k = 6))
    ggsave(file.path(figure_dir, "panel_PD_Virulence.pdf"), p, width = 4, height = 3)
  }
}

message("Analytical panel and source-data export complete.")
