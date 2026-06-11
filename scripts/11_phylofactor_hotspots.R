#!/usr/bin/env Rscript

# Phylogenetic factorisation and hotspot analysis.
# The script prepares MAG-level feature matrices and, when the phylofactor package
# is available, runs phylogenetic factorisation with the manuscript-level settings.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ape)
  library(vegan)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
config_file <- ifelse(length(args) >= 1, args[1], "config/analysis_paths_template.yml")
cfg <- yaml::read_yaml(config_file)
out_dir <- cfg$output$tables
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tree <- read.tree(cfg$input$mag_tree)
arg_matrix <- fread(cfg$input$arg_matrix) |> as.data.frame()
vf_matrix <- fread(cfg$input$vf_matrix) |> as.data.frame()

prepare_matrix <- function(feature_df) {
  rownames(feature_df) <- feature_df$mag_id
  feature_df$mag_id <- NULL
  feature_df <- as.matrix(feature_df)
  storage.mode(feature_df) <- "numeric"
  common <- intersect(rownames(feature_df), tree$tip.label)
  list(
    matrix = feature_df[common, , drop = FALSE],
    tree = keep.tip(tree, common)
  )
}

run_phylofactor_if_available <- function(feature_df, endpoint) {
  prepared <- prepare_matrix(feature_df)
  feature_matrix <- prepared$matrix
  pruned_tree <- prepared$tree

  fwrite(data.frame(mag_id = rownames(feature_matrix), feature_matrix),
         file.path(out_dir, paste0(endpoint, "_feature_matrix_for_phylofactor.tsv")), sep = "\t")

  if (!requireNamespace("phylofactor", quietly = TRUE)) {
    message("phylofactor package not installed. Exported input matrix for ", endpoint, ".")
    return(invisible(NULL))
  }

  # The following call reflects the manuscript settings. Depending on the installed
  # phylofactor version, users may need to adjust the namespace/function signature.
  pf <- tryCatch({
    phylofactor::PhyloFactor(
      data = feature_matrix,
      tree = pruned_tree,
      nfactors = 5,
      method = "max.var",
      choice = "var"
    )
  }, error = function(e) {
    message("PhyloFactor call failed for ", endpoint, ": ", conditionMessage(e))
    NULL
  })

  if (is.null(pf)) return(invisible(NULL))

  saveRDS(pf, file.path(out_dir, paste0(endpoint, "_phylofactor_model.rds")))
  message("Saved phylofactor model for ", endpoint)
  invisible(pf)
}

run_phylofactor_if_available(arg_matrix, "ARG")
run_phylofactor_if_available(vf_matrix, "Virulence")

message("Phylofactor analysis step complete.")
