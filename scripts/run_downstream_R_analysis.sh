#!/usr/bin/env bash
set -euo pipefail

# Example downstream analysis execution order. Edit input paths in each script or pass command-line options.

Rscript scripts/07_prepare_analysis_inputs.R \
  --mag_abundance data/MAG_TPM.tsv \
  --mag_arg data/MAG_ARG_annotations.tsv \
  --mag_vf data/MAG_VFDB_annotations.tsv \
  --env data/environmental_metadata.tsv

Rscript scripts/08_phylogenetic_diversity_PD_GAM.R \
  --tree data/unified_MAG_tree.nwk

Rscript scripts/09_phylogenetic_structure_GAM.R

Rscript scripts/10_environmental_screening_RF_XGBoost_SHAP.R \
  --endpoint ARG_load \
  --outdir results/environmental_models

Rscript scripts/10_environmental_screening_RF_XGBoost_SHAP.R \
  --endpoint VIR_load \
  --outdir results/environmental_models

Rscript scripts/11_phylofactor_hotspots.R \
  --tree data/unified_MAG_tree.nwk

Rscript scripts/12_hotspot_overlap_taxonomy.R \
  --taxonomy data/MAG_taxonomy.tsv

Rscript scripts/13_make_analysis_panels.R
