#!/usr/bin/env bash
set -euo pipefail

CONFIG=${1:-config/analysis_paths_template.yml}

Rscript scripts/07_prepare_analysis_inputs.R "$CONFIG"
Rscript scripts/08_phylogenetic_diversity_PD_GAM.R "$CONFIG"
Rscript scripts/09_phylogenetic_structure_GAM.R "$CONFIG" || true
Rscript scripts/10_environmental_screening_RF_XGBoost_SHAP.R "$CONFIG" || true
Rscript scripts/11_phylofactor_hotspots.R "$CONFIG" || true
Rscript scripts/12_hotspot_overlap_taxonomy.R "$CONFIG" || true
Rscript scripts/13_make_analysis_panels.R "$CONFIG" || true
