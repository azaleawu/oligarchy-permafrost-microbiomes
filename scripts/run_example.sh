#!/usr/bin/env bash
set -euo pipefail

# Example commands only. Edit configuration files before running on real data.

bash scripts/01_metagenome_to_mags.sh config/config_template.env
bash scripts/02_annotate_ARG_VFDB.sh config/config_template.env
bash scripts/03_MAG_abundance_coverm.sh config/config_template.env
bash scripts/run_downstream_R_analysis.sh config/analysis_paths_template.yml
