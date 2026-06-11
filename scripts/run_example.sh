#!/usr/bin/env bash
set -euo pipefail

# Example execution order.
# Copy/edit config/config_template.env to config/config.env before running.
# Copy/edit config/sample_table_template.tsv and config/mag_table_template.tsv with local paths.

bash scripts/01_metagenome_to_mags.sh config/sample_table_template.tsv config/config.env
bash scripts/02_annotate_ARG_VFDB.sh config/mag_table_template.tsv config/config.env
bash scripts/03_MAG_abundance_coverm.sh config/mag_table_template.tsv config/config.env

# Example conversion of ABRicate/VFDB summary to a count matrix:
# python scripts/04_parse_abricate_summary.py \
#   --input results/processed/annotations/SAMPLE_001/vfdb/SAMPLE_001.vfdb.summary.tsv \
#   --output results/processed/annotations/SAMPLE_001/vfdb/SAMPLE_001.vfdb.counts.csv
