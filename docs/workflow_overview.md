# Workflow overview

This repository is organized into three workflow layers.

## 1. Upstream metagenome processing

Scripts:

- `scripts/01_metagenome_to_mags.sh`
- `scripts/02_annotate_ARG_VFDB.sh`
- `scripts/03_MAG_abundance_coverm.sh`

Purpose:

- assemble paired-end metagenomes
- map reads back to assemblies
- reconstruct MAGs
- dereplicate MAGs
- assess MAG quality
- assign taxonomy
- predict ORFs
- annotate ARGs and virulence-associated genes
- estimate MAG abundance

These scripts are templates for Linux/HPC environments and require real input FASTQ files and local external databases.

## 2. Matrix construction

Scripts:

- `scripts/04_parse_abricate_summary.py`
- `scripts/05_merge_abricate_matrices.py`
- `scripts/06_merge_MAG_abundance_taxonomy.py`
- `scripts/07_prepare_analysis_inputs.R`

Purpose:

- parse annotation outputs
- construct MAG-level feature matrices
- merge abundance, annotation and taxonomy tables
- generate sample-level input matrices for downstream analyses

## 3. Downstream analyses

Scripts:

- `scripts/08_phylogenetic_diversity_PD_GAM.R`
- `scripts/09_phylogenetic_structure_GAM.R`
- `scripts/10_environmental_screening_RF_XGBoost_SHAP.R`
- `scripts/11_phylofactor_hotspots.R`
- `scripts/12_hotspot_overlap_taxonomy.R`
- `scripts/13_make_analysis_panels.R`

Purpose:

- calculate phylogenetic diversity and structure metrics
- fit GAMs and segmented regression models
- screen environmental predictors using RF and XGBoost/SHAP
- fit multivariable mixed-effect GAMMs
- identify phylogenetic hotspot clades
- quantify overlap between ARG and virulence-associated hotspots
- export source-data and analytical-panel tables

## Notes on figure generation

The scripts generate analytical outputs and core plotting panels. Final figure composition, panel labels and graphical refinements may have been performed in Adobe Illustrator or equivalent vector-graphics software.
