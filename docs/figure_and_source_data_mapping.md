# Figure and source-data mapping

This file maps manuscript-level outputs to the scripts that generate the corresponding analytical results or source data.

| Manuscript output | Main code modules |
|---|---|
| MAG catalogue quality and taxonomy | `01_metagenome_to_mags.sh`, `06_merge_MAG_abundance_taxonomy.py`, `07_prepare_analysis_inputs.R` |
| ARG and virulence-associated gene annotation summaries | `02_annotate_ARG_VFDB.sh`, `04_parse_abricate_summary.py`, `05_merge_abricate_matrices.py`, `07_prepare_analysis_inputs.R` |
| Faith's PD and burden/prevalence relationships | `08_phylogenetic_diversity_PD_GAM.R` |
| MPD/MNTD phylogenetic structure analyses | `09_phylogenetic_structure_GAM.R` |
| Environmental driver screening | `10_environmental_screening_RF_XGBoost_SHAP.R` |
| Multivariable environmental GAMMs | `10_environmental_screening_RF_XGBoost_SHAP.R` |
| Phylogenetic factorisation and hotspot membership | `11_phylofactor_hotspots.R` |
| Hotspot overlap and taxonomic enrichment | `12_hotspot_overlap_taxonomy.R` |
| Analytical panels and source-data exports | `13_make_analysis_panels.R` |

The code is organized by analysis module rather than by final figure number. Final composite figures may have been assembled and annotated in Adobe Illustrator. This repository provides the code and data-processing logic behind the statistical outputs, analytical panels and Source Data.
