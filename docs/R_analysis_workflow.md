# Downstream R analysis workflow

The downstream R scripts provide a cleaned and modular version of the statistical analyses used in the manuscript. They are intended to be run after the processed MAG catalogue, annotation matrices, abundance matrices, metadata and phylogenetic tree have been generated or downloaded from the associated data repository.

## Script overview

- `07_prepare_analysis_inputs.R`: prepares sample-level and MAG-level derived matrices for downstream analyses.
- `08_phylogenetic_diversity_PD_GAM.R`: calculates Faith's phylogenetic diversity and fits PD-burden GAMs, derivative-based summaries and segmented-regression sensitivity analyses.
- `09_phylogenetic_structure_GAM.R`: calculates and models deep- and tip-level phylogenetic structure using MPD/MNTD standardized effect sizes and GAM comparisons.
- `10_environmental_screening_RF_XGBoost_SHAP.R`: performs RF and XGBoost/SHAP stability screening, VIF filtering, multivariable mixed-effect GAMMs, effect contrasts and residual Moran's I tests.
- `11_phylofactor_hotspots.R`: performs phylogenetic factorisation, hotspot assignment, feature enrichment and PERMANOVA summaries.
- `12_hotspot_overlap_taxonomy.R`: tests ARG-virulence hotspot overlaps and phylum-level taxonomic enrichment.
- `13_make_analysis_panels.R`: exports source-data tables and simplified analytical panels used for final figure assembly.

## Inputs

The expected inputs are listed in `config/analysis_paths_template.yml`. Large data files are not included in this repository and should be obtained from the data repositories described in the manuscript.

## Outputs

Outputs include analysis-ready matrices, statistical summary tables, source-data files and analytical panels. Final journal-ready composite figures may include manual layout and annotation in graphics software.
