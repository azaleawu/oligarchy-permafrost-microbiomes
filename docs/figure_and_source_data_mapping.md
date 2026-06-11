# Analysis-module mapping

This repository is organized by analysis module rather than by final figure number. The scripts generate analytical outputs, source-data tables and plotting panels that can be used for final figure assembly.

## Main analysis modules

- `07_prepare_analysis_inputs.R`: sample-level and MAG-level derived matrices used by downstream analyses.
- `08_phylogenetic_diversity_PD_GAM.R`: Faith's phylogenetic diversity, PD-burden GAMs, derivative-based threshold summaries and segmented-regression sensitivity outputs.
- `09_phylogenetic_structure_GAM.R`: MPD/MNTD standardized effect sizes, tensor-product GAM comparisons and phylogenetic-structure response surfaces.
- `10_environmental_screening_RF_XGBoost_SHAP.R`: RF/XGBoost-SHAP screening, VIF filtering, multivariable mixed-effect GAMMs, effect contrasts and residual Moran's I tests.
- `11_phylofactor_hotspots.R`: phylogenetic factorisation, hotspot membership, feature enrichment and PERMANOVA summaries.
- `12_hotspot_overlap_taxonomy.R`: ARG-virulence hotspot overlaps, Fisher exact tests, Jaccard indices and taxonomic enrichment.
- `13_make_analysis_panels.R`: export of source-data tables and simplified analytical plotting panels.

## Figure assembly note

Some final composite figures may include manual layout, annotations, panel labels or graphical elements added in Adobe Illustrator or equivalent graphics software. This repository provides the computational analyses underlying the quantitative panels and tables, not a one-command reproduction of final journal-ready figure layouts.
