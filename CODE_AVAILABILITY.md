# Code availability

This repository contains cleaned workflow scripts and downstream R analysis scripts supporting the genome-resolved analysis described in the manuscript.

The repository includes code for:

- metagenome assembly, read mapping and MAG reconstruction;
- MAG dereplication, quality assessment and taxonomic assignment;
- ARG annotation using RGI/CARD;
- virulence-associated gene annotation using ABRicate/VFDB;
- MAG abundance profiling;
- construction of analysis-ready matrices;
- phylogenetic diversity and phylogenetic-structure analyses;
- environmental driver screening using random forest and XGBoost/SHAP;
- multivariable mixed-effect GAMMs;
- phylogenetic factorisation and hotspot analyses;
- hotspot overlap and taxonomic enrichment analyses;
- generation of source-data tables and analytical plotting panels.

Large input datasets and derived data tables are available from NCBI SRA, GenBank, Zenodo and/or Source Data files as described in the manuscript. Raw reads, database files and large intermediate files are not included in this repository.

The scripts are generalized and use configurable input paths. Local adaptation may be needed for software installation paths, database locations and computing-cluster environments.
