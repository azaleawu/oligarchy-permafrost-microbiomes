# Genome-resolved analysis of ARGs and virulence-associated genes in permafrost metagenomes

This repository provides a generalized computational workflow for a genome-resolved metagenomic analysis of antibiotic resistance genes (ARGs) and virulence-associated genes in permafrost metagenomes.

The repository is designed to document the analysis logic, software, parameters and scripts used to generate the main derived matrices, statistical outputs, source-data tables and analytical plotting panels. Raw sequencing reads, large intermediate files and external databases are not included in this repository.

## What this repository contains

The workflow includes scripts for:

1. read-level processing, assembly and MAG reconstruction;
2. MAG dereplication, quality assessment and taxonomic assignment;
3. ORF prediction, ARG annotation using RGI/CARD and virulence-associated gene annotation using ABRicate/VFDB;
4. MAG abundance profiling and derived matrix construction;
5. phylogenetic diversity and phylogenetic-structure analyses;
6. environmental driver screening using random forest and XGBoost/SHAP;
7. multivariable mixed-effect GAMMs;
8. phylogenetic factorisation and hotspot analyses;
9. ARG–virulence hotspot overlap and taxonomic enrichment analyses;
10. generation of source-data tables and analytical plotting panels.

The scripts are organized by analysis module rather than by final figure number. Final composite figures may be assembled and annotated in graphics software; the scripts here generate the analytical results, source-data tables and plotting panels underlying those figures.

## Repository structure

```text
config/
  config_template.env                 # shell pipeline parameters and paths; edit locally before use
  sample_table_template.tsv           # paired-end read manifest template
  mag_table_template.tsv              # MAG/read manifest template for abundance and annotation workflows
  analysis_paths_template.yml         # downstream R analysis path template

scripts/
  01_metagenome_to_mags.sh            # assembly, mapping, VAMB binning, dRep, CheckM2 and GTDB-Tk
  02_annotate_ARG_VFDB.sh             # Prodigal, RGI/CARD and ABRicate/VFDB annotation
  03_MAG_abundance_coverm.sh          # MAG abundance profiling with Bowtie2/Samtools/CoverM
  04_parse_abricate_summary.py        # convert ABRicate summary strings to gene-count matrices
  05_merge_abricate_matrices.py       # merge gene-count matrices when needed
  06_merge_MAG_abundance_taxonomy.py  # merge MAG abundance, taxonomy and annotation summaries
  07_prepare_analysis_inputs.R        # construct matrices for downstream ecological analysis
  08_phylogenetic_diversity_PD_GAM.R  # Faith's PD, GAMs, PD50 and segmented-regression sensitivity
  09_phylogenetic_structure_GAM.R     # MPD/MNTD structure and tensor-product GAMs
  10_environmental_screening_RF_XGBoost_SHAP.R # RF/XGBoost-SHAP screening and multivariable GAMMs
  11_phylofactor_hotspots.R           # phylogenetic factorisation and hotspot enrichment
  12_hotspot_overlap_taxonomy.R       # ARG–virulence hotspot overlap and taxonomic enrichment
  13_make_analysis_panels.R           # source-data and analytical-panel exports
  run_example.sh                      # example shell-pipeline execution order
  run_downstream_R_analysis.sh        # example downstream R execution order

environment/
  software_versions.txt               # software and database versions reported for the workflow
  R_packages.txt                      # R packages used by the downstream scripts

data/
  README.md                           # expected external input data; large files are not stored here

results/
  README.md                           # expected output structure; generated files are not stored here

docs/
  R_analysis_workflow.md              # overview of downstream R analysis modules
  figure_and_source_data_mapping.md   # mapping between analysis modules and manuscript outputs
```

## Data availability

Large input datasets are not stored in this repository. Raw metagenomic reads are available from NCBI SRA as described in the manuscript. Genome assemblies, MAG catalogues, annotation tables, environmental metadata, abundance matrices and derived analysis matrices are provided through the public data repositories and source-data files described in the manuscript.

Intermediate files such as raw FASTQ files, SAM/BAM files, assemblies, binning directories and local copies of CARD, VFDB, GTDB and other databases are not included.

## Usage overview

1. Edit `config/config_template.env` and save it locally as `config/config.env`.
2. Fill `config/sample_table_template.tsv` with sample IDs and paired-end read paths.
3. Run the shell workflow scripts in `scripts/01*` to `scripts/06*` for assembly, MAG reconstruction, annotation and abundance profiling.
4. Edit `config/analysis_paths_template.yml` for downstream table locations.
5. Run `scripts/run_downstream_R_analysis.sh` or the individual R scripts for statistical analyses and source-data generation.

The scripts are intended as transparent workflow templates. Local adaptation may be required for cluster schedulers, module systems, conda environments and database paths.

## Reproducibility notes

- User-specific paths, batch logs, exploratory analyses and unrelated workflows are not included.
- Scripts use generic file names and configurable input paths.
- Large input and intermediate files are referenced through external repositories rather than stored in GitHub.
- Random seeds are set in stochastic downstream analyses where applicable.
- Final composite figures may include manual layout or annotation in graphics software; analytical panels and source-data tables are generated by the scripts.

## License

This repository is released under the MIT License. See `LICENSE`.
