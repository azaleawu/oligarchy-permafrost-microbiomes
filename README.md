# This repository contains cleaned workflow scripts and downstream analysis code supporting the manuscript

The repository provides a pipeline for metagenome processing, MAG reconstruction, gene annotation, abundance profiling, matrix construction, phylogenetic analyses, environmental modelling and hotspot analyses.

Large input datasets and intermediate files are **not** included in this repository. Raw metagenomic reads, derived MAG catalogues, annotation matrices, environmental metadata and Source Data are available from the repositories described in the manuscript.

---

## Contents

- [Overview](#overview)
- [Repository contents](#repository-contents)
- [System requirements](#system-requirements)
- [Installation guide](#installation-guide)
- [Demo](#demo)
- [Instructions for use](#instructions-for-use)
- [Expected outputs](#expected-outputs)
- [Figure and source-data mapping](#figure-and-source-data-mapping)
- [Data availability](#data-availability)
- [Code availability](#code-availability)
- [License](#license)
- [Citation](#citation)

---

## Overview

The workflow includes:

1. **Metagenome processing and MAG reconstruction**
   - read trimming
   - metagenome assembly
   - read mapping
   - MAG binning
   - dereplication
   - MAG quality assessment
   - taxonomic assignment

2. **Gene annotation and abundance profiling**
   - ORF prediction
   - ARG annotation using RGI/CARD
   - virulence-associated gene annotation using ABRicate/VFDB
   - MAG abundance profiling using CoverM
   - annotation matrix construction

3. **Downstream statistical analyses**
   - sample-level matrix preparation
   - phylogenetic diversity and structure analyses
   - generalized additive models
   - segmented regression
   - environmental driver screening using random forest and XGBoost/SHAP
   - multivariable mixed-effect GAMMs
   - phylogenetic factorisation
   - hotspot overlap and enrichment analyses
   - source-data and analytical-panel export

The scripts are organized by **analysis module**, not by final figure number. Final composite figures may be assembled and annotated in Adobe Illustrator or equivalent vector-graphics software. The scripts generate the analytical tables, source data and core plotting panels that underlie the manuscript figures.

---

## Repository contents

| Path | Description |
|---|---|
| `README.md` | Main repository documentation. |
| `CODE_AVAILABILITY.md` | Code availability statement and access notes. |
| `LICENSE` | MIT License for code reuse. |
| `config/` | Template configuration files and input table templates. |
| `scripts/` | Shell, Python and R scripts for the analysis workflow. |
| `environment/` | Software versions and R package dependencies. |
| `demo/` | Small simulated dataset and a lightweight demonstration script. |
| `data/` | Documentation for external input data. No large datasets are stored here. |
| `results/` | Documentation for expected output directories. No manuscript-scale results are stored here. |
| `figures/` | Documentation for analytical panels and source-data outputs. |
| `docs/` | Additional workflow and checklist documentation. |

---

## System requirements

### Hardware requirements

The full metagenome assembly, mapping and MAG reconstruction workflow requires a Linux high-performance computing environment. Actual requirements depend on sample sequencing depth and assembly size, but typical runs require multiple CPU threads, substantial RAM and sufficient temporary storage.

Recommended resources for full-scale upstream processing:

- Operating system: Linux / HPC environment
- CPU: 16 or more threads recommended for assembly and mapping
- RAM: 128 GB or more recommended for large metagenome assemblies
- Storage: sufficient space for FASTQ, assemblies, BAM files and intermediate outputs

The downstream R analyses can be run on a standard workstation once processed matrices are available. For the complete study-scale matrices, we recommend:

- RAM: 16 GB or more
- CPU: 4 or more cores
- R: version 4.2 or higher

The small demo dataset provided in `demo/` can be run on a normal desktop computer.

### Software requirements

The workflow was developed for Linux/HPC environments. The downstream R scripts are platform-independent in principle, provided that required R packages are available.

Core command-line tools used by the upstream workflow include:

- Trimmomatic
- MEGAHIT
- Bowtie2
- Samtools
- VAMB
- dRep
- CheckM2
- GTDB-Tk
- CoverM
- Prodigal
- RGI with CARD
- ABRicate with VFDB

R package dependencies include:

- `ape`
- `picante`
- `phytools`
- `phangorn`
- `mgcv`
- `gamm4`
- `gratia`
- `segmented`
- `vegan`
- `data.table`
- `dplyr`
- `tidyr`
- `tibble`
- `ggplot2`
- `ComplexHeatmap`
- `ranger`
- `xgboost`
- `fastshap`
- `caret`
- `car`
- `spdep`
- `readr`
- `yaml`

Detailed version information is listed in:

```text
environment/software_versions.txt
environment/R_packages.txt
```

---

## Installation guide

### Command-line tools

Install the required metagenomic software using `conda`, `mamba`, modules on an HPC system, or the official installation instructions for each package. A template conda environment is provided in:

```text
environment/conda_environment.yml
```

Example installation using `mamba`:

```bash
mamba env create -f environment/conda_environment.yml
mamba activate permafrost_arg_vf_pipeline
```

Some databases, such as GTDB, CARD, VFDB and CheckM2 databases, must be downloaded separately according to the documentation of the corresponding software. These database files are not included in this repository.

Typical install time depends strongly on the computing environment and database availability. Installing the core command-line tools may take 30-120 minutes. Downloading large external databases may take several hours.

### R packages

Install required R packages from CRAN, Bioconductor and GitHub as needed. A package list is provided in:

```text
environment/R_packages.txt
```

A minimal installation command for common CRAN packages is:

```r
install.packages(c(
  "ape", "mgcv", "segmented", "vegan", "data.table", "dplyr", "tidyr",
  "tibble", "ggplot2", "ranger", "xgboost", "fastshap", "caret",
  "car", "spdep", "readr", "yaml"
))
```

Packages such as `picante`, `gamm4`, `gratia`, `ComplexHeatmap`, `phytools` and `phangorn` should be installed according to the current instructions of CRAN or Bioconductor.

Typical installation time for the R package environment is 30-60 minutes on a standard desktop computer, excluding compilation time and system-level dependencies.

---

## Demo

A small simulated dataset is provided in:

```text
demo/
```

The demo is designed to illustrate input table formats and downstream workflow usage. It is **not** intended to reproduce manuscript results.

To run the demo:

```bash
Rscript demo/run_demo.R
```

Expected runtime for the demo is less than one minute on a normal desktop computer with the required R packages installed.

Expected demo outputs are written to:

```text
results/demo_outputs/
```

The expected outputs include:

- `demo_sample_burdens.tsv`
- `demo_model_summary.tsv`
- `demo_readme.txt`

If optional phylogenetic packages are not available, the demo still reports matrix-derived summaries and a simple model summary.

---

## Instructions for use

### 1. Configure input paths

Copy and edit the template configuration files:

```text
config/config_template.env
config/analysis_paths_template.yml
config/sample_table_template.tsv
config/mag_table_template.tsv
```

Do not edit the scripts directly to hard-code personal paths. Use configuration files or command-line arguments.

### 2. Run upstream metagenome workflow

The upstream workflow scripts are provided as generic templates:

```text
scripts/01_metagenome_to_mags.sh
scripts/02_annotate_ARG_VFDB.sh
scripts/03_MAG_abundance_coverm.sh
```

Run them after editing `config/config_template.env` and sample/MAG tables. These scripts are designed for Linux/HPC environments.

### 3. Construct matrices

Use:

```text
scripts/04_parse_abricate_summary.py
scripts/05_merge_abricate_matrices.py
scripts/06_merge_MAG_abundance_taxonomy.py
scripts/07_prepare_analysis_inputs.R
```

These scripts generate downstream matrices used for phylogenetic, environmental and hotspot analyses.

### 4. Run downstream statistical analyses

Run the downstream scripts in numerical order:

```text
scripts/08_phylogenetic_diversity_PD_GAM.R
scripts/09_phylogenetic_structure_GAM.R
scripts/10_environmental_screening_RF_XGBoost_SHAP.R
scripts/11_phylofactor_hotspots.R
scripts/12_hotspot_overlap_taxonomy.R
scripts/13_make_analysis_panels.R
```

A wrapper script is provided:

```bash
bash scripts/run_downstream_R_analysis.sh
```

The wrapper assumes that paths in `config/analysis_paths_template.yml` have been edited to point to real input files.

---

## Expected outputs

Expected outputs include:

- MAG quality and taxonomy tables
- ARG and virulence-associated gene annotation matrices
- MAG abundance matrices
- sample-level burden/prevalence summaries
- Faith's PD, MPD and MNTD summaries
- GAM and segmented-regression summaries
- environmental driver screening tables
- GAMM effect-size summaries
- phylogenetic factorisation hotspot membership tables
- PERMANOVA and overlap summaries
- analytical panel tables for figure construction

Large manuscript-scale outputs are not stored in this repository. They are available through the data repositories described in the manuscript.

---

## Figure and source-data mapping

The repository is organized by analysis module. Approximate mapping to manuscript outputs is provided in:

```text
docs/figure_and_source_data_mapping.md
```

Final composite figures may have been assembled and annotated in Adobe Illustrator. This repository provides the code used to generate the analytical results, source data and core plotting panels.

---

## Data availability

This repository does not contain raw metagenomic reads, genome assemblies, BAM files, external databases or full intermediate output directories.

The data required to reproduce the full analysis are described in the manuscript Data availability statement, including:

- NCBI SRA accessions for raw metagenomic reads
- GenBank records for newly generated MAGs where applicable
- Zenodo records for processed MAG catalogues, annotations, matrices and metadata
- Source Data files submitted with the manuscript

---

## Code availability

Cleaned workflow scripts and downstream analysis code are available in this repository.

A versioned release of this repository will be archived in a DOI-minting repository, such as Zenodo, upon publication.

---

## License

This code is released under the MIT License. See:

```text
LICENSE
```

---

## Citation

If using this repository before formal publication, please cite the manuscript title and repository URL. After publication, please cite the published article and the archived code DOI.
