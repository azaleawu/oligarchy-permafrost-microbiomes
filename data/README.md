# Data directory

This repository does not store raw metagenomic reads, large assembly outputs, BAM files, external databases or full intermediate output directories.

Expected input data for full analyses include:

- sample metadata
- environmental metadata
- MAG metadata
- MAG abundance matrix
- sample × MAG presence matrix
- MAG phylogeny in Newick format
- MAG × ARG class matrix
- MAG × virulence-associated module matrix

These files are available from the data repositories and Source Data described in the manuscript.

Use `config/analysis_paths_template.yml` to point the downstream R scripts to local copies of these files.
