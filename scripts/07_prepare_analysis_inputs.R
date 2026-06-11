#!/usr/bin/env Rscript
# Prepare derived MAG- and sample-level input tables for downstream analyses.
# This script starts from cleaned annotation/abundance tables and creates:
#   1) sample x MAG presence/absence matrix
#   2) sample-level ARG and virulence-associated gene prevalence/load summaries
#   3) MAG x ARG-class and MAG x virulence-module matrices
#   4) merged sample metadata table for ecological modelling
#
# NOTE: Large source matrices are not distributed in this repository. Edit paths below
# or pass paths through command-line arguments before running.

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(data.table)
  library(readr)
})

option_list <- list(
  make_option('--mag_abundance', type='character', default='data/MAG_TPM.tsv',
              help='MAG abundance table with columns sample_id, MAG_ID and TPM, or wide sample x MAG table.'),
  make_option('--mag_arg', type='character', default='data/MAG_ARG_annotations.tsv',
              help='Per-MAG ARG annotation table after RGI/CARD filtering.'),
  make_option('--mag_vf', type='character', default='data/MAG_VFDB_annotations.tsv',
              help='Per-MAG VFDB annotation table after ABRicate/VFDB filtering.'),
  make_option('--env', type='character', default='data/environmental_metadata.tsv',
              help='Sample-level environmental metadata table.'),
  make_option('--outdir', type='character', default='results/derived_matrices',
              help='Output directory.')
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

read_table_auto <- function(path) {
  stopifnot(file.exists(path))
  ext <- tools::file_ext(path)
  if (ext %in% c('csv')) readr::read_csv(path, show_col_types=FALSE) else readr::read_tsv(path, show_col_types=FALSE)
}

# ----- MAG abundance / sample x MAG presence -----
abund <- read_table_auto(opt$mag_abundance)

if (all(c('sample_id','MAG_ID','TPM') %in% names(abund))) {
  abund_long <- abund %>% transmute(sample_id=as.character(sample_id), MAG_ID=as.character(MAG_ID), TPM=as.numeric(TPM))
} else {
  # Wide table: first column sample_id, remaining columns MAG IDs.
  abund_long <- abund %>%
    rename(sample_id = 1) %>%
    pivot_longer(-sample_id, names_to='MAG_ID', values_to='TPM') %>%
    mutate(sample_id=as.character(sample_id), MAG_ID=as.character(MAG_ID), TPM=as.numeric(TPM))
}

presence <- abund_long %>%
  mutate(presence = as.integer(!is.na(TPM) & TPM > 0)) %>%
  select(sample_id, MAG_ID, presence) %>%
  pivot_wider(names_from=MAG_ID, values_from=presence, values_fill=0) %>%
  arrange(sample_id)
write_tsv(presence, file.path(opt$outdir, 'sample_MAG_presence.tsv'))

# ----- ARG annotations -----
arg <- read_table_auto(opt$mag_arg)
# Expected minimal columns: MAG_ID, arg_class, copy_number or count.
if (!'MAG_ID' %in% names(arg)) stop('ARG table must contain MAG_ID')
if (!'arg_class' %in% names(arg)) {
  candidate <- intersect(names(arg), c('drug_class','class','ARG_class'))[1]
  if (is.na(candidate)) stop('ARG table must contain arg_class/drug_class/class')
  arg <- rename(arg, arg_class = all_of(candidate))
}
if (!'copy_number' %in% names(arg)) {
  arg$copy_number <- 1
}

arg_class_mat <- arg %>%
  mutate(MAG_ID=as.character(MAG_ID), arg_class=as.character(arg_class), copy_number=as.numeric(copy_number)) %>%
  group_by(MAG_ID, arg_class) %>%
  summarise(ARG_copies=sum(copy_number, na.rm=TRUE), .groups='drop') %>%
  pivot_wider(names_from=arg_class, values_from=ARG_copies, values_fill=0)
write_tsv(arg_class_mat, file.path(opt$outdir, 'MAG_ARG_class_matrix.tsv'))

arg_mag <- arg_class_mat %>%
  mutate(ARG_load_MAG = rowSums(across(-MAG_ID), na.rm=TRUE)) %>%
  select(MAG_ID, ARG_load_MAG)

# ----- VFDB annotations -----
vf <- read_table_auto(opt$mag_vf)
if (!'MAG_ID' %in% names(vf)) stop('VFDB table must contain MAG_ID')
if (!'vf_module' %in% names(vf)) {
  candidate <- intersect(names(vf), c('module','VF_module','virulence_module'))[1]
  if (is.na(candidate)) stop('VFDB table must contain vf_module/module')
  vf <- rename(vf, vf_module = all_of(candidate))
}
if (!'count' %in% names(vf)) vf$count <- 1

vf_module_mat <- vf %>%
  mutate(MAG_ID=as.character(MAG_ID), vf_module=as.character(vf_module), count=as.numeric(count)) %>%
  group_by(MAG_ID, vf_module) %>%
  summarise(VF_count=sum(count, na.rm=TRUE), .groups='drop') %>%
  pivot_wider(names_from=vf_module, values_from=VF_count, values_fill=0)
write_tsv(vf_module_mat, file.path(opt$outdir, 'MAG_VF_module_matrix.tsv'))

vf_mag <- vf_module_mat %>%
  mutate(VIR_load_MAG = rowSums(across(-MAG_ID), na.rm=TRUE)) %>%
  select(MAG_ID, VIR_load_MAG)

# ----- sample-level summaries -----
mag_loads <- full_join(arg_mag, vf_mag, by='MAG_ID') %>%
  replace_na(list(ARG_load_MAG=0, VIR_load_MAG=0))

sample_burden <- abund_long %>%
  left_join(mag_loads, by='MAG_ID') %>%
  mutate(ARG_load_MAG=replace_na(ARG_load_MAG, 0), VIR_load_MAG=replace_na(VIR_load_MAG, 0)) %>%
  group_by(sample_id) %>%
  summarise(
    n_MAG_detected = sum(TPM > 0, na.rm=TRUE),
    ARG_load = sum(ARG_load_MAG * as.numeric(TPM > 0), na.rm=TRUE),
    VIR_load = sum(VIR_load_MAG * as.numeric(TPM > 0), na.rm=TRUE),
    ARG_abundance_weighted = sum(ARG_load_MAG * TPM, na.rm=TRUE),
    VIR_abundance_weighted = sum(VIR_load_MAG * TPM, na.rm=TRUE),
    .groups='drop'
  )
write_tsv(sample_burden, file.path(opt$outdir, 'sample_level_burdens.tsv'))

# ----- merge with environmental metadata -----
env <- read_table_auto(opt$env) %>% rename(sample_id = 1)
master <- sample_burden %>% left_join(env, by='sample_id')
write_tsv(master, file.path(opt$outdir, 'sample_metadata_with_burdens.tsv'))

message('Wrote derived matrices to: ', opt$outdir)
