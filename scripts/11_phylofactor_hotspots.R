#!/usr/bin/env Rscript
# Phylogenetic factorisation and hotspot inference for ARG-class and virulence-module matrices.
# The script runs a fixed number of factors (default nfactors=5), maps focal tip sets to
# hotspot labels, computes enrichment summaries and PERMANOVA effect sizes.

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(ape)
  library(vegan)
  library(readr)
})

option_list <- list(
  make_option('--tree', type='character', default='data/unified_MAG_tree.nwk'),
  make_option('--arg_matrix', type='character', default='results/derived_matrices/MAG_ARG_class_matrix.tsv'),
  make_option('--vir_matrix', type='character', default='results/derived_matrices/MAG_VF_module_matrix.tsv'),
  make_option('--outdir', type='character', default='results/phylofactor_hotspots'),
  make_option('--nfactors', type='integer', default=5),
  make_option('--seed', type='integer', default=2025)
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)
set.seed(opt$seed)

read_matrix <- function(path) {
  read_tsv(path, show_col_types=FALSE) %>%
    rename(MAG_ID=1) %>%
    column_to_rownames('MAG_ID') %>%
    as.matrix()
}

prepare_matrix <- function(mat) {
  mat[is.na(mat) | mat < 0] <- 0
  min_pos <- suppressWarnings(min(mat[mat > 0], na.rm=TRUE))
  if (is.finite(min_pos)) mat[mat == 0] <- 0.65 * min_pos
  mat
}

# A light wrapper is used because phylofactor package APIs can differ across versions.
run_phylofactor <- function(mat, tree, nfactors) {
  if (!requireNamespace('phylofactor', quietly=TRUE)) {
    stop('The phylofactor package is required. Install the version used in the manuscript before running this script.')
  }
  phylofactor::PhyloFactor(
    x = mat,
    tree = tree,
    nfactors = nfactors,
    method = 'max.var',
    choice = 'var'
  )
}

extract_membership <- function(pf_obj, mat, prefix) {
  # Expected object structure follows the original phylofactor output used in the analysis:
  # pf_obj$groups[[i]][[1]] indexes focal tips for factor i.
  if (is.null(pf_obj$groups) || is.null(pf_obj$tree$tip.label)) {
    stop('Unexpected phylofactor object structure; inspect object and adapt extract_membership().')
  }
  out <- map_dfr(seq_along(pf_obj$groups), function(i) {
    tips <- pf_obj$tree$tip.label[pf_obj$groups[[i]][[1]]]
    tibble(MAG_ID=tips, hotspot=paste0(prefix, i))
  })
  all_mags <- tibble(MAG_ID=rownames(mat))
  all_mags %>% left_join(out, by='MAG_ID') %>% mutate(hotspot=replace_na(hotspot, 'Background'))
}

enrichment_by_hotspot <- function(mat, membership, prefix) {
  x <- as.data.frame(mat) %>% rownames_to_column('MAG_ID') %>% left_join(membership, by='MAG_ID')
  classes <- colnames(mat)
  map_dfr(setdiff(unique(membership$hotspot), 'Background'), function(h) {
    map_dfr(classes, function(cls) {
      focal <- x %>% filter(hotspot == h) %>% pull(all_of(cls))
      background <- x %>% filter(hotspot != h) %>% pull(all_of(cls))
      wt <- suppressWarnings(wilcox.test(focal, background, exact=FALSE))
      tibble(
        hotspot=h, feature=cls,
        mean_focal=mean(focal, na.rm=TRUE), mean_background=mean(background, na.rm=TRUE),
        log2FC=log2((mean(focal, na.rm=TRUE)+1e-9)/(mean(background, na.rm=TRUE)+1e-9)),
        P_value=wt$p.value
      )
    }) %>% mutate(P_adj=p.adjust(P_value, method='BH'))
  })
}

permanova_hotspot <- function(mat, membership) {
  X <- log10(mat + 1)
  meta <- membership %>% filter(MAG_ID %in% rownames(X)) %>% arrange(match(MAG_ID, rownames(X)))
  ad <- vegan::adonis2(X ~ hotspot, data=meta, method='euclidean', permutations=9999)
  as.data.frame(ad) %>% rownames_to_column('term')
}

tr <- read.tree(opt$tree)
arg_mat <- read_matrix(opt$arg_matrix)
vir_mat <- read_matrix(opt$vir_matrix)

common_arg <- intersect(tr$tip.label, rownames(arg_mat)); tree_arg <- keep.tip(tr, common_arg); arg_mat <- arg_mat[tree_arg$tip.label,,drop=FALSE]
common_vir <- intersect(tr$tip.label, rownames(vir_mat)); tree_vir <- keep.tip(tr, common_vir); vir_mat <- vir_mat[tree_vir$tip.label,,drop=FALSE]

arg_pf <- run_phylofactor(prepare_matrix(arg_mat), tree_arg, opt$nfactors)
vir_pf <- run_phylofactor(prepare_matrix(vir_mat), tree_vir, opt$nfactors)
saveRDS(arg_pf, file.path(opt$outdir, 'phylofactor_ARG.rds'))
saveRDS(vir_pf, file.path(opt$outdir, 'phylofactor_VIR.rds'))

arg_mem <- extract_membership(arg_pf, arg_mat, 'A')
vir_mem <- extract_membership(vir_pf, vir_mat, 'P')
write_tsv(arg_mem, file.path(opt$outdir, 'ARG_hotspot_membership.tsv'))
write_tsv(vir_mem, file.path(opt$outdir, 'VIR_hotspot_membership.tsv'))

write_tsv(enrichment_by_hotspot(arg_mat, arg_mem, 'A'), file.path(opt$outdir, 'ARG_hotspot_feature_enrichment.tsv'))
write_tsv(enrichment_by_hotspot(vir_mat, vir_mem, 'P'), file.path(opt$outdir, 'VIR_hotspot_feature_enrichment.tsv'))
write_tsv(permanova_hotspot(arg_mat, arg_mem), file.path(opt$outdir, 'PERMANOVA_ARG_hotspots.tsv'))
write_tsv(permanova_hotspot(vir_mat, vir_mem), file.path(opt$outdir, 'PERMANOVA_VIR_hotspots.tsv'))

message('Phylogenetic factorisation completed: ', opt$outdir)
