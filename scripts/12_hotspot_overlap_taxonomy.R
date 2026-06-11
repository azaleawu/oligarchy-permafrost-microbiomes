#!/usr/bin/env Rscript
# ARG-VIR hotspot overlap and taxonomic enrichment analyses.

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(readr)
})

option_list <- list(
  make_option('--arg_membership', type='character', default='results/phylofactor_hotspots/ARG_hotspot_membership.tsv'),
  make_option('--vir_membership', type='character', default='results/phylofactor_hotspots/VIR_hotspot_membership.tsv'),
  make_option('--taxonomy', type='character', default='data/MAG_taxonomy.tsv'),
  make_option('--outdir', type='character', default='results/hotspot_overlap_taxonomy')
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

arg <- read_tsv(opt$arg_membership, show_col_types=FALSE) %>% rename(ARG_hotspot=hotspot)
vir <- read_tsv(opt$vir_membership, show_col_types=FALSE) %>% rename(VIR_hotspot=hotspot)
mem <- full_join(arg, vir, by='MAG_ID') %>%
  mutate(ARG_hotspot=replace_na(ARG_hotspot, 'Background'), VIR_hotspot=replace_na(VIR_hotspot, 'Background'))

# Pairwise Jaccard and Fisher exact tests for A x P hotspot overlaps.
A_levels <- setdiff(sort(unique(mem$ARG_hotspot)), 'Background')
P_levels <- setdiff(sort(unique(mem$VIR_hotspot)), 'Background')
overlap_tbl <- crossing(ARG_hotspot=A_levels, VIR_hotspot=P_levels) %>%
  mutate(stats = map2(ARG_hotspot, VIR_hotspot, function(a, p){
    A <- mem$ARG_hotspot == a
    P <- mem$VIR_hotspot == p
    tab <- table(factor(A, levels=c(FALSE, TRUE)), factor(P, levels=c(FALSE, TRUE)))
    shared <- sum(A & P)
    union_n <- sum(A | P)
    tibble(shared_MAGs=shared, Jaccard=ifelse(union_n > 0, shared/union_n, NA_real_), P_value=fisher.test(tab)$p.value)
  })) %>% unnest(stats) %>% mutate(P_adj=p.adjust(P_value, method='BH'))
write_tsv(overlap_tbl, file.path(opt$outdir, 'ARG_VIR_hotspot_overlap.tsv'))

# Taxonomic composition and enrichment.
tax <- read_tsv(opt$taxonomy, show_col_types=FALSE) %>% rename(MAG_ID=1)
if (!'phylum' %in% names(tax)) {
  cand <- intersect(names(tax), c('phylum_clean','Phylum','phylum_grp'))[1]
  if (is.na(cand)) stop('Taxonomy table must contain phylum/phylum_clean/Phylum/phylum_grp')
  tax <- rename(tax, phylum=all_of(cand))
}
mem_tax <- mem %>% left_join(tax %>% select(MAG_ID, phylum), by='MAG_ID') %>% mutate(phylum=replace_na(phylum, 'Unclassified'))

tax_test <- function(df, hotspot_col, hotspot_value) {
  focal <- df[[hotspot_col]] == hotspot_value
  tab <- table(factor(focal, levels=c(TRUE, FALSE)), df$phylum)
  use_mc <- sum(tab) < 200
  ct <- suppressWarnings(chisq.test(tab, simulate.p.value=use_mc, B=ifelse(use_mc, 100000, 0)))
  tibble(hotspot=hotspot_value, chisq=unname(ct$statistic), df=unname(ct$parameter), P_value=ct$p.value)
}

tax_arg <- map_dfr(A_levels, ~tax_test(mem_tax, 'ARG_hotspot', .x)) %>% mutate(endpoint='ARG')
tax_vir <- map_dfr(P_levels, ~tax_test(mem_tax, 'VIR_hotspot', .x)) %>% mutate(endpoint='VIR')
tax_tests <- bind_rows(tax_arg, tax_vir) %>% group_by(endpoint) %>% mutate(P_adj=p.adjust(P_value, method='BH')) %>% ungroup()
write_tsv(tax_tests, file.path(opt$outdir, 'hotspot_taxonomic_enrichment_tests.tsv'))

tax_comp <- bind_rows(
  mem_tax %>% filter(ARG_hotspot != 'Background') %>% count(endpoint='ARG', hotspot=ARG_hotspot, phylum, name='n'),
  mem_tax %>% filter(VIR_hotspot != 'Background') %>% count(endpoint='VIR', hotspot=VIR_hotspot, phylum, name='n')
) %>% group_by(endpoint, hotspot) %>% mutate(prop=n/sum(n)) %>% ungroup()
write_tsv(tax_comp, file.path(opt$outdir, 'hotspot_taxonomic_composition.tsv'))

message('Hotspot overlap and taxonomy analyses completed: ', opt$outdir)
