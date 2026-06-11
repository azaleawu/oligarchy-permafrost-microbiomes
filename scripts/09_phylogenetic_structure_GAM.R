#!/usr/bin/env Rscript
# Deep- and tip-level phylogenetic structure analyses using MPD_Z and MNTD_Z.
# Compares single-predictor, additive and tensor-product GAMs for ARG and virulence endpoints.

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(mgcv)
  library(gratia)
  library(readr)
})

option_list <- list(
  make_option('--sample_data', type='character', default='results/phylogenetic_diversity/sample_data_with_PD_MPD_MNTD.tsv'),
  make_option('--outdir', type='character', default='results/phylogenetic_structure'),
  make_option('--seed', type='integer', default=2025)
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)
set.seed(opt$seed)

dat <- read_tsv(opt$sample_data, show_col_types=FALSE)

fit_candidates <- function(dat, response) {
  dat <- dat %>% filter(!is.na(MPD_Z), !is.na(MNTD_Z), !is.na(.data[[response]])) %>%
    mutate(Y = log1p(.data[[response]]))
  list(
    MPD = gam(Y ~ s(MPD_Z, k=6, bs='cs'), data=dat, method='REML', select=TRUE),
    MNTD = gam(Y ~ s(MNTD_Z, k=6, bs='cs'), data=dat, method='REML', select=TRUE),
    additive = gam(Y ~ s(MPD_Z, k=6, bs='cs') + s(MNTD_Z, k=6, bs='cs'), data=dat, method='REML', select=TRUE),
    tensor = gam(Y ~ te(MPD_Z, MNTD_Z, k=c(6,6), bs=c('tp','tp')), data=dat, method='REML')
  )
}

summarise_model <- function(model, endpoint, model_name) {
  sm <- summary(model)
  tibble(
    endpoint=endpoint,
    model=model_name,
    n=nobs(model),
    AIC=AIC(model),
    adj_R2=sm$r.sq,
    deviance_explained=sm$dev.expl,
    total_edf=sum(sm$s.table[,'edf']),
    min_P=min(sm$s.table[,'p-value'], na.rm=TRUE)
  )
}

endpoints <- c(ARG='ARG_load', Virulence='VIR_load')
model_summaries <- list(); term_summaries <- list()
for (nm in names(endpoints)) {
  response <- endpoints[[nm]]
  mods <- fit_candidates(dat, response)
  saveRDS(mods, file.path(opt$outdir, paste0('candidate_GAMs_', nm, '.rds')))
  model_summaries[[nm]] <- imap_dfr(mods, ~summarise_model(.x, nm, .y))
  term_summaries[[nm]] <- imap_dfr(mods, function(m, model_name){
    st <- as.data.frame(summary(m)$s.table) %>% rownames_to_column('term')
    st %>% transmute(endpoint=nm, model=model_name, term, edf=edf, F=F, P_value=`p-value`)
  })
}
write_tsv(bind_rows(model_summaries), file.path(opt$outdir, 'phylogenetic_structure_model_selection.tsv'))
write_tsv(bind_rows(term_summaries), file.path(opt$outdir, 'phylogenetic_structure_term_significance.tsv'))

# Pairwise rank correlations for reporting.
cor_tbl <- dat %>%
  select(MPD_Z, MNTD_Z, ARG_load, VIR_load) %>%
  summarise(
    rho_MPD_MNTD = cor(MPD_Z, MNTD_Z, method='spearman', use='complete.obs'),
    rho_MPD_ARG = cor(MPD_Z, ARG_load, method='spearman', use='complete.obs'),
    rho_MNTD_ARG = cor(MNTD_Z, ARG_load, method='spearman', use='complete.obs'),
    rho_MPD_VIR = cor(MPD_Z, VIR_load, method='spearman', use='complete.obs'),
    rho_MNTD_VIR = cor(MNTD_Z, VIR_load, method='spearman', use='complete.obs')
  )
write_tsv(cor_tbl, file.path(opt$outdir, 'phylogenetic_structure_spearman_correlations.tsv'))

message('Phylogenetic structure GAM analyses completed: ', opt$outdir)
