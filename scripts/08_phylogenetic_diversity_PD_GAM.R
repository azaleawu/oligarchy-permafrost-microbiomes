#!/usr/bin/env Rscript
# Faith's phylogenetic diversity and diversity-burden GAM analyses.
# Generates PD/MPD/MNTD summaries, PD-burden GAMs, derivative-based PD50 estimates,
# segmented-regression sensitivity analyses and diagnostic/source tables.

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(ape)
  library(picante)
  library(mgcv)
  library(gratia)
  library(segmented)
  library(boot)
  library(readr)
})

option_list <- list(
  make_option('--tree', type='character', default='data/unified_MAG_tree.nwk'),
  make_option('--presence', type='character', default='results/derived_matrices/sample_MAG_presence.tsv'),
  make_option('--sample_data', type='character', default='results/derived_matrices/sample_metadata_with_burdens.tsv'),
  make_option('--outdir', type='character', default='results/phylogenetic_diversity'),
  make_option('--bootstrap', type='integer', default=999),
  make_option('--seed', type='integer', default=2025)
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)
set.seed(opt$seed)

# ----- helpers -----
fix_tree <- function(tr) {
  if (is.null(tr$edge.length)) tr$edge.length <- rep(1, nrow(tr$edge))
  pos <- tr$edge.length[tr$edge.length > 0]
  eps <- if (length(pos)) min(pos) * 1e-6 else 1e-6
  tr$edge.length[is.na(tr$edge.length) | tr$edge.length <= 0] <- eps
  ape::multi2di(tr)
}

fit_pd_gam <- function(df, response, k=8) {
  mgcv::gam(as.formula(paste0('log1p(', response, ') ~ s(PD, k=', k, ', bs="cs")')),
            data=df, method='REML', select=TRUE)
}

pd50_from_gam <- function(model) {
  dv <- gratia::derivatives(model, term='s(PD)')
  d0 <- max(dv$.derivative, na.rm=TRUE)
  out <- dv %>% filter(.derivative <= 0.5*d0) %>% slice_head(n=1)
  if (nrow(out) == 0) return(NA_real_)
  if ('PD' %in% names(out)) return(out$PD[1])
  out$data[1]
}

bootstrap_pd50 <- function(df, response, R=999, seed=2025) {
  set.seed(seed)
  boot_fun <- function(d, i) {
    m <- fit_pd_gam(d[i, , drop=FALSE], response=response, k=8)
    pd50_from_gam(m)
  }
  b <- boot::boot(df, boot_fun, R=R)
  vals <- na.omit(as.numeric(b$t[,1]))
  tibble(
    median = median(vals),
    ci_low = quantile(vals, 0.025),
    ci_high = quantile(vals, 0.975),
    n_success = length(vals),
    n_bootstrap = R
  )
}

segmented_pd <- function(df, response, psi_start) {
  dat <- df %>% mutate(y = log1p(.data[[response]])) %>% filter(is.finite(y), is.finite(PD))
  lm0 <- lm(y ~ PD, data=dat)
  sg <- segmented::segmented(lm0, seg.Z=~PD, psi=psi_start)
  dav <- segmented::davies.test(lm0, seg.Z=~PD)
  list(model=sg, davies=dav)
}

# ----- input -----
tr <- read.tree(opt$tree) |> fix_tree()
comm <- read_tsv(opt$presence, show_col_types=FALSE) %>%
  rename(sample_id = 1) %>%
  column_to_rownames('sample_id') %>%
  as.matrix()
mode(comm) <- 'numeric'

common_tips <- intersect(colnames(comm), tr$tip.label)
comm <- comm[, common_tips, drop=FALSE]
tr <- keep.tip(tr, common_tips)
comm <- comm[, tr$tip.label, drop=FALSE]

# ----- PD, MPD and MNTD -----
pd_df <- picante::pd(comm, tr, include.root=FALSE) %>%
  rownames_to_column('sample_id') %>%
  transmute(sample_id, PD=PD, SR=SR)

dist_mat <- cophenetic(tr)
mpd_df <- picante::ses.mpd(comm, dist_mat, null.model='taxa.labels', runs=999) %>%
  rownames_to_column('sample_id') %>%
  transmute(sample_id, MPD_Z=mpd.obs.z)
mntd_df <- picante::ses.mntd(comm, dist_mat, null.model='taxa.labels', runs=999) %>%
  rownames_to_column('sample_id') %>%
  transmute(sample_id, MNTD_Z=mntd.obs.z)

div <- reduce(list(pd_df, mpd_df, mntd_df), left_join, by='sample_id')
write_tsv(div, file.path(opt$outdir, 'sample_phylogenetic_diversity.tsv'))

sample_data <- read_tsv(opt$sample_data, show_col_types=FALSE) %>% rename(sample_id = 1)
master <- sample_data %>% left_join(div, by='sample_id')
write_tsv(master, file.path(opt$outdir, 'sample_data_with_PD_MPD_MNTD.tsv'))

# ----- PD-burden GAMs -----
endpoints <- c(ARG='ARG_load', Virulence='VIR_load')
summary_list <- list(); pd50_list <- list(); seg_list <- list()
for (nm in names(endpoints)) {
  response <- endpoints[[nm]]
  dat <- master %>% filter(!is.na(PD), !is.na(.data[[response]]))
  model <- fit_pd_gam(dat, response=response, k=8)
  saveRDS(model, file.path(opt$outdir, paste0('GAM_PD_', nm, '.rds')))
  sm <- summary(model)
  summary_list[[nm]] <- tibble(
    endpoint = nm,
    n = nrow(dat),
    adj_R2 = sm$r.sq,
    deviance_explained = sm$dev.expl,
    edf = sm$s.table[1, 'edf'],
    F = sm$s.table[1, 'F'],
    P_value = sm$s.table[1, 'p-value']
  )
  p50 <- bootstrap_pd50(dat, response=response, R=opt$bootstrap, seed=opt$seed)
  pd50_list[[nm]] <- p50 %>% mutate(endpoint=nm, .before=1)
  sg <- segmented_pd(dat, response=response, psi_start=p50$median[1])
  seg_list[[nm]] <- tibble(
    endpoint=nm,
    breakpoint=as.numeric(sg$model$psi[,'Est.']),
    Davies_P=sg$davies$p.value
  )
}
write_tsv(bind_rows(summary_list), file.path(opt$outdir, 'PD_GAM_model_summary.tsv'))
write_tsv(bind_rows(pd50_list), file.path(opt$outdir, 'PD50_bootstrap_summary.tsv'))
write_tsv(bind_rows(seg_list), file.path(opt$outdir, 'segmented_PD_summary.tsv'))

message('Phylogenetic diversity analyses completed: ', opt$outdir)
