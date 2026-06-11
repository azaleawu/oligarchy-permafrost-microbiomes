#!/usr/bin/env Rscript
# Generate analysis panels and source-data tables used for final figures.
# Final composite figures may be assembled and annotated in Adobe Illustrator; this script
# documents the analytical panels generated from the data.

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(ggplot2)
  library(readr)
})

option_list <- list(
  make_option('--sample_data', type='character', default='results/phylogenetic_diversity/sample_data_with_PD_MPD_MNTD.tsv'),
  make_option('--env_effects_ARG', type='character', default='results/environmental_models/ARG_load_GAMM_10_to_90_effects.tsv'),
  make_option('--env_effects_VIR', type='character', default='results/environmental_models/VIR_load_GAMM_10_to_90_effects.tsv'),
  make_option('--outdir', type='character', default='figures/analysis_panels')
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

theme_publication <- function(base_size=9) {
  theme_bw(base_size=base_size) +
    theme(panel.grid=element_blank(), strip.background=element_rect(fill='grey95'), legend.position='right')
}

if (file.exists(opt$sample_data)) {
  dat <- read_tsv(opt$sample_data, show_col_types=FALSE)
  p1 <- ggplot(dat, aes(PD, log1p(ARG_load))) +
    geom_point(alpha=0.55, size=1.5) +
    geom_smooth(method='gam', formula=y ~ s(x, k=8, bs='cs'), se=TRUE) +
    labs(x="Faith's phylogenetic diversity", y='log(1 + ARG load)') + theme_publication()
  ggsave(file.path(opt$outdir, 'panel_PD_ARG.pdf'), p1, width=4, height=3, device=cairo_pdf)

  p2 <- ggplot(dat, aes(PD, log1p(VIR_load))) +
    geom_point(alpha=0.55, size=1.5) +
    geom_smooth(method='gam', formula=y ~ s(x, k=8, bs='cs'), se=TRUE) +
    labs(x="Faith's phylogenetic diversity", y='log(1 + virulence-associated gene load)') + theme_publication()
  ggsave(file.path(opt$outdir, 'panel_PD_VIR.pdf'), p2, width=4, height=3, device=cairo_pdf)
}

plot_effects <- function(path, endpoint) {
  if (!file.exists(path)) return(invisible(NULL))
  eff <- read_tsv(path, show_col_types=FALSE) %>% arrange(delta_log)
  p <- ggplot(eff, aes(delta_log, reorder(var, delta_log))) +
    geom_vline(xintercept=0, linetype='dashed', linewidth=0.3) +
    geom_errorbarh(aes(xmin=delta_log_low, xmax=delta_log_high), height=0.2) +
    geom_point(size=1.8) +
    labs(x='10th to 90th percentile contrast on log scale', y=NULL, title=endpoint) + theme_publication()
  ggsave(file.path(opt$outdir, paste0('panel_environmental_effects_', endpoint, '.pdf')), p, width=4.8, height=3.5, device=cairo_pdf)
}
plot_effects(opt$env_effects_ARG, 'ARG')
plot_effects(opt$env_effects_VIR, 'Virulence')

message('Analysis panels written to: ', opt$outdir)
