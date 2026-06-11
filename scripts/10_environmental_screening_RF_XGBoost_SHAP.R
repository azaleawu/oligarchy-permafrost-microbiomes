#!/usr/bin/env Rscript
# Two-stage environmental driver analysis: RF/XGBoost-SHAP screening followed by
# multivariable mixed-effect GAMMs. This script is used for both ARG and virulence endpoints.

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(data.table)
  library(caret)
  library(ranger)
  library(xgboost)
  library(fastshap)
  library(car)
  library(mgcv)
  library(gamm4)
  library(spdep)
  library(readr)
})

option_list <- list(
  make_option('--input', type='character', default='results/phylogenetic_diversity/sample_data_with_PD_MPD_MNTD.tsv'),
  make_option('--endpoint', type='character', default='VIR_load', help='Response column, e.g. ARG_load or VIR_load.'),
  make_option('--region', type='character', default='region'),
  make_option('--longitude', type='character', default='longitude'),
  make_option('--latitude', type='character', default='latitude'),
  make_option('--candidate_vars', type='character', default='GST,ALT,Elevation,MAT,MAP,NDSI,NDVI,NPP,srad,wind,ai,aridity,vapr,bdod,cfvo,nitrogen,ocd,pH,sand,silt,soc,latitude'),
  make_option('--outdir', type='character', default='results/environmental_models'),
  make_option('--seed', type='integer', default=2025),
  make_option('--nfold', type='integer', default=5),
  make_option('--repeats', type='integer', default=10),
  make_option('--topn', type='integer', default=10),
  make_option('--stable_rate', type='double', default=0.8),
  make_option('--vif_threshold', type='double', default=7)
)
opt <- parse_args(OptionParser(option_list=option_list))
set.seed(opt$seed)
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

dt <- as.data.table(read_tsv(opt$input, show_col_types=FALSE))
if (!opt$endpoint %in% names(dt)) stop('Endpoint not found: ', opt$endpoint)
dt[, Y := log1p(get(opt$endpoint))]

candidate_vars <- trimws(strsplit(opt$candidate_vars, ',')[[1]])
candidate_vars <- intersect(candidate_vars, names(dt))
# Keep numeric predictors only.
candidate_vars <- candidate_vars[sapply(dt[, ..candidate_vars], is.numeric)]

# Remove rows with missing response, predictors, region and coordinates.
needed <- unique(c('Y', candidate_vars, opt$region, opt$longitude, opt$latitude))
needed <- intersect(needed, names(dt))
dt <- dt[complete.cases(dt[, ..needed])]

# Transform strongly skewed predictors with log1p where values are non-negative, then z-standardize.
skewness <- function(x) mean((x - mean(x, na.rm=TRUE))^3, na.rm=TRUE) / sd(x, na.rm=TRUE)^3
skew_tbl <- tibble(var=candidate_vars, skewness=sapply(candidate_vars, function(v) skewness(dt[[v]])))
vars_to_log <- skew_tbl %>% filter(abs(skewness) > 1) %>% pull(var)
vars_to_log <- vars_to_log[sapply(dt[, ..vars_to_log], function(x) all(x >= 0, na.rm=TRUE))]
for (v in vars_to_log) dt[[v]] <- log1p(dt[[v]])
for (v in candidate_vars) dt[[v]] <- as.numeric(scale(dt[[v]]))
write_tsv(skew_tbl, file.path(opt$outdir, paste0(opt$endpoint, '_predictor_skewness.tsv')))

X_all <- as.data.frame(dt[, ..candidate_vars])
y_all <- dt$Y
folds <- caret::createMultiFolds(y_all, k=opt$nfold, times=opt$repeats)
rf_freq <- setNames(integer(length(candidate_vars)), candidate_vars)
xgb_freq <- setNames(integer(length(candidate_vars)), candidate_vars)

xgb_params <- list(
  booster='gbtree', objective='reg:squarederror', eta=0.05, max_depth=6,
  subsample=0.8, colsample_bytree=0.8, min_child_weight=3, gamma=0
)

for (i in seq_along(folds)) {
  idx <- folds[[i]]
  df_rf <- X_all[idx, , drop=FALSE]
  df_rf$Y <- y_all[idx]
  rf_mod <- ranger::ranger(Y ~ ., data=df_rf, importance='permutation', num.trees=1200, seed=opt$seed+i)
  top_rf <- names(sort(ranger::importance(rf_mod), decreasing=TRUE))[seq_len(min(opt$topn, length(candidate_vars)))]
  rf_freq[top_rf] <- rf_freq[top_rf] + 1

  dtr <- xgb.DMatrix(as.matrix(X_all[idx, , drop=FALSE]), label=y_all[idx])
  xgb_mod <- xgb.train(params=xgb_params, data=dtr, nrounds=600, verbose=0)
  X_df <- X_all[idx, , drop=FALSE]
  colnames(X_df) <- xgb_mod$feature_names
  shap_vals <- fastshap::explain(
    object=xgb_mod, X=X_df,
    pred_wrapper=function(m, newdata) predict(m, xgb.DMatrix(as.matrix(newdata))),
    nsim=50, adjust=TRUE
  )
  top_xgb <- names(sort(colMeans(abs(shap_vals)), decreasing=TRUE))[seq_len(min(opt$topn, length(candidate_vars)))]
  xgb_freq[top_xgb] <- xgb_freq[top_xgb] + 1
}

n_resamples <- length(folds)
rate_tbl <- tibble(var=candidate_vars, RF_rate=rf_freq/n_resamples, XGB_SHAP_rate=xgb_freq/n_resamples) %>%
  arrange(desc(pmax(RF_rate, XGB_SHAP_rate)))
write_tsv(rate_tbl, file.path(opt$outdir, paste0(opt$endpoint, '_RF_XGB_SHAP_stability.tsv')))
stable_vars <- rate_tbl %>% filter(RF_rate >= opt$stable_rate | XGB_SHAP_rate >= opt$stable_rate) %>% pull(var)
if (length(stable_vars) == 0) stop('No stable predictors selected. Consider lowering stable_rate after documenting the change.')

# Iterative VIF filtering.
drop_high_vif <- function(df, y, vars, threshold) {
  dat <- as.data.frame(df[, c(y, vars), with=FALSE])
  repeat {
    if (length(vars) <= 1) break
    vf <- car::vif(lm(as.formula(paste0(y, ' ~ .')), data=dat))
    if (max(vf, na.rm=TRUE) <= threshold) break
    drop_var <- names(which.max(vf))
    vars <- setdiff(vars, drop_var)
    dat <- dat[, c(y, vars), drop=FALSE]
  }
  final_vif <- if (length(vars) > 1) car::vif(lm(as.formula(paste0(y, ' ~ .')), data=dat)) else setNames(1, vars)
  list(vars=vars, vif=tibble(var=names(final_vif), VIF=as.numeric(final_vif)))
}
vif_res <- drop_high_vif(dt, 'Y', stable_vars, opt$vif_threshold)
final_vars <- vif_res$vars
write_tsv(vif_res$vif, file.path(opt$outdir, paste0(opt$endpoint, '_VIF_final_predictors.tsv')))
write_lines(final_vars, file.path(opt$outdir, paste0(opt$endpoint, '_final_predictors.txt')))

# Multivariable mixed-effect GAMM: all retained predictors enter simultaneously as smooth terms.
fml <- as.formula(paste0('Y ~ ', paste0('s(', final_vars, ', k=6, bs="cs")', collapse=' + ')))
rand_fml <- as.formula(paste0('~(1|', opt$region, ')'))
mgam <- gamm4::gamm4(fml, random=rand_fml, data=as.data.frame(dt))
saveRDS(mgam, file.path(opt$outdir, paste0(opt$endpoint, '_multivariable_GAMM.rds')))

sm <- summary(mgam$gam)
term_tbl <- as.data.frame(sm$s.table) %>% rownames_to_column('term') %>%
  transmute(term, edf=edf, F=F, P_value=`p-value`)
model_tbl <- tibble(endpoint=opt$endpoint, n=nrow(dt), adj_R2=sm$r.sq, deviance_explained=sm$dev.expl)
write_tsv(term_tbl, file.path(opt$outdir, paste0(opt$endpoint, '_GAMM_terms.tsv')))
write_tsv(model_tbl, file.path(opt$outdir, paste0(opt$endpoint, '_GAMM_model_fit.tsv')))

# 10th -> 90th percentile contrasts, holding other predictors at medians.
effect_delta <- function(gam_obj, data, var, q1=0.10, q2=0.90) {
  vars <- all.vars(formula(gam_obj))[-1]
  new1 <- data.frame(lapply(data[, vars, drop=FALSE], median, na.rm=TRUE))
  new2 <- new1
  new1[[var]] <- as.numeric(quantile(data[[var]], q1, na.rm=TRUE))
  new2[[var]] <- as.numeric(quantile(data[[var]], q2, na.rm=TRUE))
  X1 <- predict(gam_obj, newdata=new1, type='lpmatrix')
  X2 <- predict(gam_obj, newdata=new2, type='lpmatrix')
  beta <- coef(gam_obj); Vb <- vcov(gam_obj)
  diff_eta <- drop((X2-X1) %*% beta)
  se <- sqrt(drop((X2-X1) %*% Vb %*% t(X2-X1)))
  tibble(var=var,
         q10=new1[[var]], q90=new2[[var]], delta_log=diff_eta,
         delta_log_low=diff_eta-1.96*se, delta_log_high=diff_eta+1.96*se,
         fold_change=exp(diff_eta), fold_low=exp(diff_eta-1.96*se), fold_high=exp(diff_eta+1.96*se))
}
eff_tbl <- map_dfr(final_vars, ~effect_delta(mgam$gam, as.data.frame(dt), .x))
write_tsv(eff_tbl, file.path(opt$outdir, paste0(opt$endpoint, '_GAMM_10_to_90_effects.tsv')))

# Residual spatial autocorrelation.
if (all(c(opt$longitude, opt$latitude) %in% names(dt))) {
  coords <- as.matrix(dt[, c(opt$longitude, opt$latitude), with=FALSE])
  nb <- spdep::knn2nb(spdep::knearneigh(coords, k=5))
  lw <- spdep::nb2listw(nb, style='W')
  mi <- spdep::moran.test(residuals(mgam$gam), lw)
  capture.output(mi, file=file.path(opt$outdir, paste0(opt$endpoint, '_MoranI_residuals.txt')))
}

message('Environmental screening and GAMM completed for endpoint: ', opt$endpoint)
