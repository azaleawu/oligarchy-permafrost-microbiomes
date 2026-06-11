#!/usr/bin/env Rscript

# Environmental driver screening and multivariable GAMM modelling.
# This script implements the manuscript-level logic: repeated RF/XGBoost screening,
# VIF filtering, multivariable mixed-effect GAMM fitting and effect-size summaries.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(caret)
  library(ranger)
  library(xgboost)
  library(fastshap)
  library(car)
  library(mgcv)
  library(gamm4)
  library(spdep)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
config_file <- ifelse(length(args) >= 1, args[1], "config/analysis_paths_template.yml")
cfg <- yaml::read_yaml(config_file)
seed <- ifelse(is.null(cfg$analysis$random_seed), 2025, cfg$analysis$random_seed)
set.seed(seed)
out_dir <- cfg$output$tables
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

df <- fread(file.path(out_dir, "sample_level_burden_summary.tsv")) |> as.data.frame()

df <- df |> mutate(ARG_log = log1p(ARG_burden), VF_log = log1p(VF_burden))

candidate_vars <- intersect(
  c("GST", "ALT", "Elevation", "MAT", "MAP", "NDSI", "NDVI", "NPP", "srad", "wind",
    "ai", "aridity", "vapr", "bdod", "cfvo", "nitrogen", "ocd", "pH", "sand", "silt", "soc", "latitude"),
  names(df)
)

screen_variables <- function(data, response, candidate_vars, n_top = 10, threshold = 0.80) {
  data <- data[complete.cases(data[, c(response, candidate_vars)]), c(response, candidate_vars), drop = FALSE]
  folds <- caret::createMultiFolds(data[[response]], k = 5, times = 10)
  rf_freq <- setNames(integer(length(candidate_vars)), candidate_vars)
  xgb_freq <- setNames(integer(length(candidate_vars)), candidate_vars)
  X <- data[, candidate_vars, drop = FALSE]
  y <- data[[response]]

  for (i in seq_along(folds)) {
    idx <- folds[[i]]
    train <- data[idx, , drop = FALSE]

    rf <- ranger(
      formula = as.formula(paste(response, "~", paste(candidate_vars, collapse = "+"))),
      data = train,
      importance = "permutation",
      num.trees = 1000,
      seed = seed + i
    )
    top_rf <- names(sort(ranger::importance(rf), decreasing = TRUE))[seq_len(min(n_top, length(candidate_vars)))]
    rf_freq[top_rf] <- rf_freq[top_rf] + 1

    dtrain <- xgb.DMatrix(as.matrix(X[idx, , drop = FALSE]), label = y[idx])
    xgb <- xgb.train(
      params = list(objective = "reg:squarederror", eta = 0.05, max_depth = 6,
                    subsample = 0.8, colsample_bytree = 0.8, min_child_weight = 3, gamma = 0),
      data = dtrain,
      nrounds = 600,
      verbose = 0
    )
    shap <- fastshap::explain(
      object = xgb,
      X = X[idx, , drop = FALSE],
      pred_wrapper = function(object, newdata) predict(object, xgb.DMatrix(as.matrix(newdata))),
      nsim = 50,
      adjust = TRUE
    )
    top_xgb <- names(sort(colMeans(abs(shap)), decreasing = TRUE))[seq_len(min(n_top, length(candidate_vars)))]
    xgb_freq[top_xgb] <- xgb_freq[top_xgb] + 1
  }

  n_folds <- length(folds)
  rates <- data.frame(
    variable = candidate_vars,
    rf_rate = as.numeric(rf_freq) / n_folds,
    xgb_rate = as.numeric(xgb_freq) / n_folds
  ) |>
    mutate(stable = rf_rate >= threshold | xgb_rate >= threshold)
  rates
}

vif_filter <- function(data, response, vars, threshold = 7) {
  vars <- vars[vars %in% names(data)]
  dat <- data[complete.cases(data[, c(response, vars)]), c(response, vars), drop = FALSE]
  repeat {
    if (length(vars) <= 1) break
    fit <- lm(as.formula(paste(response, "~", paste(vars, collapse = "+"))), data = dat)
    vf <- car::vif(fit)
    if (max(vf) <= threshold) break
    drop_var <- names(which.max(vf))
    vars <- setdiff(vars, drop_var)
    dat <- data[complete.cases(data[, c(response, vars)]), c(response, vars), drop = FALSE]
  }
  vars
}

fit_multivariable_gamm <- function(data, response, vars) {
  smooth_terms <- paste0("s(", vars, ", k = 6, bs = 'cs')", collapse = " + ")
  formula <- as.formula(paste(response, "~", smooth_terms))
  gamm4(formula, random = ~(1 | region), data = data)
}

effect_delta <- function(gam_obj, data, var, q1 = 0.10, q2 = 0.90) {
  model_vars <- all.vars(formula(gam_obj))[-1]
  new1 <- as.data.frame(lapply(data[, model_vars, drop = FALSE], median, na.rm = TRUE))
  new2 <- new1
  new1[[var]] <- quantile(data[[var]], q1, na.rm = TRUE)
  new2[[var]] <- quantile(data[[var]], q2, na.rm = TRUE)
  X1 <- predict(gam_obj, newdata = new1, type = "lpmatrix")
  X2 <- predict(gam_obj, newdata = new2, type = "lpmatrix")
  beta <- coef(gam_obj)
  Vb <- vcov(gam_obj)
  diff_eta <- drop((X2 - X1) %*% beta)
  se_diff <- sqrt(drop((X2 - X1) %*% Vb %*% t(X2 - X1)))
  data.frame(
    variable = var,
    delta_log = diff_eta,
    fold_change = exp(diff_eta),
    lower_fold = exp(diff_eta - 1.96 * se_diff),
    upper_fold = exp(diff_eta + 1.96 * se_diff)
  )
}

run_endpoint <- function(response, endpoint) {
  screening <- screen_variables(df, response, candidate_vars,
                                n_top = cfg$analysis$top_n_features,
                                threshold = cfg$analysis$stable_frequency_threshold)
  fwrite(screening, file.path(out_dir, paste0(endpoint, "_environmental_screening.tsv")), sep = "\t")
  stable_vars <- screening$variable[screening$stable]
  final_vars <- vif_filter(df, response, stable_vars, threshold = 7)
  model_data <- df[complete.cases(df[, c(response, final_vars, "region")]), c(response, final_vars, "region", "longitude", "latitude"), drop = FALSE]
  fit <- fit_multivariable_gamm(model_data, response, final_vars)
  eff <- bind_rows(lapply(final_vars, function(v) effect_delta(fit$gam, model_data, v)))
  eff$endpoint <- endpoint
  fwrite(eff, file.path(out_dir, paste0(endpoint, "_GAMM_effect_sizes.tsv")), sep = "\t")
  invisible(fit)
}

if (length(candidate_vars) >= 2 && "region" %in% names(df)) {
  run_endpoint("ARG_log", "ARG")
  run_endpoint("VF_log", "Virulence")
} else {
  warning("Insufficient environmental variables or missing region column; skipping environmental GAMM analysis.")
}
