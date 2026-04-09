library(dplyr)
library(glmnet)
library(pROC)
library(ggplot2)
library(Matrix)
library(stringr)
library(sharp)
library(parallelly)
library(tidyr)

# =========================
# Global settings
# =========================
setwd("/rds/general/project/hda_25-26/live/TDS/TDS_Group1")

base_out_dir <- "modelling_script/stability_analysis/subsample_lasso_sex"
dir.create(base_out_dir, recursive = TRUE, showWarnings = FALSE)

get_requested_cores <- function() {
  candidates <- c(
    Sys.getenv("PBS_NCPUS", unset = ""),
    Sys.getenv("PBS_NP", unset = ""),
    Sys.getenv("NCPUS", unset = ""),
    Sys.getenv("MY_NCORES", unset = "")
  )
  
  candidates <- suppressWarnings(as.integer(candidates))
  candidates <- candidates[!is.na(candidates) & candidates > 0]
  
  if (length(candidates) == 0) return(1L)
  candidates[1]
}

requested_cores <- get_requested_cores()
available_cores <- parallelly::availableCores()
n_cores <- min(requested_cores, available_cores)

nfolds_cv <- 10
n_boot <- 100

# =========================
# Load data
# =========================
train_df <- readRDS("split_imputed_data/ukb_G1_train_imputed.rds")
test_df  <- readRDS("split_imputed_data/ukb_G1_test_imputed.rds")

outcome <- "cvd_incident"

# =========================
# Define predictors (sex has been stratified)
# =========================
confounders <- c("age", "ethnicity_5cat")

exposome_list <- c(
  "tv_group",
  "total_met_group",
  "smoking_3cat",
  "alcohol_combined",
  "redwine_group",
  "salt_3cat",
  "education_2",
  "employment_2cat",
  "imd_quintile",
  "sf_score",
  "stress_group_2yr"
)

pattern_regex <- paste(
  c(
    "fruit", "fish", "tea", "coffee",
    "sleep", "neurotic", "satis",
    "no2", "pm10", "pm2_5", "noise",
    "greenspace", "temp",
    "living", "urban", "work_hours"
  ),
  collapse = "|"
)

pattern_vars <- grep(pattern_regex, names(train_df), value = TRUE)

vars <- unique(c(exposome_list, pattern_vars))
vars <- vars[vars %in% names(train_df) & vars %in% names(test_df)]

preds <- unique(c(confounders, vars))

# =========================
# Helper functions
# =========================
make_pretty_names <- function(x) {
  x |>
    str_replace_all("_", " ") |>
    str_squish() |>
    str_replace_all("pm2 5", "PM2.5") |>
    str_replace_all("pm10", "PM10") |>
    str_replace_all("no2", "NO2") |>
    str_replace_all("imd", "IMD") |>
    str_replace_all("bmi", "BMI") |>
    str_replace_all("tv", "TV") |>
    str_replace_all("cvd", "CVD")
}

safe_auc <- function(y_true, p_hat) {
  tryCatch(
    as.numeric(auc(roc(y_true, p_hat, quiet = TRUE))),
    error = function(e) NA_real_
  )
}

cv_auc_from_selected <- function(X, y, selected_terms, fold_id) {
  if (length(selected_terms) == 0) return(NA_real_)
  
  pred_all <- rep(NA_real_, length(y))
  
  for (k in sort(unique(fold_id))) {
    train_idx <- which(fold_id != k)
    valid_idx <- which(fold_id == k)
    
    X_train <- as.data.frame(X[train_idx, selected_terms, drop = FALSE])
    X_valid <- as.data.frame(X[valid_idx, selected_terms, drop = FALSE])
    
    safe_names <- make.names(colnames(X_train), unique = TRUE)
    colnames(X_train) <- safe_names
    colnames(X_valid) <- safe_names
    
    df_train <- cbind(y = y[train_idx], X_train)
    df_valid <- cbind(y = y[valid_idx], X_valid)
    
    fit <- tryCatch(
      suppressWarnings(glm(y ~ ., data = df_train, family = binomial())),
      error = function(e) NULL
    )
    
    if (is.null(fit)) return(NA_real_)
    
    pred_k <- tryCatch(
      suppressWarnings(predict(fit, newdata = df_valid, type = "response")),
      error = function(e) rep(NA_real_, length(valid_idx))
    )
    
    pred_all[valid_idx] <- pred_k
  }
  
  if (anyNA(pred_all)) return(NA_real_)
  safe_auc(y, pred_all)
}

# =================================================
# Data preprocessing and design matrix preparation
# =================================================
prepare_xy <- function(train_sub, test_sub, preds) {
  x_tr <- train_sub[, preds, drop = FALSE]
  x_te <- test_sub[, preds, drop = FALSE]
  
  temp_vars <- grep("^temp", names(x_tr), value = TRUE)
  for (v in temp_vars) {
    if (is.factor(x_tr[[v]])) {
      x_tr[[v]] <- as.numeric(as.character(x_tr[[v]]))
      x_te[[v]] <- as.numeric(as.character(x_te[[v]]))
    } else if (is.character(x_tr[[v]])) {
      x_tr[[v]] <- as.numeric(x_tr[[v]])
      x_te[[v]] <- as.numeric(x_te[[v]])
    }
  }
  
  char_vars <- names(x_tr)[sapply(x_tr, is.character)]
  for (v in char_vars) {
    x_tr[[v]] <- factor(x_tr[[v]])
    x_te[[v]] <- factor(x_te[[v]], levels = levels(x_tr[[v]]))
  }
  
  factor_vars <- names(x_tr)[sapply(x_tr, is.factor)]
  for (v in factor_vars) {
    x_te[[v]] <- factor(x_te[[v]], levels = levels(x_tr[[v]]))
  }
  
  num_vars <- preds[sapply(x_tr[preds], is.numeric)]
  for (v in num_vars) {
    mu <- mean(x_tr[[v]], na.rm = TRUE)
    sdv <- sd(x_tr[[v]], na.rm = TRUE)
    if (is.na(sdv) || sdv == 0) {
      x_tr[[v]] <- 0
      x_te[[v]] <- 0
    } else {
      x_tr[[v]] <- (x_tr[[v]] - mu) / sdv
      x_te[[v]] <- (x_te[[v]] - mu) / sdv
    }
  }
  
  X_tr <- sparse.model.matrix(~ . - 1, data = x_tr)
  X_te <- sparse.model.matrix(~ . - 1, data = x_te)
  
  miss <- setdiff(colnames(X_tr), colnames(X_te))
  if (length(miss) > 0) {
    X_te <- cbind(
      X_te,
      Matrix(
        0,
        nrow = nrow(X_te),
        ncol = length(miss),
        sparse = TRUE,
        dimnames = list(NULL, miss)
      )
    )
  }
  
  X_te <- X_te[, colnames(X_tr), drop = FALSE]
  
  list(
    X_tr_mat = as.matrix(X_tr),
    X_te_mat = as.matrix(X_te)
  )
}

# ===========================================
# Run stability selection within a sex stratum
# ===========================================
run_stratified_stability <- function(train_sub, test_sub, sex_label, preds, confounders,
                                     outcome, n_boot, nfolds_cv, n_cores, base_out_dir) {
  
  out_dir <- file.path(base_out_dir, tolower(sex_label))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  y_tr <- as.integer(train_sub[[outcome]])
  y_te <- as.integer(test_sub[[outcome]])
  
  xy <- prepare_xy(train_sub, test_sub, preds)
  X_tr_mat <- xy$X_tr_mat
  X_te_mat <- xy$X_te_mat
  
  pf <- rep(1, ncol(X_tr_mat))
  names(pf) <- colnames(X_tr_mat)
  
  for (v in confounders) {
    pf[grep(paste0("^", v), colnames(X_tr_mat))] <- 0
  }
  
  always_in_terms <- names(pf)[pf == 0]
  penalized_terms <- names(pf)[pf > 0]
  
  cat("\n==============================\n")
  cat("Running sex-stratified model:", sex_label, "\n")
  cat("==============================\n")
  cat("N train =", nrow(train_sub), "| N test =", nrow(test_sub), "\n")
  cat("Cases in train =", sum(y_tr, na.rm = TRUE), "| Cases in test =", sum(y_te, na.rm = TRUE), "\n")
  cat("Total design-matrix columns:", ncol(X_tr_mat), "\n")
  cat("Unpenalized terms:", length(always_in_terms), "\n")
  cat("Penalized terms:", length(penalized_terms), "\n")
  cat("Number of subsamples =", n_boot, "\n")
  cat("Sampling proportion =", 0.5, "\n")
  cat("Cores passed to sharp =", n_cores, "\n")
  flush.console()
  
  set.seed(2026)
  
  stab <- VariableSelection(
    xdata = X_tr_mat,
    ydata = y_tr,
    family = "binomial",
    implementation = PenalisedRegression,
    penalty.factor = pf,
    alpha = 1,
    standardize = FALSE,
    verbose = TRUE,
    K = n_boot,
    resampling = "subsampling",
    tau = 0.5,
    n_cores = n_cores,
    n_cat = 3
  )
  
  saveRDS(
    stab,
    file.path(out_dir, paste0("stab_subsample_", tolower(sex_label), ".rds"))
  )
  
  best_par <- Argmax(stab)
  best_lambda <- best_par[1, 1]
  best_pi     <- best_par[1, 2]
  
  selected_vec <- SelectedVariables(stab)
  selected_terms_all <- names(selected_vec)[as.vector(selected_vec) == 1]
  
  best_selected_penalized <- setdiff(selected_terms_all, always_in_terms)
  final_terms <- c(always_in_terms, best_selected_penalized)
  
  sel_prop <- SelectionProportions(stab)
  
  sel_prop_df <- data.frame(
    term = names(sel_prop),
    selection_proportion = as.numeric(sel_prop),
    stringsAsFactors = FALSE
  ) |>
    arrange(desc(selection_proportion))
  
  write.csv(
    sel_prop_df,
    file.path(out_dir, paste0("selection_proportions_all_variables_", tolower(sex_label), ".csv")),
    row.names = FALSE
  )
  
  sel_prop_best_pen <- sel_prop_df |>
    filter(term %in% penalized_terms) |>
    mutate(
      selected = selection_proportion >= best_pi,
      term_plot = make_pretty_names(term),
      sex_group = sex_label
    ) |>
    arrange(desc(selection_proportion))
  
  write.csv(
    sel_prop_best_pen,
    file.path(out_dir, paste0("selection_proportions_penalized_", tolower(sex_label), ".csv")),
    row.names = FALSE
  )
  
  set.seed(999)
  fold_id <- sample(rep(seq_len(nfolds_cv), length.out = length(y_tr)))
  
  final_cv_auc <- cv_auc_from_selected(
    X = X_tr_mat,
    y = y_tr,
    selected_terms = final_terms,
    fold_id = fold_id
  )
  
  final_test_auc <- NA_real_
  test_prob <- NULL
  
  X_tr_sel <- as.data.frame(X_tr_mat[, final_terms, drop = FALSE])
  X_te_sel <- as.data.frame(X_te_mat[, final_terms, drop = FALSE])
  
  safe_names <- make.names(colnames(X_tr_sel), unique = TRUE)
  colnames(X_tr_sel) <- safe_names
  colnames(X_te_sel) <- safe_names
  
  final_fit <- tryCatch(
    suppressWarnings(glm(
      y_tr ~ .,
      data = cbind(y_tr = y_tr, X_tr_sel),
      family = binomial()
    )),
    error = function(e) NULL
  )
  
  if (!is.null(final_fit)) {
    test_prob <- tryCatch(
      suppressWarnings(predict(final_fit, newdata = X_te_sel, type = "response")),
      error = function(e) rep(NA_real_, nrow(X_te_sel))
    )
    
    if (!anyNA(test_prob)) {
      final_test_auc <- safe_auc(y_te, test_prob)
    }
  }
  
  png(
    filename = file.path(out_dir, paste0("tuning_heatmap_sharp_subsample_", tolower(sex_label), ".png")),
    width = 12,
    height = 7.5,
    units = "in",
    res = 320
  )
  par(mar = c(6, 5, 4, 2) + 0.1)
  CalibrationPlot(stab, xlab = expression(lambda), ylab = expression(pi))
  dev.off()
  
  plot_df <- sel_prop_best_pen |>
    arrange(desc(selection_proportion))
  
  plot_df$term_plot <- factor(plot_df$term_plot, levels = plot_df$term_plot)
  
  p_sel <- ggplot(plot_df, aes(x = term_plot, y = selection_proportion, color = selected)) +
    geom_segment(
      aes(x = term_plot, xend = term_plot, y = 0, yend = selection_proportion),
      linewidth = 0.9
    ) +
    geom_point(size = 2.2) +
    geom_hline(
      yintercept = best_pi,
      linetype = "dashed",
      linewidth = 0.8,
      color = "#C44E52"
    ) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey75")) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      title = paste0("Selection proportions: ", sex_label),
      subtitle = paste0(
        "Best lambda = ", format(best_lambda, scientific = TRUE),
        "; best pi = ", sprintf("%.2f", best_pi),
        "; penalized selected = ", length(best_selected_penalized),
        "; total terms in final model = ", length(final_terms)
      ),
      x = NULL,
      y = "Selection proportion"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 18),
      plot.subtitle = element_text(size = 11),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.position = "none",
      panel.grid.minor = element_blank()
    )
  
  ggsave(
    file.path(out_dir, paste0("selection_proportion_subsample_", tolower(sex_label), ".png")),
    p_sel,
    width = 14,
    height = 8.5,
    dpi = 320
  )
  
  if (!is.null(test_prob) && !anyNA(test_prob)) {
    roc_obj <- roc(y_te, test_prob, quiet = TRUE)
    
    roc_df <- data.frame(
      specificity = roc_obj$specificities,
      sensitivity = roc_obj$sensitivities
    ) |>
      mutate(fpr = 1 - specificity)
    
    p_roc <- ggplot(roc_df, aes(x = fpr, y = sensitivity)) +
      geom_path(linewidth = 1.2, lineend = "round") +
      geom_abline(
        slope = 1, intercept = 0,
        linetype = "dashed", linewidth = 0.8, color = "grey65"
      ) +
      coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
      scale_x_continuous(breaks = seq(0, 1, 0.2)) +
      scale_y_continuous(breaks = seq(0, 1, 0.2)) +
      labs(
        title = paste0("ROC curve on test set: ", sex_label),
        subtitle = paste0("AUC = ", sprintf("%.3f", final_test_auc)),
        x = "False positive rate",
        y = "True positive rate"
      ) +
      theme_classic(base_size = 13) +
      theme(
        plot.title = element_text(face = "bold", size = 16, hjust = 0),
        plot.subtitle = element_text(size = 11, color = "grey30"),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        axis.line = element_line(linewidth = 0.6),
        plot.margin = margin(12, 12, 12, 12)
      )
    
    ggsave(
      file.path(out_dir, paste0("roc_test_subsample_", tolower(sex_label), ".png")),
      p_roc,
      width = 6.5,
      height = 6.5,
      dpi = 320
    )
  }
  
  summary_lines <- c(
    paste0("Sex stratum: ", sex_label),
    "Method: sharp stability selection with subsampling",
    paste0("Outcome: ", outcome),
    paste0("Number of subsamples: ", n_boot),
    paste0("Sampling proportion: 0.5"),
    paste0("Number of CV folds: ", nfolds_cv),
    paste0("Best lambda value: ", format(best_lambda, scientific = TRUE)),
    paste0("Best pi threshold: ", sprintf("%.2f", best_pi)),
    paste0("Number of selected penalized variables: ", length(best_selected_penalized)),
    paste0("Final refit CV AUC: ", round(final_cv_auc, 4)),
    paste0("Final refit test AUC: ", round(final_test_auc, 4))
  )
  
  writeLines(
    summary_lines,
    file.path(out_dir, paste0("summary_", tolower(sex_label), ".txt"))
  )
  
  list(
    sex_label = sex_label,
    sel_prop_best_pen = sel_prop_best_pen,
    best_pi = best_pi,
    best_lambda = best_lambda,
    best_selected_penalized = best_selected_penalized,
    final_terms = final_terms,
    final_cv_auc = final_cv_auc,
    final_test_auc = final_test_auc
  )
}

# =========================
# Sex stratification
# =========================
train_female <- train_df[train_df$sex == "Female", , drop = FALSE]
train_male   <- train_df[train_df$sex == "Male", , drop = FALSE]

test_female  <- test_df[test_df$sex == "Female", , drop = FALSE]
test_male    <- test_df[test_df$sex == "Male", , drop = FALSE]

# =========================
# Run female model
# =========================
res_female <- run_stratified_stability(
  train_sub = train_female,
  test_sub = test_female,
  sex_label = "Female",
  preds = preds,
  confounders = confounders,
  outcome = outcome,
  n_boot = n_boot,
  nfolds_cv = nfolds_cv,
  n_cores = n_cores,
  base_out_dir = base_out_dir
)

# =========================
# Run male model
# =========================
res_male <- run_stratified_stability(
  train_sub = train_male,
  test_sub = test_male,
  sex_label = "Male",
  preds = preds,
  confounders = confounders,
  outcome = outcome,
  n_boot = n_boot,
  nfolds_cv = nfolds_cv,
  n_cores = n_cores,
  base_out_dir = base_out_dir
)

# =========================
# Combine female/male selection proportions
# =========================
plot_compare <- full_join(
  res_female$sel_prop_best_pen |>
    select(term, term_plot, selection_proportion, selected) |>
    rename(
      selection_female = selection_proportion,
      selected_female = selected
    ),
  res_male$sel_prop_best_pen |>
    select(term, term_plot, selection_proportion, selected) |>
    rename(
      selection_male = selection_proportion,
      selected_male = selected
    ),
  by = c("term", "term_plot")
) |>
  mutate(
    selection_female = ifelse(is.na(selection_female), 0, selection_female),
    selection_male   = ifelse(is.na(selection_male), 0, selection_male),
    selected_female  = ifelse(is.na(selected_female), FALSE, selected_female),
    selected_male    = ifelse(is.na(selected_male), FALSE, selected_male),
    max_selection    = pmax(selection_female, selection_male)
  ) |>
  filter(!(selection_female == 0 & selection_male == 0)) |>
  arrange(desc(max_selection))

write.csv(
  plot_compare,
  file.path(base_out_dir, "selection_proportions_female_vs_male.csv"),
  row.names = FALSE
)

combined_summary <- c(
  "Sex-stratified sharp stability selection with subsampling",
  paste0("Outcome: ", outcome),
  paste0("Number of subsamples: ", n_boot),
  paste0("Sampling proportion: 0.5"),
  paste0("Number of CV folds: ", nfolds_cv),
  paste0("Female best lambda: ", format(res_female$best_lambda, scientific = TRUE)),
  paste0("Female best pi: ", sprintf("%.2f", res_female$best_pi)),
  paste0("Female selected penalized variables: ", length(res_female$best_selected_penalized)),
  paste0("Female test AUC: ", round(res_female$final_test_auc, 4)),
  paste0("Male best lambda: ", format(res_male$best_lambda, scientific = TRUE)),
  paste0("Male best pi: ", sprintf("%.2f", res_male$best_pi)),
  paste0("Male selected penalized variables: ", length(res_male$best_selected_penalized)),
  paste0("Male test AUC: ", round(res_male$final_test_auc, 4))
)

writeLines(
  combined_summary,
  file.path(base_out_dir, "summary_female_vs_male.txt")
)