library(dplyr)
library(glmnet)
library(pROC)
library(ggplot2)
library(Matrix)
library(stringr)
library(forcats)

# =========================
# Set paths and load data
# =========================
setwd("/rds/general/project/hda_25-26/live/TDS/TDS_Group1")

out_dir <- "modelling_script/multivariate_finding/lasso_sensitivity"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

train_df <- readRDS("split_imputed_data/ukb_G1_train_imputed.rds")
test_df  <- readRDS("split_imputed_data/ukb_G1_test_imputed.rds")

outcome <- "cvd_incident"
y_tr <- as.integer(train_df[[outcome]])
y_te <- as.integer(test_df[[outcome]])

# =========================
# Define predictors
# =========================
confounders <- c("age", "sex", "ethnicity_5cat")

external_exposome_vars <- c(
  "tv_group",
  "total_met_group",
  "smoking_3cat",
  "alcohol_combined",
  "salt_3cat",
  "education_2",
  "employment_2cat",
  "imd_quintile",
  "sf_score",
  "redwine_group",
  "sbp_mean",
  "dbp_mean"
)

pattern_regex <- paste(
  c(
    "fruit", "fish", "tea", "coffee",
    "sleep", "stress_group", "neurotic", "satis",
    "no2", "pm10", "pm2_5", "noise",
    "greenspace", "temp",
    "living", "urban", "work_hours"
  ),
  collapse = "|"
)

pattern_vars <- grep(pattern_regex, names(train_df), value = TRUE)

vars <- unique(c(external_exposome_vars, pattern_vars))
vars <- vars[vars %in% names(train_df) & vars %in% names(test_df)]

preds <- c(confounders, vars)

cat("Predictors used:", length(preds), "\n")

# =========================
# Prepare predictors and standardize numeric variables
# =========================
x_tr <- train_df[, preds, drop = FALSE]
x_te <- test_df[, preds, drop = FALSE]

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

if (length(temp_vars) > 0) {
  cat("Temp variables checked:", paste(temp_vars, collapse = ", "), "\n")
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

# =========================
# Create design matrices
# =========================
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

# =========================
# Set penalty factors
# =========================
pf <- rep(1, ncol(X_tr))
names(pf) <- colnames(X_tr)

for (v in confounders) {
  pf[grep(paste0("^", v), colnames(X_tr))] <- 0
}

# =========================
# Fit LASSO model
# =========================
set.seed(123)

cvfit <- cv.glmnet(
  x = X_tr,
  y = y_tr,
  family = "binomial",
  alpha = 1,
  nfolds = 10,
  type.measure = "auc",
  penalty.factor = pf
)

lambda_min <- cvfit$lambda.min
lambda_1se <- cvfit$lambda.1se

cat("lambda.min =", lambda_min, "\n")
cat("lambda.1se =", lambda_1se, "\n")

# =========================
# Define result extraction function
# =========================
run_lasso_pipeline <- function(lambda_value, lambda_name) {
  coef_mat <- as.matrix(coef(cvfit, s = lambda_value))
  
  coef_df <- data.frame(
    term = rownames(coef_mat),
    beta = as.numeric(coef_mat[, 1])
  ) |>
    filter(term != "(Intercept)", beta != 0) |>
    mutate(
      OR = exp(beta),
      abs_beta = abs(beta)
    ) |>
    arrange(desc(abs_beta))
  
  write.csv(
    coef_df,
    file.path(out_dir, paste0("lasso_selected_variables_", lambda_name, ".csv")),
    row.names = FALSE
  )
  
  p_te <- predict(cvfit, newx = X_te, s = lambda_value, type = "response")
  roc_obj <- roc(y_te, as.numeric(p_te), quiet = TRUE)
  auc_te <- as.numeric(auc(roc_obj))
  
  roc_df <- data.frame(
    tpr = roc_obj$sensitivities,
    fpr = 1 - roc_obj$specificities
  )
  
  roc_plot <- ggplot(roc_df, aes(fpr, tpr)) +
    geom_area(fill = "#DCE6F2", alpha = 0.8) +
    geom_line(linewidth = 1.6, color = "#2C7FB8") +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      linewidth = 0.9,
      color = "#D95F5F"
    ) +
    annotate(
      "label",
      x = 0.78,
      y = 0.12,
      label = paste0("AUC = ", sprintf("%.3f", auc_te)),
      size = 5,
      fill = "white",
      color = "#1A1A1A",
      label.size = 0.3,
      label.r = unit(0.15, "lines")
    ) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(
      title = paste0("ROC curve (", lambda_name, ")"),
      subtitle = "Test-set discrimination performance",
      x = "False positive rate",
      y = "True positive rate"
    ) +
    coord_fixed(ratio = 1) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 20, color = "#1A1A1A"),
      plot.subtitle = element_text(size = 13, color = "#4D4D4D"),
      axis.text = element_text(size = 11, color = "#222222"),
      axis.title = element_text(size = 13, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey86", linewidth = 0.45),
      plot.margin = margin(15, 15, 15, 15)
    )
  
  ggsave(
    file.path(out_dir, paste0("roc_curve_", lambda_name, ".png")),
    roc_plot,
    width = 7.5,
    height = 7.5,
    dpi = 320
  )
  
  if (nrow(coef_df) == 0) {
    writeLines(
      c(
        paste0("Outcome: ", outcome),
        paste0("Lambda setting: ", lambda_name),
        paste0("Test AUC: ", round(auc_te, 4)),
        paste0("Lambda value: ", format(lambda_value, scientific = TRUE)),
        "Selected predictors: 0",
        "No non-zero predictors selected."
      ),
      file.path(out_dir, paste0("model_summary_", lambda_name, ".txt"))
    )
    
    return(
      list(
        coef_df = coef_df,
        forest_df = NULL,
        auc_te = auc_te
      )
    )
  }
  
  X_sel <- as.data.frame(as.matrix(X_tr[, coef_df$term, drop = FALSE]))
  
  safe_names <- make.names(colnames(X_sel), unique = TRUE)
  
  name_map <- data.frame(
    term_safe = safe_names,
    term = colnames(X_sel),
    stringsAsFactors = FALSE
  )
  
  colnames(X_sel) <- safe_names
  
  refit_df <- cbind(y = y_tr, X_sel)
  refit <- glm(y ~ ., data = refit_df, family = binomial())
  
  s <- summary(refit)$coef
  
  forest_df <- data.frame(
    term_safe = rownames(s)[-1],
    beta = s[-1, 1],
    se = s[-1, 2],
    stringsAsFactors = FALSE
  ) |>
    left_join(name_map, by = "term_safe") |>
    mutate(
      OR = exp(beta),
      low = exp(beta - 1.96 * se),
      high = exp(beta + 1.96 * se),
      label = sprintf("%.2f (%.2f, %.2f)", OR, low, high)
    )
  
  forest_df$term_plot <- forest_df$term |>
    str_replace_all("_", " ") |>
    str_squish()
  
  forest_df$term_plot <- forest_df$term_plot |>
    str_replace_all("pm2 5", "PM2.5") |>
    str_replace_all("pm10", "PM10") |>
    str_replace_all("no2", "NO2") |>
    str_replace_all("imd", "IMD") |>
    str_replace_all("bmi", "BMI") |>
    str_replace_all("tv", "TV") |>
    str_replace_all("cvd", "CVD") |>
    str_replace_all("sbp", "SBP") |>
    str_replace_all("dbp", "DBP")
  
  forest_df$penalty_group <- ifelse(
    grepl("^age$", forest_df$term) |
      grepl("^sex", forest_df$term) |
      grepl("^ethnicity_5cat", forest_df$term),
    "Unpenalized",
    "Penalized"
  )
  
  forest_df <- forest_df |>
    arrange(desc(abs(log(OR))))
  
  forest_df$term_plot <- factor(
    forest_df$term_plot,
    levels = rev(forest_df$term_plot)
  )
  
  x_min <- min(forest_df$low, na.rm = TRUE) * 0.9
  x_max <- max(forest_df$high, na.rm = TRUE) * 2.6
  text_x <- max(forest_df$high, na.rm = TRUE) * 1.35
  
  forest_plot <- ggplot(
    forest_df,
    aes(x = OR, y = term_plot,
        shape = penalty_group,
        color = penalty_group)
  ) +
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      linewidth = 0.7,
      color = "grey55"
    ) +
    geom_errorbarh(
      aes(xmin = low, xmax = high),
      height = 0.16,
      linewidth = 0.8,
      color = "grey45"
    ) +
    geom_point(
      size = 4,
      stroke = 1.2
    ) +
    geom_text(
      aes(x = text_x, label = label),
      hjust = 0,
      size = 3.7,
      color = "#333333"
    ) +
    scale_x_log10() +
    scale_color_manual(
      values = c(
        "Penalized" = "#2C7FB8",
        "Unpenalized" = "#E6550D"
      )
    ) +
    scale_shape_manual(
      values = c(
        "Penalized" = 16,
        "Unpenalized" = 17
      )
    ) +
    labs(
      title = paste0("LASSO-selected predictors (", lambda_name, ")"),
      subtitle = "Odds ratios and 95% confidence intervals from refitted logistic model",
      x = "Odds ratio (log scale)",
      y = NULL,
      shape = "Variable type",
      color = "Variable type"
    ) +
    coord_cartesian(xlim = c(x_min, x_max), clip = "off") +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "grey88"),
      legend.position = "top",
      plot.margin = margin(15, 160, 15, 15)
    )
  
  ggsave(
    file.path(out_dir, paste0("forest_plot_", lambda_name, ".png")),
    forest_plot,
    width = 12,
    height = max(5.5, 0.52 * nrow(forest_df) + 1.8),
    dpi = 320
  )
  
  write.csv(
    forest_df,
    file.path(out_dir, paste0("forest_plot_data_", lambda_name, ".csv")),
    row.names = FALSE
  )
  
  writeLines(
    c(
      paste0("Outcome: ", outcome),
      paste0("Lambda setting: ", lambda_name),
      paste0("Test AUC: ", round(auc_te, 4)),
      paste0("Lambda value: ", format(lambda_value, scientific = TRUE)),
      paste0("Selected predictors: ", nrow(coef_df))
    ),
    file.path(out_dir, paste0("model_summary_", lambda_name, ".txt"))
  )
  
  return(
    list(
      coef_df = coef_df,
      forest_df = forest_df,
      auc_te = auc_te
    )
  )
}

# =========================
# Run pipeline for lambda.min
# =========================
res_min <- run_lasso_pipeline(lambda_min, "lambda_min")

# =========================
# Run pipeline for lambda.1se
# =========================
res_1se <- run_lasso_pipeline(lambda_1se, "lambda_1se")

# =========================
# Save overall summary
# =========================
writeLines(
  c(
    paste0("Outcome: ", outcome),
    paste0("lambda.min: ", format(lambda_min, scientific = TRUE)),
    paste0("lambda.1se: ", format(lambda_1se, scientific = TRUE)),
    paste0("Test AUC (lambda.min): ", round(res_min$auc_te, 4)),
    paste0("Test AUC (lambda.1se): ", round(res_1se$auc_te, 4)),
    paste0("Selected predictors (lambda.min): ", nrow(res_min$coef_df)),
    paste0("Selected predictors (lambda.1se): ", nrow(res_1se$coef_df))
  ),
  file.path(out_dir, "model_summary_overall.txt")
)

cat("Done.\n")
cat("Outputs saved to:", normalizePath(out_dir), "\n")