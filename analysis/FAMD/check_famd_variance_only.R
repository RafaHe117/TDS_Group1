# analysis/FAMD/check_famd_variance_only.R
suppressPackageStartupMessages({})

# -----------------------------
# Paths
# -----------------------------
PATH_FINAL <- "ukb_G1_imputed_final.rds"
PATH_TRAIN <- file.path("split_imputed_data", "ukb_G1_train_imputed.rds")

# -----------------------------
# Variable sets: pseudo_famd1a
# -----------------------------
v1a_cont <- c(
  "age", "bmi", "sbp_mean", "dbp_mean", "work_hours_week_clean"
)

v1a_cat <- c(
  "sex", "ethnicity_5cat", "imd_quintile", "urban_rural_2cat",
  "education_2", "employment_2cat", "living_with_partner",
  "diabetes_bin", "famhx_cvd", "smoking_3cat"
)

# -----------------------------
# Helpers
# -----------------------------
assert_has_vars <- function(df, vars, context = "") {
  miss <- setdiff(vars, names(df))
  if (length(miss) > 0) {
    stop(sprintf("[%s] Missing variables: %s",
                 context, paste(miss, collapse = ", ")),
         call. = FALSE)
  }
}

prep_mixed_matrix <- function(df, cont_vars, cat_vars) {
  dat <- df[, c(cont_vars, cat_vars), drop = FALSE]
  
  # continuous
  Xc <- NULL
  if (length(cont_vars) > 0) {
    for (v in cont_vars) dat[[v]] <- suppressWarnings(as.numeric(dat[[v]]))
    Xc <- scale(dat[, cont_vars, drop = FALSE])
    Xc <- as.matrix(Xc)
    colnames(Xc) <- cont_vars
  }
  
  # categorical
  Xf <- NULL
  if (length(cat_vars) > 0) {
    for (v in cat_vars) {
      if (is.logical(dat[[v]])) {
        dat[[v]] <- factor(dat[[v]], levels = c(FALSE, TRUE), labels = c("No", "Yes"))
      } else {
        dat[[v]] <- as.factor(dat[[v]])
      }
      dat[[v]] <- droplevels(dat[[v]])
    }
    
    form <- as.formula(paste("~", paste(cat_vars, collapse = " + "), "-1"))
    Xf <- model.matrix(form, data = dat)
  }
  
  if (is.null(Xc) && is.null(Xf)) stop("No variables provided.", call. = FALSE)
  
  if (is.null(Xc)) {
    X <- Xf
  } else if (is.null(Xf)) {
    X <- Xc
  } else {
    X <- cbind(Xc, Xf)
  }
  
  ok <- complete.cases(X)
  X <- X[ok, , drop = FALSE]
  
  X
}

check_one_dataset <- function(path, dataset_name, cont_vars, cat_vars) {
  cat("\n==============================\n")
  cat("Dataset:", dataset_name, "\n")
  cat("Path   :", path, "\n")
  cat("==============================\n")
  
  df <- readRDS(path)
  vars_needed <- unique(c(cont_vars, cat_vars))
  assert_has_vars(df, vars_needed, context = dataset_name)
  
  X <- prep_mixed_matrix(df, cont_vars, cat_vars)
  pca <- prcomp(X, center = FALSE, scale. = FALSE)
  
  var_expl <- 100 * (pca$sdev^2) / sum(pca$sdev^2)
  
  n_show <- min(10, length(var_expl))
  out <- data.frame(
    PC = paste0("PC", seq_len(n_show)),
    Variance_Explained_Percent = round(var_expl[seq_len(n_show)], 6)
  )
  
  print(out, row.names = FALSE)
  
  cat("\nRounded display used in slides:\n")
  cat(sprintf("PC1 = %.1f%%\n", var_expl[1]))
  cat(sprintf("PC2 = %.1f%%\n", var_expl[2]))
  
  cat("\nMore precise values:\n")
  cat(sprintf("PC1 = %.6f%%\n", var_expl[1]))
  cat(sprintf("PC2 = %.6f%%\n", var_expl[2]))
  cat(sprintf("PC3 = %.6f%%\n", var_expl[3]))
}

# -----------------------------
# Run checks
# -----------------------------
check_one_dataset(PATH_FINAL, "final", v1a_cont, v1a_cat)
check_one_dataset(PATH_TRAIN, "train", v1a_cont, v1a_cat)