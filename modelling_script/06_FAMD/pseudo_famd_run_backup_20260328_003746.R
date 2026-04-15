# analysis/FAMD/pseudo_famd_run.R
# ------------------------------------------------------------
# PURPOSE
# Pseudo-FAMD / mixed-data PCA for TWO datasets:
#   1) final imputed: ukb_G1_imputed_final.rds
#   2) train imputed: split_imputed_data/ukb_G1_train_imputed.rds
#
# For each dataset, run 3 versions:
#   - pseudo_famd1a: baseline main, no alcohol
#   - pseudo_famd1b: baseline + alcohol sensitivity
#   - pseudo_famd2: biological supplementary
#
# For each version, export 5 plots:
#   1. scree
#   2. individuals by outcome
#   3. individuals by sex
#   4. contribution PC1
#   5. contribution PC2
#
# Total output = 2 datasets x 3 versions x 5 plots = 30 PNGs
#
# IMPORTANT
# - Scatter axes are labelled as PC1 (xx.x%), PC2 (yy.y%)
# - Output folders are overwritten each run
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
})

# -----------------------------
# Paths
# -----------------------------
PATH_FINAL <- "ukb_G1_imputed_final.rds"
PATH_TRAIN <- file.path("split_imputed_data", "ukb_G1_train_imputed.rds")

OUTROOT <- file.path("analysis", "FAMD", "FAMD_output")

# -----------------------------
# Colouring variables
# -----------------------------
COLOR_BY_OUTCOME <- "cvd_incident"
COLOR_BY_SEX     <- "sex"

# -----------------------------
# Variable sets
# -----------------------------
# Version 1a: baseline main, no alcohol
v1a_cont <- c(
  "age", "bmi", "sbp_mean", "dbp_mean", "work_hours_week_clean"
)

v1a_cat <- c(
  "sex", "ethnicity_5cat", "imd_quintile", "urban_rural_2cat",
  "education_2", "employment_2cat", "living_with_partner",
  "diabetes_bin", "famhx_cvd", "smoking_3cat"
)

# Version 1b: baseline + alcohol
v1b_cont <- v1a_cont
v1b_cat  <- c(v1a_cat, "alcohol_combined")

# Version 2: biological supplementary
v2_cont <- c(
  "age", "bmi",
  "crp", "hba1c", "glucose",
  "hdl_cholesterol", "ldl_cholesterol", "total_triglyceride",
  "creatinine", "blood_vitamin_d"
)

v2_cat <- c(
  "sex", "diabetes_bin"
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

save_png <- function(p, outpath, width = 8, height = 6, dpi = 300) {
  dir.create(dirname(outpath), recursive = TRUE, showWarnings = FALSE)
  ggsave(outpath, plot = p, width = width, height = height, dpi = dpi)
  message("Saved: ", outpath)
}

clean_group_vector <- function(x, fallback_name = "group") {
  if (is.logical(x)) {
    x <- factor(x, levels = c(FALSE, TRUE), labels = c("No", "Yes"))
  } else {
    x <- as.factor(x)
  }
  x <- droplevels(x)
  x
}

# Build mixed matrix:
# - continuous: numeric + z-score
# - categorical: factor + one-hot
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
  mm_cols <- character(0)
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
    mm_cols <- colnames(Xf)
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
  
  list(
    X = X,
    ok = ok,
    cont_vars = cont_vars,
    cat_vars = cat_vars,
    mm_cols = mm_cols
  )
}

aggregate_contrib <- function(loadings_vec, cont_vars, cat_vars, mm_colnames) {
  contrib <- list()
  
  # continuous: squared loading
  if (length(cont_vars) > 0) {
    for (v in cont_vars) {
      if (v %in% names(loadings_vec)) {
        contrib[[v]] <- loadings_vec[[v]]^2
      }
    }
  }
  
  # categorical: sum squared loadings across dummy columns starting with var name
  if (length(cat_vars) > 0 && length(mm_colnames) > 0) {
    for (v in cat_vars) {
      idx <- startsWith(mm_colnames, v)
      cols <- mm_colnames[idx]
      cols <- intersect(cols, names(loadings_vec))
      if (length(cols) > 0) {
        contrib[[v]] <- sum(loadings_vec[cols]^2)
      }
    }
  }
  
  out <- data.frame(
    variable = names(contrib),
    contrib  = as.numeric(unlist(contrib)),
    stringsAsFactors = FALSE
  )
  
  out[order(out$contrib, decreasing = TRUE), , drop = FALSE]
}

plot_scree <- function(pca, prefix) {
  var_expl <- (pca$sdev^2) / sum(pca$sdev^2)
  dfv <- data.frame(
    PC = seq_along(var_expl),
    explained = 100 * var_expl
  )
  
  ggplot(dfv[1:min(10, nrow(dfv)), ], aes(x = PC, y = explained)) +
    geom_col() +
    scale_x_continuous(breaks = dfv$PC[1:min(10, nrow(dfv))]) +
    labs(
      title = paste0(prefix, " - Scree Plot"),
      x = "Principal Component",
      y = "Variance Explained (%)"
    ) +
    theme_bw()
}

plot_individuals <- function(scores, color_vec, prefix, label, pc1_pct, pc2_pct) {
  dfi <- data.frame(
    PC1 = scores[, 1],
    PC2 = scores[, 2],
    group = clean_group_vector(color_vec, label)
  )
  
  ggplot(dfi, aes(x = PC1, y = PC2, color = group)) +
    geom_point(alpha = 0.25, size = 0.7) +
    labs(
      title = paste0(prefix, " - Individuals (colored by ", label, ")"),
      x = sprintf("PC1 (%.1f%%)", pc1_pct),
      y = sprintf("PC2 (%.1f%%)", pc2_pct),
      color = label
    ) +
    theme_bw()
}

plot_top_contrib <- function(contrib_df, prefix, pc_label, top_n = 20) {
  top <- head(contrib_df, top_n)
  top$variable <- factor(top$variable, levels = rev(top$variable))
  
  ggplot(top, aes(x = variable, y = contrib)) +
    geom_col() +
    coord_flip() +
    labs(
      title = paste0(prefix, " - Top Contributions to ", pc_label),
      x = "",
      y = "Sum of squared loadings"
    ) +
    theme_bw()
}

run_one_version <- function(df, cont_vars, cat_vars, prefix, outdir_dataset) {
  vars_needed <- unique(c(cont_vars, cat_vars, COLOR_BY_OUTCOME, COLOR_BY_SEX))
  assert_has_vars(df, vars_needed, context = prefix)
  
  prep <- prep_mixed_matrix(df, cont_vars, cat_vars)
  X <- prep$X
  ok <- prep$ok
  
  # PCA on prepared mixed matrix
  pca <- prcomp(X, center = FALSE, scale. = FALSE)
  
  var_expl <- 100 * (pca$sdev^2) / sum(pca$sdev^2)
  pc1_pct <- var_expl[1]
  pc2_pct <- var_expl[2]
  
  scores <- pca$x
  outcome_vec <- df[[COLOR_BY_OUTCOME]][ok]
  sex_vec     <- df[[COLOR_BY_SEX]][ok]
  
  # 1) scree
  save_png(
    plot_scree(pca, prefix),
    file.path(outdir_dataset, paste0(prefix, "_scree.png"))
  )
  
  # 2) outcome scatter
  save_png(
    plot_individuals(scores, outcome_vec, prefix, COLOR_BY_OUTCOME, pc1_pct, pc2_pct),
    file.path(outdir_dataset, paste0(prefix, "_indiv_by_outcome.png"))
  )
  
  # 3) sex scatter
  save_png(
    plot_individuals(scores, sex_vec, prefix, COLOR_BY_SEX, pc1_pct, pc2_pct),
    file.path(outdir_dataset, paste0(prefix, "_indiv_by_sex.png"))
  )
  
  # 4-5) contributions PC1 / PC2
  load1 <- pca$rotation[, 1]
  load2 <- pca$rotation[, 2]
  
  c1 <- aggregate_contrib(load1, prep$cont_vars, prep$cat_vars, prep$mm_cols)
  c2 <- aggregate_contrib(load2, prep$cont_vars, prep$cat_vars, prep$mm_cols)
  
  save_png(
    plot_top_contrib(c1, prefix, "PC1", top_n = 20),
    file.path(outdir_dataset, paste0(prefix, "_contrib_pc1.png")),
    width = 9, height = 6
  )
  
  save_png(
    plot_top_contrib(c2, prefix, "PC2", top_n = 20),
    file.path(outdir_dataset, paste0(prefix, "_contrib_pc2.png")),
    width = 9, height = 6
  )
  
  invisible(pca)
}

run_dataset <- function(inpath, dataset_name) {
  message("\n==============================")
  message("Reading dataset: ", dataset_name)
  message("Path: ", inpath)
  message("==============================")
  
  df <- readRDS(inpath)
  
  outdir_dataset <- file.path(OUTROOT, dataset_name)
  
  # overwrite this dataset folder on each run
  if (dir.exists(outdir_dataset)) {
    unlink(outdir_dataset, recursive = TRUE, force = TRUE)
  }
  dir.create(outdir_dataset, recursive = TRUE, showWarnings = FALSE)
  
  # 1a
  run_one_version(
    df = df,
    cont_vars = v1a_cont,
    cat_vars  = v1a_cat,
    prefix    = "pseudo_famd1a_baseline_main_no_alcohol",
    outdir_dataset = outdir_dataset
  )
  
  # 1b
  run_one_version(
    df = df,
    cont_vars = v1b_cont,
    cat_vars  = v1b_cat,
    prefix    = "pseudo_famd1b_baseline_plus_alcohol_sensitivity",
    outdir_dataset = outdir_dataset
  )
  
  # 2
  run_one_version(
    df = df,
    cont_vars = v2_cont,
    cat_vars  = v2_cat,
    prefix    = "pseudo_famd2_biological_supplementary",
    outdir_dataset = outdir_dataset
  )
  
  message("Done dataset: ", dataset_name, " -> ", outdir_dataset)
}

# -----------------------------
# Main
# -----------------------------
dir.create(OUTROOT, recursive = TRUE, showWarnings = FALSE)

run_dataset(PATH_FINAL, "final")
run_dataset(PATH_TRAIN, "train")

message("\nALL DONE. Outputs under: ", OUTROOT)