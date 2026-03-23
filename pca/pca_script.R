# =============================================================================
# PCA on Continuous Variables 
# =============================================================================
# Run AFTER: imputation/ukb_G1_imputed.rds has been created
# Input:  imputation/ukb_G1_imputed_final.rds
# Outputs:
#   pca/pca_result.rds            — full PCA object
#   pca/pca_scaled_data.rds       — standardised data matrix
#   pca/pca_scree_plot.png        — variance explained per PC
#   pca/pca_loadings_heatmap.png  — clustered loadings heatmap
# =============================================================================

library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)
library(dplyr)
library(ggplot2)

setwd("/rds/general/project/hda_25-26/live/TDS/TDS_Group1")

out_dir <- "pca"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
ukb <- readRDS("ukb_G1_imputed_final.rds")

# -----------------------------------------------------------------------------
# Auto-detect continuous variables
# -----------------------------------------------------------------------------
# Automatically selects numeric columns, then excludes:
#   - ID / date / outcome variables (same protected list as imputation script)
#   - Binary 0/1 flags (encoded categoricals, not truly continuous)
#   - Derived group/category variables
#   - Known redundant variables (yob when age is present; raw BP when means exist)

id_like      <- grep("^eid$|id$",            names(ukb), ignore.case = TRUE, value = TRUE)
date_like    <- grep("date|death|time|dod",  names(ukb), ignore.case = TRUE, value = TRUE)
outcome_like <- grep("^cvd",                 names(ukb), ignore.case = TRUE, value = TRUE)
group_like   <- grep("_group$|_cat$|_2cat$|_5cat$|_status$|_quintile$",
                     names(ukb), ignore.case = TRUE, value = TRUE)

# Variables to manually exclude (redundant or not interpretable as continuous)
manual_exclude <- c(
  "yob",                  # redundant with age
  "sys_bp_automatic",     # replaced by sbp_mean
  "sbp_manual",           # replaced by sbp_mean
  "dys_bp_automatic",     # replaced by dbp_mean
  "dbp_manual",           # replaced by dbp_mean
  "imd_source",           # character label for which nation IMD came from
  "imd_final"             # raw scores not comparable across nations; use imd_quintile as covariate instead
)

protected <- unique(c(id_like, date_like, outcome_like, group_like, manual_exclude))

# Select numeric columns not in the protected list, drop binary vars
all_numeric     <- names(ukb)[sapply(ukb, is.numeric)]
candidate_vars  <- setdiff(all_numeric, protected)

is_binary <- function(x) {
  vals <- unique(na.omit(x))
  length(vals) <= 2 && all(vals %in% c(0, 1))
}

binary_vars     <- candidate_vars[sapply(candidate_vars, function(v) is_binary(ukb[[v]]))]
continuous_vars <- setdiff(candidate_vars, binary_vars)

cat("Variables excluded as binary flags:", length(binary_vars), "\n")
cat(paste(binary_vars, collapse = ", "), "\n\n")
cat("Continuous variables selected for PCA:", length(continuous_vars), "\n")
cat(paste(continuous_vars, collapse = ", "), "\n")

# -----------------------------------------------------------------------------
# Continuous variable data frame
# -----------------------------------------------------------------------------

pca_data <- ukb %>%
  select(all_of(continuous_vars)) %>%
  mutate(across(everything(), as.numeric))

cat("\nDimensions of PCA data frame:", nrow(pca_data), "x", ncol(pca_data), "\n")

# Confirm no missing values
n_missing <- sum(is.na(pca_data))
if (n_missing > 0) {
  warning(n_missing, " missing values remain. Dropping incomplete rows.")
  pca_data <- na.omit(pca_data)
} else {
  cat("No missing values — proceeding with all", nrow(pca_data), "participants.\n")
}

# -----------------------------------------------------------------------------
# Standardise df
# -----------------------------------------------------------------------------

pca_scaled <- scale(pca_data)

cat("\nVerification — means ~0, SDs ~1:\n")
cat("Mean of column means:", round(mean(colMeans(pca_scaled)), 6), "\n")
cat("Mean of column SDs:  ", round(mean(apply(pca_scaled, 2, sd)), 6), "\n")

saveRDS(as.data.frame(pca_scaled), file.path(out_dir, "pca_scaled_data.rds"))

# -----------------------------------------------------------------------------
# PCA
# -----------------------------------------------------------------------------

pca_result <- prcomp(pca_scaled, center = FALSE, scale. = FALSE)
saveRDS(pca_result, file.path(out_dir, "pca_result.rds"))

# -----------------------------------------------------------------------------
# Scree plot
# -----------------------------------------------------------------------------

var_pct <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100
scree_df <- data.frame(
  PC      = seq_along(var_pct),
  VarPct  = var_pct
)

scree_plot <- ggplot(scree_df[1:20, ], aes(x = PC, y = VarPct)) +
  geom_col(fill = "#2C7BB6", colour = "white") +
  geom_text(aes(label = sprintf("%.1f%%", VarPct)), vjust = -0.4, size = 3) +
  scale_x_continuous(breaks = 1:20) +
  labs(title = "Scree Plot — Variance Explained by Each PC",
       x = "Principal Component", y = "% Variance Explained") +
  theme_bw()

ggsave(file.path(out_dir, "pca_scree_plot.png"),
       scree_plot, width = 9, height = 5, dpi = 120)

# Print variance explained
var_explained <- summary(pca_result)$importance
cat("\nVariance explained by first 10 PCs:\n")
print(round(var_explained[, 1:10], 3))

cumvar  <- var_explained["Cumulative Proportion", ]
n_pc_70 <- which(cumvar >= 0.70)[1]
n_pc_80 <- which(cumvar >= 0.80)[1]
cat("\nPCs needed for 70% cumulative variance:", n_pc_70, "\n")
cat("PCs needed for 80% cumulative variance:", n_pc_80, "\n")

# -----------------------------------------------------------------
# Examine loadings + heatmap for PC1-10
# --------------------------------
loadings_df <- as.data.frame(pca_result$rotation[, 1:n_pcs])
for (pc in paste0("PC", 1:n_pcs)) {
  cat("\n===", pc, "===\n")
  abs_load <- abs(loadings[[pc]])
  names(abs_load) <- rownames(loadings)
  top10_names <- names(sort(abs_load, decreasing = TRUE))[1:10]
  result <- loadings[top10_names, pc, drop = FALSE]
  print(round(result, 3))
}

n_pcs <- 10   # number of PCs to display 

loading_mat <- pca_result$rotation[, 1:n_pcs]

# Hierarchical clustering on rows — groups variables that load on the same PCs
row_order <- hclust(dist(loading_mat))$order
var_order  <- rownames(loading_mat)[row_order]

loading_long <- as.data.frame(loading_mat) %>%
  rownames_to_column("Variable") %>%
  pivot_longer(-Variable, names_to = "PC", values_to = "Loading") %>%
  mutate(
    Variable = factor(Variable, levels = var_order),      
    PC       = factor(PC, levels = paste0("PC", 1:n_pcs))  
  )

heatmap_plot <- ggplot(loading_long, aes(x = PC, y = Variable, fill = Loading)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient2(low  = "#D7191C", mid = "#FAFAF0", high = "#2C7BB6",
                       midpoint = 0, limits = c(-1, 1)) +
  labs(title = "PCA Loadings Heatmap", x = NULL, y = NULL) +
  theme_bw(base_size = 10) +
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 9))

ggsave(file.path(out_dir, "pca_loadings_heatmap.png"),
       heatmap_plot, width = 8, height = 12, dpi = 120)

