# Paths
in_dir  <- "/rds/general/project/hda_25-26/live/TDS/TDS_Group1"
out_dir <- file.path(in_dir, "modelling_script")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
setwd(in_dir)

UK1_1 <- readRDS("ukb_G1_imputed_final.rds")

# Keep numeric vars
UK_num <- UK1_1[, vapply(UK1_1, is.numeric, logical(1)), drop = FALSE]
UK_num <- as.data.frame(lapply(UK_num, as.numeric))

# Correlation matrix
corr_matrix <- cor(UK_num, use = "pairwise.complete.obs")

# ---------- Full correlation table ----------
upper_full <- corr_matrix
upper_full[lower.tri(upper_full, diag = TRUE)] <- NA

all_pairs <- na.omit(as.data.frame(as.table(upper_full)))
names(all_pairs) <- c("Variable_1", "Variable_2", "Correlation")
all_pairs <- all_pairs[order(-abs(all_pairs$Correlation)), ]

write.csv(
  all_pairs,
  file = file.path(out_dir, "all_correlations.csv"),
  row.names = FALSE
)

# ---------- Strong correlation table ----------
thr <- 0.7
strong_pairs <- subset(all_pairs, abs(Correlation) >= thr & abs(Correlation) < 1)

write.csv(
  strong_pairs,
  file = file.path(out_dir, "strong_pairs.csv"),
  row.names = FALSE
)

suppressPackageStartupMessages({
  library(reshape2)
  library(ggplot2)
})

# ---------- Full heatmap ----------
full_melt <- melt(corr_matrix)

full_plot <- ggplot(full_melt, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6)
  ) +
  labs(title = "Full Correlation Heatmap", x = "", y = "")

ggsave(
  filename = file.path(out_dir, "full_correlation_heatmap.png"),
  plot = full_plot,
  width = 12,
  height = 12,
  dpi = 300
)

# ---------- Strong heatmap ----------
strong_matrix <- corr_matrix
strong_matrix[abs(strong_matrix) < thr] <- NA

strong_melt <- melt(strong_matrix)

strong_plot <- ggplot(strong_melt, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0, na.value = "white"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6)
  ) +
  labs(title = "Strong Correlations (|r| ≥ 0.7)", x = "", y = "")

ggsave(
  filename = file.path(out_dir, "strong_correlations_heatmap.png"),
  plot = strong_plot,
  width = 12,
  height = 12,
  dpi = 300
)