library(ggplot2)
library(reshape2)

setwd("/rds/general/project/hda_25-26/live/TDS/TDS_Group1")

ukb <- readRDS("ukb_G1_preprocessed.rds")

############################################
# Print variables with high missingness
############################################
print_high_missingness <- function(df, threshold = 0.30) {
  
  na_props <- colMeans(is.na(df))
  
  high_na_cols <- na_props[na_props >= threshold]
  high_na_cols <- sort(high_na_cols, decreasing = TRUE)
  
  result_table <- data.frame(
    Variable = names(high_na_cols),
    Missing_Proportion = unname(high_na_cols)
  )
  
  print(result_table, row.names = FALSE)
}

print_high_missingness(ukb)

############################################
# Produce missingness heatmap
############################################
missingness_heatmap <- function(df, bin_size = 100, fig_name) {
  
  # Boolean transformation
  na_matrix <- is.na(df) * 1
  
  # Define output dimensions and binning
  n_bins <- ceiling(nrow(na_matrix) / bin_size)
  row_indices <- rep(1:n_bins, each = bin_size, length.out = nrow(na_matrix))
  
  # Density calculation
  agg_matrix <- rowsum(na_matrix, row_indices) / bin_size
  
  # Data reshaping
  long_data <- melt(agg_matrix, varnames = c("Row_Bin", "Feature"))
  
  p <- ggplot(long_data, aes(x = Feature, y = Row_Bin, fill = value)) +
    geom_raster() +
    scale_fill_gradient(low = "grey90", high = "darkred", 
                        name = "Missingness Proportion",
                        limits = c(0, 1)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_discrete(expand = c(0,0)) +
    labs(title = "Missingness Density Heatmap",
         x = "Variables",
         y = sprintf("Aggregated Row Bins (1 bin = %d observations)", bin_size)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 9), 
          axis.ticks.x = element_line(color = "grey80"),
          panel.grid = element_blank(),
          legend.title = element_text(margin = margin(b = 15)))
  
  ggsave(filename = file.path("figure", fig_name), plot = p, 
         width = 20, height = 12, dpi = 450)
}

missingness_heatmap(ukb, bin_size = 1000, fig_name = "heatmap_pre_exclusion.png")
missingness_heatmap(ukb_G1_cleaned, bin_size = 1000, fig_name = "heatmap_pre_imputation.png")
