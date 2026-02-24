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

