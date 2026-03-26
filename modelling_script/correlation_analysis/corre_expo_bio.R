library(dplyr)

# Paths
in_dir  <- "/rds/general/user/rh725/projects/hda_25-26/live/TDS/TDS_Group1"
out_dir <- file.path(in_dir, "modelling_script")
corr_dir <- file.path(out_dir, "correlation_analysis")

# Read correlation matrix
corr_file <- file.path(corr_dir, "correlation_matrix.csv")

corr_df <- read.csv(
  corr_file,
  row.names = 1,
  check.names = FALSE
)

corr_mat <- as.matrix(corr_df)

# Exposure list
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
  "sf_score"
)

pattern_regex <- paste(
  c(
    "fruit", "fish",
    "tea", "coffee",
    "sleep",
    "stress_group",
    "neurotic", "satis",
    "no2", "pm10", "pm2_5",
    "noise",
    "greenspace",
    "temp",
    "living",
    "urban",
    "work_hours"
  ),
  collapse = "|"
)

pattern_vars <- grep(pattern_regex, rownames(corr_mat), value = TRUE)
exposure_vars <- unique(c(exposome_list, pattern_vars))
exposure_vars <- exposure_vars[exposure_vars %in% rownames(corr_mat)]

# Biomarker list
biomarker_map <- data.frame(
  variable = c(
    "wbc_count","rbc_count","hemoglobin","mcv","platelet_count","mpv",
    "lymphocyte_count","monocyte_count","neutrophil_count",
    "eosinophil_count","basophil_count","nrbc_count",
    "creatinine","cystatin_c","urea","urate","microalbumin",
    "albumin","alkaline_phosphatase",
    "alanine_aminotransferase","aspartate_aminotransferase",
    "direct_bilirubin","total_bilirubin","total_blood_protein",
    "cholesterol","hdl_cholesterol","ldl_cholesterol",
    "lipoprotein_a","apolipoprotein_a","apolipoprotein_b",
    "total_triglyceride","glucose","hba1c","igf_1",
    "crp","blood_vitamin_d","calcium","phosphate"
  ),
  stringsAsFactors = FALSE
)

biomarker_vars <- biomarker_map$variable
biomarker_vars <- biomarker_vars[biomarker_vars %in% rownames(corr_mat)]

# BP variables
bp_vars <- c("sbp_mean", "dbp_mean", "sbp_avg", "dbp_avg")
bp_vars <- bp_vars[bp_vars %in% rownames(corr_mat)]

# Helper functions
extract_cross_pairs <- function(mat, x_vars, y_vars) {
  expand.grid(
    var1 = x_vars,
    var2 = y_vars,
    stringsAsFactors = FALSE
  ) |>
    rowwise() |>
    mutate(
      correlation = as.numeric(mat[var1, var2]),
      abs_correlation = abs(correlation)
    ) |>
    ungroup() |>
    arrange(desc(abs_correlation))
}

extract_all_pairs <- function(mat) {
  vars <- rownames(mat)
  comb <- t(combn(vars, 2))
  
  data.frame(
    var1 = comb[, 1],
    var2 = comb[, 2],
    stringsAsFactors = FALSE
  ) |>
    rowwise() |>
    mutate(
      correlation = as.numeric(mat[var1, var2]),
      abs_correlation = abs(correlation)
    ) |>
    ungroup() |>
    arrange(desc(abs_correlation))
}

strong_threshold <- 0.7

# 1. Exposure x Biomarker (ALL pairs)
exposure_biomarker_pairs <- extract_cross_pairs(
  corr_mat,
  exposure_vars,
  biomarker_vars
)

write.csv(
  exposure_biomarker_pairs,
  file.path(corr_dir, "exposure_biomarker_pairs.csv"),
  row.names = FALSE
)

# 2. Exposure x SBP/DBP (ALL pairs)
exposure_bp_pairs <- extract_cross_pairs(
  corr_mat,
  exposure_vars,
  bp_vars
)

write.csv(
  exposure_bp_pairs,
  file.path(corr_dir, "exposure_bp_pairs.csv"),
  row.names = FALSE
)

# 3. All strong pairs overall only
strong_pairs_overall <- extract_all_pairs(corr_mat) |>
  filter(abs_correlation >= strong_threshold)

write.csv(
  strong_pairs_overall,
  file.path(corr_dir, "strong_pairs_overall.csv"),
  row.names = FALSE
)

# Summary
summary_lines <- c(
  paste0("Exposure variables used: ", length(exposure_vars)),
  paste0("Biomarker variables used: ", length(biomarker_vars)),
  paste0("BP variables used: ", length(bp_vars)),
  paste0("Exposure-biomarker ALL pairs: ", nrow(exposure_biomarker_pairs)),
  paste0("Exposure-BP ALL pairs: ", nrow(exposure_bp_pairs)),
  paste0("Strong threshold for overall pairs only: |r| >= ", strong_threshold),
  paste0("Overall strong pairs: ", nrow(strong_pairs_overall))
)

writeLines(
  summary_lines,
  file.path(corr_dir, "correlation_summary.txt")
)

cat("Done. Files saved in:", corr_dir, "\n")