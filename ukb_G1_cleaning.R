library(dplyr)

setwd("/rds/general/project/hda_25-26/live/TDS/TDS_Group1")

ukb <- readRDS("ukb_G1_preprocessed.rds")
dir.create("cleaning.output", showWarnings = FALSE, recursive = TRUE)

# =========================================================================================
# Logging objects
# =========================================================================================
stage_summary <- data.frame(
  stage = character(),
  n_before = integer(),
  n_after = integer(),
  n_excluded = integer(),
  stringsAsFactors = FALSE
)

append_stage <- function(stage_name, n_before, n_after) {
  data.frame(
    stage = stage_name,
    n_before = n_before,
    n_after = n_after,
    n_excluded = n_before - n_after,
    stringsAsFactors = FALSE
  )
}

# =========================================================================================
# Exclusion Criteria
# =========================================================================================

# ---------------------------------------------------------------------------------------------------------
# 1. Exclusion - Prevalent CVD at recruitment
#    Remove participants with CVD before baseline to define incident cohort
# ---------------------------------------------------------------------------------------------------------

n1_before <- nrow(ukb)

ukb <- ukb %>%
  filter(cvd_prevalent == FALSE)

n1_after <- nrow(ukb)

stage_summary <- bind_rows(
  stage_summary,
  append_stage("Stage 1 - Prevalent CVD exclusion", n1_before, n1_after)
)

cat("\nStage 1 - Prevalent CVD exclusion:\n")
cat("N before:", n1_before, "\n")
cat("N after :", n1_after, "\n")
cat("Excluded:", n1_before - n1_after, "\n")

# ---------------------------------------------------------------------------------------------------------
# 2. Exclusion - Cholesterol lowering medications
# ---------------------------------------------------------------------------------------------------------

n2_before <- nrow(ukb)

# 2a. Remove female-only questions for male participants (and vice versa)
ukb <- ukb %>%
  mutate(
    # female-only medication items -> set NA for males
    across(
      c(
        chol_bp_diabetes_hormone_medication.0.0,
        chol_bp_diabetes_hormone_medication.0.1,
        chol_bp_diabetes_hormone_medication.0.2,
        chol_bp_diabetes_hormone_medication.0.3
      ),
      ~ ifelse(as.character(sex) == "Male", NA, as.character(.))
    ),
    # male-only medication items -> set NA for females
    across(
      c(
        chol_bp_diabetes_medication.0.0,
        chol_bp_diabetes_medication.0.1,
        chol_bp_diabetes_medication.0.2
      ),
      ~ ifelse(as.character(sex) == "Female", NA, as.character(.))
    )
  )

# 2b. Remove all participants on cholesterol lowering medication
ukb <- ukb %>%
  filter(
    !if_any(
      c(
        chol_bp_diabetes_hormone_medication.0.0,
        chol_bp_diabetes_hormone_medication.0.1,
        chol_bp_diabetes_hormone_medication.0.2,
        chol_bp_diabetes_hormone_medication.0.3,
        chol_bp_diabetes_medication.0.0,
        chol_bp_diabetes_medication.0.1,
        chol_bp_diabetes_medication.0.2
      ),
      ~ as.character(.) == "Cholesterol lowering medication" & !is.na(.)
    )
  )

n2_after <- nrow(ukb)

stage_summary <- bind_rows(
  stage_summary,
  append_stage("Stage 2 - Cholesterol lowering medication exclusion", n2_before, n2_after)
)

cat("\nStage 2 - Cholesterol lowering medication exclusion:\n")
cat("N before:", n2_before, "\n")
cat("N after :", n2_after, "\n")
cat("Excluded:", n2_before - n2_after, "\n")

# Drop cholesterol medication variables (used only for exclusion)
ukb <- ukb %>%
  select(
    -chol_bp_diabetes_hormone_medication.0.0,
    -chol_bp_diabetes_hormone_medication.0.1,
    -chol_bp_diabetes_hormone_medication.0.2,
    -chol_bp_diabetes_hormone_medication.0.3,
    -chol_bp_diabetes_medication.0.0,
    -chol_bp_diabetes_medication.0.1,
    -chol_bp_diabetes_medication.0.2
  )

# ---------------------------------------------------------------------------------------------------------
# 3. Exclusion - Missing ethnicity / urban-rural (non-imputable structural missingness)
# ---------------------------------------------------------------------------------------------------------

n3_before <- nrow(ukb)

ukb <- ukb %>%
  filter(
    !is.na(ethnicity_5cat),
    !is.na(urban_rural_2cat)
  )

n3_after <- nrow(ukb)

stage_summary <- bind_rows(
  stage_summary,
  append_stage("Stage 3 - Missing ethnicity / urban-rural exclusion", n3_before, n3_after)
)

cat("\nStage 3 - Missing ethnicity / urban-rural exclusion:\n")
cat("N before:", n3_before, "\n")
cat("N after :", n3_after, "\n")
cat("Excluded:", n3_before - n3_after, "\n")

# ---------------------------------------------------------------------------------------------------------
# 4. Automatically drop variables with ≥30% missingness
#    Retain key outcome/time variables (e.g., dod, cvd_first_date)
# ---------------------------------------------------------------------------------------------------------

n4_before_var <- ncol(ukb)

thr <- 0.30

na_prop   <- colMeans(is.na(ukb))   # proportion missing per variable
high_miss <- names(na_prop[na_prop >= thr])

protected <- intersect(c("dod", "cvd_first_date"), names(ukb))
drop_auto <- setdiff(high_miss, protected)

# Table of dropped variables and missingness
missingness_dropped <- data.frame(
  variable = drop_auto,
  missing_n = colSums(is.na(ukb[, drop_auto, drop = FALSE])),
  missing_prop = as.numeric(na_prop[drop_auto]),
  missing_percent = round(as.numeric(na_prop[drop_auto]) * 100, 2),
  stringsAsFactors = FALSE
) %>%
  arrange(desc(missing_prop))

cat("\nDropping variables with missingness >=", thr * 100, "% (protected kept):\n")
print(missingness_dropped)

ukb <- ukb %>% 
  select(-any_of(drop_auto))

n4_after_var <- ncol(ukb)

stage_summary <- bind_rows(
  stage_summary,
  data.frame(
    stage = paste0("Stage 4 - Drop variables with missingness >= ", thr * 100, "%"),
    n_before = n4_before_var,
    n_after = n4_after_var,
    n_excluded = n4_before_var - n4_after_var,
    stringsAsFactors = FALSE
  )
)

cat("\nStage 4 - Missingness-based variable exclusion:\n")
cat("Variables before:", n4_before_var, "\n")
cat("Variables after :", n4_after_var, "\n")
cat("Variables dropped:", n4_before_var - n4_after_var, "\n")

# =========================================================================================
# Save outputs
# =========================================================================================

# save cleaned dataset
saveRDS(ukb, "ukb_G1_cleaned.rds")

# save stage summary
write.csv(
  stage_summary,
  file = "cleaning.output/cleaning_stage_summary.csv",
  row.names = FALSE
)

# save missingness-dropped variables summary
write.csv(
  missingness_dropped,
  file = "cleaning.output/missingness_dropped_variables.csv",
  row.names = FALSE
)

cat("\nSaved cleaned dataset to: ukb_G1_cleaned.rds\n")
cat("Saved stage summary to: cleaning.output/cleaning_stage_summary.csv\n")
cat("Saved dropped-variable summary to: cleaning.output/missingness_dropped_variables.csv\n")