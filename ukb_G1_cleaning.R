# Libraries
library(dplyr)

# path
setwd("/rds/general/project/hda_25-26/live/TDS/TDS_Group1")

# Read data
ukb <- readRDS("ukb_G1_preprocessed.rds")

#=========================================================================================================
# Exclusion Criteria
#=========================================================================================================

# ---------------------------------------------------------------------------------------------------------
# 1. Exclusion - Prevalent CVD at recruitment
#    Remove participants with CVD before baseline to define incident cohort
# ---------------------------------------------------------------------------------------------------------

n1_before <- nrow(ukb)

ukb <- ukb %>%
  filter(cvd_prevalent == FALSE)

n1_after <- nrow(ukb)

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

cat("\nStage 3 - Missing ethnicity / urban-rural exclusion:\n")
cat("N before:", n3_before, "\n")
cat("N after :", n3_after, "\n")
cat("Excluded:", n3_before - n3_after, "\n")

# -------------------------------------------------------------------------
# Drop variables with >30% missingness and variables not used for CVD prediction
# -------------------------------------------------------------------------
drop_manual <- c(
  "shift_work",
  "mixed_shift",
  "living_alone",
  "oral_contraception",
  "household_size_clean",
  "job_shift_work_clean",
  "nrbc_count",
  "microalbumin",
  "hrt",
  "self_health_bin",
  "gp_anxdep",
  "salt_intake",
  "sodium_in_urine",
  "work_hours_unified_cat",
  "cancer_reported"
)

ukb <- ukb %>% 
  select(-any_of(drop_manual))

# save
saveRDS(ukb, "ukb_G1_cleaned.rds")