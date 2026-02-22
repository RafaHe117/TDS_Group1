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
# 2. Exclusion - Specified ICD Congenital Diseases
# ---------------------------------------------------------------------------------------------------------

n2_before <- nrow(ukb)

# 2a. ICD Congenital Codes
exclude_codes_ICD10 <- c(
  "Q99", "Q98", "Q97", "Q96", "Q95", "Q93", "Q92", "Q91", "Q89", "Q87",
  "Q86", "Q85", "Q81", "Q80", "Q79", "Q78", "Q77", "Q76", "Q67", "Q64",
  "Q63", "Q62", "Q61", "Q60", "Q56", "Q45", "Q44", "Q43", "Q42", "Q41",
  "Q34", "Q33", "Q32", "Q28", "Q27", "Q26", "Q25", "Q24", "Q23", "Q22",
  "Q21", "Q20", "Q07", "Q06", "Q05", "Q04", "Q03", "Q02", "Q01", "Q00"
)

exclude_codes_ICD9 <- c(
  "237", "740", "741", "742", "745", "746", "747", "748", "751", "752",
  "753", "754", "755", "756", "757", "758", "759", "760"
)

# 2b. Remove participants with specified ICD10 codes
ukb <- ukb %>%
  filter(
    !if_any(
      c(starts_with("primary_diagnoses_icd10.0."),
        starts_with("secondary_diagnoses_icd10.0.")),
      ~ substr(as.character(.), 1, 3) %in% exclude_codes_ICD10
    )
  )

# 2c. Remove participants with specified ICD9 codes
ukb <- ukb %>%
  filter(
    !if_any(
      c(starts_with("primary_diagnoses_icd9.0."),
        starts_with("secondary_diagnoses_icd9.0.")),
      ~ substr(as.character(.), 1, 3) %in% exclude_codes_ICD9
    )
  )

n2_after <- nrow(ukb)

cat("\nStage 2 - Congenital ICD exclusion:\n")
cat("N before:", n2_before, "\n")
cat("N after :", n2_after, "\n")
cat("Excluded:", n2_before - n2_after, "\n")

# ---------------------------------------------------------------------------------------------------------
# 3. Exclusion - Cholesterol lowering medications
# ---------------------------------------------------------------------------------------------------------

n3_before <- nrow(ukb)

# 3a. Remove female-only questions for male participants (and vice versa)
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

# 3b. Remove all participants on cholesterol lowering medication
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

n3_after <- nrow(ukb)

cat("\nStage 3 - Cholesterol lowering medication exclusion:\n")
cat("N before:", n3_before, "\n")
cat("N after :", n3_after, "\n")
cat("Excluded:", n3_before - n3_after, "\n")

# ---------------------------------------------------------------------------------------------------------
# 4. Exclusion - Missing ethnicity / urban-rural (non-imputable structural missingness)
# ---------------------------------------------------------------------------------------------------------
n4_before <- nrow(ukb)

ukb <- ukb %>%
  filter(
    !is.na(ethnicity_5cat),
    !is.na(urban_rural_2cat)
  )

n4_after <- nrow(ukb)

cat("\nStage 4 - Missing ethnicity / urban-rural exclusion:\n")
cat("N before:", n4_before, "\n")
cat("N after :", n4_after, "\n")
cat("Excluded:", n4_before - n4_after, "\n")



# save
saveRDS(ukb, "ukb_G1_cleaned.rds")