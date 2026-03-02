# analysis/table1_config.R
# ------------------------
# Config for Table 1 pipeline (variables, labels, rules)

# Data paths (repo root)
PATH_CLEAN   <- "ukb_G1_cleaned.rds"
PATH_IMPUTED <- "ukb_G1_imputed_final.rds"

# Output directory
OUTDIR <- file.path("analysis", "output")

# Collapse rule for high-cardinality categorical variables (Appendix tables)
TOPK_LEVELS <- 10

# ------------------------
# Variable sets
# ------------------------

# Main table (paper): outcome-stratified
vars_main <- c(
  # Demographics / SES
  "age", "sex", "ethnicity_5cat",
  "imd_quintile", "urban_rural_2cat",
  "education_2", "employment_2cat", "living_with_partner",
  "work_hours_week_clean",
  # Clinical (core)
  "bmi", "sbp_mean", "dbp_mean",
  "diabetes_bin", "famhx_cvd", "birth_weight_clean",
  # Lifestyle (core)
  "smoking_3cat", "alcohol_combined"
)

# Appendix table (supplement): outcome-stratified, large set (+ medication_category)
vars_appendix <- c(
  # Lifestyle detail
  "total_met_group", "tv_group", "sleep_quality_score",
  "fruit_intake_fresh", "oily_fish_3cat", "salt_3cat",
  "coffee_intake", "tea_intake", "redwine_group",
  # Psychosocial
  "neuroticism_score", "sf_score", "mh_satis_mean_score",
  "stress_group_2yr", "mood_disorder",
  # Clinical (moved)
  "medication_category",
  # Environment
  "pm2_5_2010", "pm10_2010", "no2_2010", "noise_24h", "greenspace_pct_1000m", "temp_average",
  # Biomarkers
  "crp", "hba1c", "glucose",
  "hdl_cholesterol", "ldl_cholesterol", "cholesterol", "total_triglyceride",
  "creatinine", "cystatin_c", "urea", "urate",
  "albumin", "total_bilirubin", "direct_bilirubin",
  "alkaline_phosphatase", "alanine_aminotransferase", "aspartate_aminotransferase",
  "calcium", "phosphate", "total_blood_protein",
  "blood_vitamin_d",
  "apolipoprotein_a", "apolipoprotein_b", "lipoprotein_a",
  # Blood counts
  "wbc_count", "rbc_count", "hemoglobin",
  "platelet_count", "mpv", "mcv",
  "lymphocyte_count", "monocyte_count", "neutrophil_count",
  "eosinophil_count", "basophil_count",
  # ICD flags
  "ICD_Diabetes", "ICD_Lipidemia", "ICD_Hypertension", "ICD_AF",
  "ICD_CKD", "ICD_Mental_Health", "ICD_Migraine", "ICD_Atopy",
  "ICD_Autoimmune"
)

# Sex-stratified biological table: curated 14 variables (core)
vars_bio <- c(
  "hba1c", "glucose",
  "hdl_cholesterol", "ldl_cholesterol", "total_triglyceride",
  "crp",
  "creatinine", "cystatin_c",
  "alanine_aminotransferase", "alkaline_phosphatase",
  "albumin", "blood_vitamin_d",
  "hemoglobin", "wbc_count"
)

# Safety excludes (never in tables)
vars_exclude_always <- c(
  "eid", "yob", "date_recr",
  "dod", "age_of_death",
  "cvd_first_date", "cvd_event",
  "cvd_prevalent" # exclusion criterion, not baseline description row
)

# ------------------------
# Labels for stratification
# ------------------------
label_outcome_0 <- "No incident CVD"
label_outcome_1 <- "Incident CVD"

label_sex_f <- "Female"
label_sex_m <- "Male"