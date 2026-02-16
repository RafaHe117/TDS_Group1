# Libraries
library(dplyr)

# path
setwd("/rds/general/project/hda_25-26/live/TDS/TDS_Group1")

# Read data
ukb <- readRDS("ukb_G1_Raw.rds")

##############################################################################
# Blood pressure cleanup + mean BP
##############################################################################

# remove extreme systolic values
if ("sys_bp_automatic" %in% names(ukb)) {
  ukb$sys_bp_automatic[ukb$sys_bp_automatic > 300] <- NA
}

if ("sbp_manual" %in% names(ukb)) {
  ukb$sbp_manual[ukb$sbp_manual > 300] <- NA
}

# mean BP
if (all(c("sys_bp_automatic","sbp_manual") %in% names(ukb))) {
  ukb$sbp_mean <- rowMeans(cbind(ukb$sys_bp_automatic, ukb$sbp_manual), na.rm = TRUE)
}

if (all(c("dys_bp_automatic","dbp_manual") %in% names(ukb))) {
  ukb$dbp_mean <- rowMeans(cbind(ukb$dys_bp_automatic, ukb$dbp_manual), na.rm = TRUE)
}

if ("sys_bp_automatic" %in% names(ukb)) print(summary(ukb$sys_bp_automatic))
if ("sbp_mean" %in% names(ukb)) print(summary(ukb$sbp_mean))
if ("dbp_mean" %in% names(ukb)) print(summary(ukb$dbp_mean))

##############################################################################
# recode ethnicity -> 5 groups (White / South Asian / Chinese / Black / Other)
###############################################################################
x <- ukb$ethnicity

eth5 <- case_when(
  x %in% c("Prefer not to answer", "Do not know") ~ NA_character_,
  
  x %in% c("White", "British", "Irish", "Any other white background") ~ "White",
  
  x %in% c("Asian or Asian British",
           "Indian", "Pakistani", "Bangladeshi",
           "Any other Asian background") ~ "South Asian",
  
  x == "Chinese" ~ "Chinese",
  
  x %in% c("Black or Black British",
           "Caribbean", "African",
           "Any other Black background") ~ "Black",
  
  x %in% c("Mixed",
           "White and Black Caribbean",
           "White and Black African",
           "White and Asian",
           "Any other mixed background",
           "Other ethnic group") ~ "Other",
  
  TRUE ~ NA_character_
)

ukb$ethnicity_5cat <- factor(eth5, levels = c("White", "South Asian", "Chinese", "Black", "Other"))

print(table(ukb$ethnicity_5cat, useNA = "ifany"))
print(table(ukb$ethnicity, ukb$ethnicity_5cat, useNA = "ifany"))

##############################################################################
# recode urban_rural -> 2 groups (Urban / Rural)
##############################################################################
x <- ukb$urban_rural

ur2 <- case_when(
  grepl("Postcode not linkable", x) ~ NA_character_,
  grepl("Urban", x) ~ "Urban",
  grepl("Town and Fringe|Village|Hamlet|Rural|Small Town", x) ~ "Rural",
  TRUE ~ NA_character_
)

ukb$urban_rural_2cat <- factor(ur2, levels = c("Urban", "Rural"))

print(table(ukb$urban_rural_2cat, useNA = "ifany"))
print(table(ukb$urban_rural, ukb$urban_rural_2cat, useNA = "ifany"))

##############################################################################
# Air pollution cleanup (2010)
##############################################################################

if ("no_2010" %in% names(ukb)) {
  ukb$no_2010 <- NULL
}

if ("no2_2010" %in% names(ukb)) {
  ukb$no2_2010[ukb$no2_2010 < 0] <- NA
}

if ("pm2_5_2010" %in% names(ukb)) {
  ukb$pm2_5_2010[ukb$pm2_5_2010 < 0] <- NA
}

if ("pm10_2010" %in% names(ukb)) {
  ukb$pm10_2010[ukb$pm10_2010 < 0] <- NA
}

if ("no2_2010" %in% names(ukb)) print(summary(ukb$no2_2010))
if ("pm2_5_2010" %in% names(ukb)) print(summary(ukb$pm2_5_2010))
if ("pm10_2010" %in% names(ukb)) print(summary(ukb$pm10_2010))

##############################################################################
# Greenspace / Noise cleanup (1000m + 24h)
##############################################################################

if ("greenspace_pct_1000m" %in% names(ukb)) {
  ukb$greenspace_pct_1000m[ukb$greenspace_pct_1000m < 0] <- NA
  ukb$greenspace_pct_1000m[ukb$greenspace_pct_1000m > 100] <- NA
}

if ("noise_24h" %in% names(ukb)) {
  ukb$noise_24h[ukb$noise_24h < 0] <- NA
}

if ("greenspace_pct_1000m" %in% names(ukb)) print(summary(ukb$greenspace_pct_1000m))
if ("noise_24h" %in% names(ukb)) print(summary(ukb$noise_24h))

##############################################################################
# Biomarker cleanup
##############################################################################
if ("blood_sample_attempted" %in% names(ukb)) {
  ukb$blood_sample_attempted <- factor(ukb$blood_sample_attempted)
  print(table(ukb$blood_sample_attempted, useNA = "ifany"))
}

if ("urine_device" %in% names(ukb)) {
  print(table(ukb$urine_device, useNA = "ifany"))
  if (length(unique(na.omit(ukb$urine_device))) <= 1) ukb$urine_device <- NULL
}

biomarkers_num <- c(
  "sodium_in_urine","creatinine","microalbumin",
  "wbc_count","rbc_count","hemoglobin","mcv","platelet_count","mpv",
  "lymphocyte_count","monocyte_count","neutrophil_count",
  "eosinophil_count","basophil_count","nrbc_count",
  "albumin","alkaline_phosphatase","alanine_aminotransferase",
  "aspartate_aminotransferase","direct_bilirubin","urea","calcium",
  "cystatin_c","glucose","hba1c","igf_1","phosphate",
  "total_bilirubin","total_blood_protein","total_triglyceride",
  "urate","blood_vitamin_d","crp","lipoprotein_a",
  "ldl_cholesterol","hdl_cholesterol"
)

biomarkers_num <- biomarkers_num[biomarkers_num %in% names(ukb)]

for (v in biomarkers_num) {
  if (is.character(ukb[[v]])) ukb[[v]] <- as.numeric(ukb[[v]])
}

for (v in biomarkers_num) {
  if (is.numeric(ukb[[v]])) ukb[[v]][ukb[[v]] < 0] <- NA
}

print(sapply(biomarkers_num, function(v) mean(is.na(ukb[[v]]))))
print(lapply(biomarkers_num, function(v) summary(ukb[[v]])))

##############################################################################
# Family history (father / mother / sibling) -> first-degree CVD (binary)
##############################################################################
fam_vars <- grep("^(father|mother|sibling)_illness\\.0\\.", names(ukb), value = TRUE)

noninfo <- c(
  "Do not know (group 1)", "Do not know (group 2)",
  "Prefer not to answer (group 1)", "Prefer not to answer (group 2)"
)

cvd_terms <- c("Heart disease", "Stroke")

if (length(fam_vars) > 0) {
  fam_mat <- sapply(fam_vars, function(v) as.character(ukb[[v]]))
  fam_mat[fam_mat %in% noninfo] <- NA
  
  famhx_cvd <- apply(fam_mat, 1, function(row) {
    row <- row[!is.na(row)]
    if (length(row) == 0) return(NA_integer_)
    if (any(row %in% cvd_terms)) return(1L)
    0L
  })
  
  ukb$famhx_cvd <- factor(famhx_cvd, levels = c(0, 1), labels = c("No", "Yes"))
  print(table(ukb$famhx_cvd, useNA = "ifany"))
}

##############################################################################
# Bipolar / Depression -> Mood disorder (binary)
##############################################################################

if ("bipolar_depression_status" %in% names(ukb)) {
  
  ukb$mood_disorder <- ifelse(
    ukb$bipolar_depression_status == "No Bipolar or Depression",
    0, 1
  )
  
  ukb$mood_disorder <- factor(
    ukb$mood_disorder,
    levels = c(0, 1),
    labels = c("No", "Yes")
  )
  
  print(table(ukb$mood_disorder, useNA = "ifany"))
}

##############################################################################
# Major life events (past 2 years) -> Stress (count + 3-level group)
##############################################################################
iibs_vars <- grep("^iibs_2yr\\.0\\.", names(ukb), value = TRUE)

noninfo <- c("Prefer not to answer")
event_terms <- c(
  "Serious illness, injury or assault to yourself",
  "Serious illness, injury or assault of a close relative",
  "Death of a close relative",
  "Death of a spouse or partner",
  "Marital separation/divorce",
  "Financial difficulties"
)

if (length(iibs_vars) > 0) {
  
  iibs_mat <- sapply(iibs_vars, function(v) as.character(ukb[[v]]))
  
  stress_count <- apply(iibs_mat, 1, function(row) {
    row <- row[!row %in% noninfo]
    if (length(row) == 0) return(NA_integer_)   # all PNA
    sum(row %in% event_terms)
  })
  
  ukb$stress_count_2yr <- stress_count
  
  ukb$stress_group_2yr <- factor(
    ifelse(is.na(stress_count), NA,
           ifelse(stress_count <= 2, "0–2",
                  ifelse(stress_count <= 4, "3–4", "5–6"))),
    levels = c("0–2", "3–4", "5–6")
  )
  
  print(table(ukb$stress_count_2yr, useNA = "ifany"))
  print(table(ukb$stress_group_2yr, useNA = "ifany"))
}

##############################################################################
# GP anxiety / depression -> binary
##############################################################################

if ("gp_anxiety_depression" %in% names(ukb)) {
  
  ukb$gp_anxdep <- ifelse(
    ukb$gp_anxiety_depression == "Yes", 1,
    ifelse(ukb$gp_anxiety_depression == "No", 0, NA)
  )
  
  ukb$gp_anxdep <- factor(
    ukb$gp_anxdep,
    levels = c(0,1),
    labels = c("No","Yes")
  )
  
  print(table(ukb$gp_anxdep, useNA = "ifany"))
}

##############################################################################
# Birth weight -> clean using UKB plausible range (0.45–10 kg)
##############################################################################

if ("birth_weight" %in% names(ukb)) {
  
  ukb$birth_weight_clean <- ukb$birth_weight
  
  ukb$birth_weight_clean[
    ukb$birth_weight_clean < 0.45 | ukb$birth_weight_clean > 10
  ] <- NA
  
  cat("\nBirth weight (raw):\n")
  print(summary(ukb$birth_weight))
  
  cat("\nBirth weight (clean, 0.45–10kg):\n")
  print(summary(ukb$birth_weight_clean))
  
  cat("\nN set to NA by cleaning:",
      sum(is.na(ukb$birth_weight_clean)) - sum(is.na(ukb$birth_weight)),
      "\n")
}

##############################################################################
# Extract the highest qualification -> education
##############################################################################
highest_edu <- function(df, edu_col_base = "education") {
  
  edu_levels <- c(
    "Prefer not to answer",
    "None of the above",
    "CSEs or equivalent",
    "O levels/GCSEs or equivalent",
    "A levels/AS levels or equivalent",
    "NVQ or HND or HNC or equivalent",
    "Other professional qualifications eg: nursing, teaching",
    "College or University degree"
  )
  
  naming_pattern <- paste0("^", edu_col_base, "\\.0\\.[0-5]$")
  edu_cols <- grep(naming_pattern, names(df), value = TRUE)
  
  temp_df <- df[edu_cols]
  temp_df[] <- lapply(temp_df, factor, levels = edu_levels, ordered = TRUE)
  
  df$education <- do.call(pmax, c(temp_df, na.rm = TRUE))
  
  return(df)
}

ukb <- highest_edu(ukb)
print(table(ukb$education, useNA = "ifany"))

##############################################################################
# recode education -> education_2
##############################################################################
recode_edu <- function(df, edu_col_name = "education") {
  df$education_2 <- NA_character_
  
  df$education_2[df[[edu_col_name]] %in% c("None of the above",
                                           "Prefer not to answer",
                                           "O levels/GCSEs or equivalent",
                                           "CSEs or equivalent")] <- "10-11Yrs"
  
  df$education_2[df[[edu_col_name]] %in% c("A levels/AS levels or equivalent",
                                           "NVQ or HND or HNC or equivalent")] <- "12-14Yrs"
  
  df$education_2[df[[edu_col_name]] %in% c("College or University degree",
                                           "Other professional qualifications eg: nursing, teaching")] <- ">15Yrs"
  
  df$education_2 <- factor(df$education_2, levels = c("10-11Yrs", "12-14Yrs", ">15Yrs"))
  df
}

ukb <- recode_edu(ukb)
print(table(ukb$education_2, useNA = "ifany"))
print(table(ukb$education, ukb$education_2, useNA = "ifany"))

##############################################################################
# combine alcohol_status and alcohol_freq -> alcohol_combined
##############################################################################
recode_alcohol <- function(df, alcohol_status = "alcohol_status", alcohol_freq = "alcohol_freq") {
  df <- df %>%
    mutate(alcohol_combined = case_when(
      .data[[alcohol_status]] == "Never" & .data[[alcohol_freq]] == "Never" ~ "Never",
      .data[[alcohol_status]] == "Never" & .data[[alcohol_freq]] != "Never" ~ "Never",
      .data[[alcohol_status]] != "Never" & .data[[alcohol_freq]] == "Never" ~ "Never",
      TRUE ~ paste(.data[[alcohol_status]], .data[[alcohol_freq]], sep = ": ")
    ))
  
  df
}

ukb <- recode_alcohol(ukb)
print(table(ukb$alcohol_status, ukb$alcohol_freq, useNA = "always"))
print(table(ukb$alcohol_combined, useNA = "ifany"))

##############################################################################
# Calculate Saturated Fat (sf_score)
##############################################################################
recode_saturated_fat <- function(df, 
                                 beef_freq = "beef", 
                                 lamb_freq = "lamb", 
                                 pork_freq = "pork", 
                                 proc_freq = "processed_meat_intake", 
                                 cheese_freq = "cheese", 
                                 milk_type = "milk") {
  
  freq_map <- c(
    "Never"=0, "Less than once a week"=1, "Once a week"=2,
    "2-4 times a week"=3, "5-6 times a week"=4, "Once or more daily"=5,
    "Prefer not to answer" = NA, "Do not know" = NA
  )
  
  milk_weights <- c("Full cream"=1.9, "Semi-skimmed"=1.2, "Skimmed"=0.162, 
                    "Soya"=0.206, "Never/rarely have milk"=0, "Other type of milk"=NA,
                    "Prefer not to answer" = NA, "Do not know" = NA)
  
  meat_weights <- c("beef"=6.34, "lamb"=3.49, "pork"=6.28, "proc"=7.67, "cheese"=19.2)
  
  df <- df %>%
    mutate(
      sf_score = (meat_weights["beef"] * freq_map[as.character(df[[beef_freq]])]) +
        (meat_weights["lamb"] * freq_map[as.character(df[[lamb_freq]])]) +
        (meat_weights["pork"] * freq_map[as.character(df[[pork_freq]])]) +
        (meat_weights["proc"] * freq_map[as.character(df[[proc_freq]])]) +
        (meat_weights["cheese"] * freq_map[as.character(df[[cheese_freq]])]) +
        (milk_weights[as.character(df[[milk_type]])])
    )
  
  return(df)
}


ukb <- recode_saturated_fat(ukb)
summary(ukb$sf_score)

verification_subset <- ukb %>%
  # filter(!is.na(sf_score)) %>%
  select(beef, lamb, pork, processed_meat_intake, cheese, milk, sf_score) %>%
  head(20)

print(verification_subset)

rm(verification_subset)

##############################################################################
# Calculate mental health satisfaction score (mh_satis_mean_score)
##############################################################################
recode_mental_health <- function(df, 
                                 happiness = "happiness", work = "work_satis", 
                                 family = "family_satis", friend = "friend_satis", 
                                 finance = "fin_satis") {
  
  cat_map <- c(
    "Extremely unhappy" = 0, "I am not employed" = 0, 
    "Very unhappy" = 1, "Moderately unhappy" = 2, "Moderately happy" = 3,
    "Very happy" = 4, "Extremely happy" = 5, "Extremely happy/Other" = 6,
    "Prefer not to answer" = NA, "Do not know" = NA
  )

  cat_matrix <- cbind(
    cat_map[as.character(df[[happiness]])], 
    cat_map[as.character(df[[work]])], 
    cat_map[as.character(df[[family]])],
    cat_map[as.character(df[[friend]])], 
    cat_map[as.character(df[[finance]])]
  )
  
  # Calculate how many questions each person answered
  n_answered <- rowSums(!is.na(cat_matrix))
  
  df <- df %>%
    mutate(
      mh_n_answered = n_answered,
      mh_satis_mean_score = ifelse(n_answered >= 4, rowMeans(cat_matrix, na.rm = TRUE), NA)
    )
  
  return(df)
}

ukb <- recode_mental_health(ukb)

table(ukb$mh_n_answered, is.na(ukb$mh_satis_mean_score), useNA = "always")

##############################################################################
# recode had menopause (Field 2724) -> 2 groups (Yes / No) ；Indeterminate categories -> NA
# Restrict to females; males set to NA
# Menopause is strongly age-dependent.Age adjustment or stratification required in downstream models.
##############################################################################

recode_menopause <- function(df,
                             var = "menopause",
                             sex_var = "sex") {
  
  df$menopause_bin <- NA_character_
  
  # Define only among females
  df$menopause_bin[df[[sex_var]] == "Female" & df[[var]] == "Yes"] <- "Yes"
  df$menopause_bin[df[[sex_var]] == "Female" & df[[var]] == "No"]  <- "No"
  
  df$menopause_bin <- factor(df$menopause_bin,
                             levels = c("No", "Yes"))
  
  df
}

ukb <- recode_menopause(ukb)

# quick checks
print(table(ukb$sex, ukb$menopause_bin, useNA = "ifany"))
print(table(ukb$menopause, ukb$menopause_bin, useNA = "ifany"))

##############################################################################
# recode self-rated health (Field 2178) -> 2 groups (Good vs Poor)
# "Prefer not to answer" and "Do not know" -> missing (NA)
# Self-rated health may act as a mediator or collider
##############################################################################

recode_self_health <- function(df, var = "self_health_rating") {
  
  df$self_health_bin <- NA_character_
  
  df$self_health_bin[df[[var]] %in% c("Excellent", "Good")] <- "Good"
  df$self_health_bin[df[[var]] %in% c("Fair", "Poor")]      <- "Poor"
  
  df$self_health_bin <- factor(df$self_health_bin,
                               levels = c("Good", "Poor"))
  
  df
}

ukb <- recode_self_health(ukb)

# quick checks
print(table(ukb$self_health_bin, useNA = "ifany"))
print(table(ukb$self_health_rating, ukb$self_health_bin, useNA = "ifany"))


# NOTE:
# Household structure is conceptually distinct from SES.
# Household income will be used as primary socioeconomic indicator.

##############################################################################
# recode living with partner (Field 6141)
# Define Yes if "Husband, wife or partner" ；Prefer not to answer -> NA
##############################################################################

recode_partner <- function(df, var = "household_relationship") {
  
  df$partner_status <- NA_character_
  
  # Yes
  df$partner_status[df[[var]] == "Husband, wife or partner"] <- "Yes"
  
  # No (but exclude explicit non-response)
  df$partner_status[df[[var]] != "Husband, wife or partner" &
                      df[[var]] != "Prefer not to answer" &
                      !is.na(df[[var]])] <- "No"
  
  df$partner_status <- factor(df$partner_status,
                              levels = c("No", "Yes"))
  
  df
}

ukb <- recode_partner(ukb)

# quick checks
print(table(ukb$partner_status, useNA="ifany"))
print(table(ukb$household_relationship,
            ukb$partner_status,
            useNA="ifany"))

##############################################################################
# Clean household_size (Field 709) - Create living_alone (1 vs >=2)
# Convert to numeric ；Set special codes to NA (-9999, 999, etc.) ；Remove implausible values (>20)
##############################################################################

clean_household_size <- function(df, var = "household_size") {
  
  # 1️⃣ Convert to character first (in case factor)
  x <- as.character(df[[var]])
  
  # 2️⃣ Set explicit missing strings to NA
  x[x %in% c("Do not know",
             "Prefer not to answer",
             "",
             "NA")] <- NA
  
  # 3️⃣ Convert to numeric
  x_num <- suppressWarnings(as.numeric(x))
  
  # 4️⃣ Set known UKB missing codes to NA
  x_num[x_num %in% c(-9999, -999, -1)] <- NA
  
  # 5️⃣ Remove implausible household sizes
  # UKB typical max household size < 15
  x_num[x_num > 20] <- NA
  
  # 6️⃣ Save cleaned variable
  df$household_size_clean <- x_num
  
  # 7️⃣ Create living_alone
  df$living_alone <- NA_character_
  
  df$living_alone[x_num == 1] <- "Yes"
  df$living_alone[x_num >= 2] <- "No"
  
  df$living_alone <- factor(df$living_alone,
                            levels = c("No", "Yes"))
  
  df
}

# Apply
ukb <- clean_household_size(ukb)

# QUICK CHECKS

# Check distribution of cleaned size
summary(ukb$household_size_clean)

# Check frequency table
table(ukb$household_size_clean, useNA = "ifany")

# Check living alone
table(ukb$living_alone, useNA = "ifany")

# Cross-check
table(ukb$household_size_clean, ukb$living_alone, useNA="ifany")


############################################################
# recode employment_status (Field 6142) -> 2 groups (In paid，Not in paid)
# Exclude: Unable to work due to sickness/disability
############################################################
x <- ukb$employment_status

emp2 <- case_when(
  
  # 1. Non-response
  x %in% c("None of the above", "Prefer not to answer") ~ NA_character_,
  
  # 2. Exclusion group (tutor comment)
  grepl("Unable to work", x, ignore.case = TRUE) ~ NA_character_,
  
  # 3. Priority group: In paid employment
  grepl("In paid employment|self-employed", x, ignore.case = TRUE) 
  ~ "In paid employment",
  
  # 4. Not in paid employment
  grepl("Retired", x, ignore.case = TRUE) ~ "Not in paid employment",
  grepl("Unemployed", x, ignore.case = TRUE) ~ "Not in paid employment",
  grepl("Looking after home", x, ignore.case = TRUE) ~ "Not in paid employment",
  grepl("Doing unpaid", x, ignore.case = TRUE) ~ "Not in paid employment",
  grepl("student", x, ignore.case = TRUE) ~ "Not in paid employment",
  
  TRUE ~ NA_character_
)

ukb$employment_2cat <- factor(
  emp2,
  levels = c("Not in paid employment", "In paid employment")
)

# Quick checks
print(table(ukb$employment_2cat, useNA = "ifany"))

print(table(
  ukb$employment_status,
  ukb$employment_2cat,
  useNA = "ifany"
))








# save
saveRDS(ukb, "ukb_G1_preprocessed.rds")


