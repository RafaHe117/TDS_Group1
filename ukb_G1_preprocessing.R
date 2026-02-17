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
# Diabetes diagnosed -> binary
##############################################################################
if ("diabetes_diagnosed" %in% names(ukb)) {
  
  ukb$diabetes_bin <- ifelse(
    ukb$diabetes_diagnosed == "Yes", 1,
    ifelse(ukb$diabetes_diagnosed == "No", 0, NA)
  )
  
  ukb$diabetes_bin <- factor(
    ukb$diabetes_bin,
    levels = c(0, 1),
    labels = c("No", "Yes")
  )
  
  print(table(ukb$diabetes_bin, useNA = "ifany"))
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
# TV duration -> numeric + group
##############################################################################

if ("tv_duration" %in% names(ukb)) {
  
  noninfo <- c("Do not know", "Prefer not to answer")
  
  tv_raw <- as.character(ukb$tv_duration)
  tv_raw[tv_raw %in% noninfo] <- NA
  tv_num <- ifelse(tv_raw == "Less than an hour a day", 0.5, tv_raw)
  tv_num <- suppressWarnings(as.numeric(tv_num))
  tv_num[tv_num < 0 | tv_num > 24] <- NA
  
  ukb$tv_hours_day <- tv_num
  
  ukb$tv_group <- factor(
    ifelse(is.na(tv_num), NA,
           ifelse(tv_num < 1, "<1h",
                  ifelse(tv_num <= 3, "1–3h",
                         ifelse(tv_num <= 5, "3–5h", ">5h")))),
    levels = c("<1h", "1–3h", "3–5h", ">5h")
  )
  
  print(table(ukb$tv_group, useNA = "ifany"))
}

##############################################################################
# Total MET minutes weekly -> cleaned + group
##############################################################################

if ("total_met_minutes_weekly" %in% names(ukb)) {
  
  met_clean <- ukb$total_met_minutes_weekly
  met_clean[met_clean < 0] <- NA_real_
  
  ukb$total_met_min_wk <- met_clean
  
  q <- quantile(met_clean, probs = c(1/3, 2/3), na.rm = TRUE, type = 2)
  
  ukb$total_met_group <- factor(
    ifelse(is.na(met_clean), NA,
           ifelse(met_clean <= q[1], "Low",
                  ifelse(met_clean <= q[2], "Moderate", "High"))),
    levels = c("Low", "Moderate", "High")
  )
  
  print(table(ukb$total_met_group, useNA = "ifany"))
}

##############################################################################
# Sleep quality score (duration + insomnia + snoring)
##############################################################################
tmp <- ukb$sleep_duration
tmp[tmp %in% c("Do not know", "Prefer not to answer")] <- NA
sleep_num <- as.numeric(tmp)

sleep_num[sleep_num < 1 | sleep_num > 23] <- NA

sleep_duration_pt <- ifelse(
  is.na(sleep_num), NA_integer_,
  ifelse(sleep_num %in% c(7, 8), 1L, 0L)
)

insomnia_pt <- ifelse(
  ukb$insomnia == "Never/rarely", 2L,
  ifelse(ukb$insomnia == "Sometimes", 1L,
         ifelse(ukb$insomnia == "Usually", 0L, NA_integer_))
)

snoring_pt <- ifelse(
  ukb$snoring == "No", 1L,
  ifelse(ukb$snoring == "Yes", 0L, NA_integer_)
)

pt_mat <- cbind(sleep_duration_pt, insomnia_pt, snoring_pt)
n_obs  <- rowSums(!is.na(pt_mat))

ukb$sleep_quality_score <- ifelse(
  n_obs < 2, NA_integer_,
  rowSums(pt_mat, na.rm = TRUE)
)

print(summary(ukb$sleep_quality_score))
print(table(ukb$sleep_quality_score, useNA = "ifany"))

##############################################################################
# Smoking -> combined 3-level variable & drop others
##############################################################################

if ("smoking_status" %in% names(ukb)) {
  
  smk <- ukb$smoking_status
  smk[smk == "Prefer not to answer"] <- NA
  
  ukb$smoking_3cat <- factor(
    smk,
    levels = c("Never", "Previous", "Current")
  )
  
  print(table(ukb$smoking_3cat, useNA = "ifany"))
  
  ukb <- ukb[, !names(ukb) %in% 
               c("pack_years_smoked", "cigarettes_per_day_prev")]
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
  
  return(df)
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
  
  return(df)
}

ukb <- recode_alcohol(ukb)
print(table(ukb$alcohol_status, ukb$alcohol_freq, useNA = "always"))
print(table(ukb$alcohol_combined, useNA = "ifany"))

##############################################################################
# Bin average weekly red wine consumption
##############################################################################
bin_redwine <- function(df, redwine_consumption = "redwine") {
  
  df <- df %>%
    mutate(
      redwine = as.numeric(.data[[redwine_consumption]]),
      
      redwine_group = case_when(
        is.na(redwine) ~ NA_character_,
        redwine == 0 ~ "Non-drinker",
        redwine >= 1 & redwine <= 14 ~ "Light drinker",
        redwine >= 15 & redwine <= 21 ~ "Moderate drinker",
        redwine >= 22 ~ "Heavy drinker"
      )
    )
  
  df$redwine_group <- factor(df$redwine_group, levels = c("Non-drinker", "Light drinker", 
                                                          "Moderate drinker", "Heavy drinker"), ordered = TRUE)
  return(df)
}

ukb <- bin_redwine(ukb)

print(table(ukb$redwine_group, useNA = "always"))
tapply(ukb$redwine, ukb$redwine_group, summary)

##############################################################################
# Calculate Saturated Fat (sf_score)
##############################################################################
recode_saturated_fat <- function(df, 
                                 beef_freq = "beef", lamb_freq = "lamb", 
                                 pork_freq = "pork", proc_freq = "processed_meat_intake", 
                                 cheese_freq = "cheese", milk_type = "milk") {
  
  freq_map <- c(
    "Never"=0, "Less than once a week"=1, "Once a week"=2,
    "2-4 times a week"=3, "5-6 times a week"=4, "Once or more daily"=5,
    "Prefer not to answer" = NA, "Do not know" = NA
  )
  
  milk_weights <- c("Full cream"=1.9, "Semi-skimmed"=1.2, "Skimmed"=0.162, 
                    "Soya"=0.206, "Never/rarely have milk"=0, "Other type of milk"=NA,
                    "Prefer not to answer" = NA, "Do not know" = NA)
  
  meat_weights <- c("beef"=6.34, "lamb"=3.49, "pork"=6.28, "proc"=7.67, "cheese"=19.2)
  
  max_vals <- c(
    beef = unname(meat_weights["beef"]) * 5, 
    lamb = unname(meat_weights["lamb"]) * 5,
    pork = unname(meat_weights["pork"]) * 5, 
    proc = unname(meat_weights["proc"]) * 5,
    cheese = unname(meat_weights["cheese"]) * 5, 
    milk = 1.9
  )
  
  # Total maximum possible score if all 6 items were maxed out
  total_max <- sum(max_vals)
  
  df <- df %>%
    # Calculate individual component scores
    mutate(
      beef_val = unname(meat_weights["beef"]) * freq_map[as.character(.data[[beef_freq]])],
      lamb_val = unname(meat_weights["lamb"]) * freq_map[as.character(.data[[lamb_freq]])],
      pork_val = unname(meat_weights["pork"]) * freq_map[as.character(.data[[pork_freq]])],
      proc_val = unname(meat_weights["proc"]) * freq_map[as.character(.data[[proc_freq]])],
      cheese_val = unname(meat_weights["cheese"]) * freq_map[as.character(.data[[cheese_freq]])],
      milk_val = milk_weights[as.character(.data[[milk_type]])]
    ) %>%
    
    mutate(
      n_valid = rowSums(!is.na(cbind(beef_val, lamb_val, pork_val, proc_val, cheese_val, milk_val))),
      
      # Sum the valid components
      sum_achieved = rowSums(cbind(beef_val, lamb_val, pork_val, proc_val, cheese_val, milk_val), na.rm = TRUE),
      
      # Sum the maximum theoretical scores only for the questions the participant answered
      sum_max_possible = 
        (max_vals["beef"] * !is.na(beef_val)) +
        (max_vals["lamb"] * !is.na(lamb_val)) +
        (max_vals["pork"] * !is.na(pork_val)) +
        (max_vals["proc"] * !is.na(proc_val)) +
        (max_vals["cheese"] * !is.na(cheese_val)) +
        (max_vals["milk"] * !is.na(milk_val)),
    
      sf_score = if_else(n_valid >= 3, (sum_achieved / sum_max_possible) * total_max, NA_real_)
    ) %>%
    
    # Remove the calculation columns
    select(-ends_with("_val"), -n_valid, -sum_achieved, -sum_max_possible)
  
  return(df)
}


ukb <- recode_saturated_fat(ukb)
summary(ukb$sf_score)

verification_subset <- ukb %>%
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
                                 finance = "fin_satis", employment_col_base = "employment_status") {
  
  cat_map <- c(
    "Extremely unhappy" = 0, "I am not employed" = NA, 
    "Very unhappy" = 1, "Moderately unhappy" = 2, "Moderately happy" = 3,
    "Very happy" = 4, "Extremely happy" = 5, "Extremely happy/Other" = 6,
    "Prefer not to answer" = NA, "Do not know" = NA
  )

  cat_matrix <- cbind(
    happiness = cat_map[as.character(df[[happiness]])], 
    work      = cat_map[as.character(df[[work]])], 
    family    = cat_map[as.character(df[[family]])],
    friend    = cat_map[as.character(df[[friend]])], 
    finance   = cat_map[as.character(df[[finance]])]
  )
  
  # Extract the employment status columns
  employment_col <- grep(paste0("^", employment_col_base), names(df), value = TRUE)
  employment_matrix <- as.matrix(df[, employment_col])
  
  # Check unemployed or retired
  is_retired <- rowSums(employment_matrix == "Retired", na.rm = TRUE) > 0
  is_unemployed <- rowSums(employment_matrix == "Unemployed", na.rm = TRUE) > 0
  
  # Overwrite the work score for unemployed and for retired
  cat_matrix[is_retired, "work"] <- 3
  cat_matrix[is_unemployed, "work"] <- 0
  
  # Calculate how many questions each person answered
  n_answered <- rowSums(!is.na(cat_matrix))
  
  df <- df %>%
    mutate(
      mh_n_answered = n_answered,
      mh_satis_mean_score = ifelse(n_answered >= 3, rowMeans(cat_matrix, na.rm = TRUE), NA)
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
  levels = c("Not in paid employment", "In paid employment"))


# Quick checks
print(table(ukb$employment_2cat, useNA = "ifany"))

print(table(
  ukb$employment_status,
  ukb$employment_2cat,
  useNA = "ifany"
))


# =================================================================
# Cleaning night shift variable
# INPUTS:
#   - ukb : job_shift_work, employment_status
#
# OUTPUT:
#  - job_shift_work_clean : factor version, non-informative values = NA, 
#                           including only participants who are employed
#   - ever_nights : (OPTIONAL, commented out) simplified binary version of 
#       night shift work
# =================================================================

# 1) Cleaning by setting non-informative values to NA

ukb <- ukb %>%
  mutate(
    job_shift_work_clean = as.character(job_shift_work),
    job_shift_work_clean = ifelse(
      job_shift_work_clean %in% c("Prefer not to answer", "Do not know"),
      NA_character_,
      job_shift_work_clean
    ),
    
# 2) Only employed participants should have exposure values, other set to NA
    job_shift_work_clean = ifelse(
      employment_status == "In paid employment or self-employed",
      job_shift_work_clean,
      NA_character_
    ),
    
# 3) Convert final result to factor
    job_shift_work_clean = factor(job_shift_work_clean)
  )


# ------------------------------------------------------------
# OPTION: Simplify to a binary indicator ever_nights
# ------------------------------------------------------------
# ukb <- ukb %>%
#   mutate(
#     ever_nights = case_when(
#       is.na(job_shift_work_clean) ~ NA,
#       job_shift_work_clean == "Never/rarely" ~ 0,
#       TRUE ~ 1
#     )
#   )



# ============================================================
# Weekly work hours cleaning + unified categorical variable
# INPUTS:
#   - work_hours_week : exact hours/week (numeric)
#   - work_hours_cat  : categorical hours/week
#
# Outputs added to ukb:
#   - work_hours_week_clean  : exact hours with values < 0 recoded to NA
#   - work_hours_unified_cat : FINAL variable used for analysis/imputation
#       ("15–<20", "20–<30", "30–40", ">40") as an ordered factor, including only currently employed participants
# ============================================================

# 1) Cleaning exact hours: UKB negative codes -> NA
ukb <- ukb %>%
  mutate(
    work_hours_week_clean = if_else(work_hours_week < 0, NA_real_, work_hours_week)
  )

# 2) Recode lumped categories to match final bins
ukb <- ukb %>%
  mutate(
    work_hours_cat_clean = case_when(
      as.character(work_hours_cat) == "15 to less-than-20 hours" ~ "15–<20",
      as.character(work_hours_cat) == "20 to less-than-30 hours" ~ "20–<30",
      as.character(work_hours_cat) == "30 to 40 hours"           ~ "30–40",
      as.character(work_hours_cat) == "Over 40 hours"            ~ ">40",
      TRUE ~ NA_character_
    )
  )

# 3) Bin exact hours into same categories
ukb <- ukb %>%
  mutate(
    work_hours_exact_binned = case_when(
      !is.na(work_hours_week_clean) & work_hours_week_clean < 20  ~ "15–<20",
      !is.na(work_hours_week_clean) & work_hours_week_clean < 30  ~ "20–<30",
      !is.na(work_hours_week_clean) & work_hours_week_clean <= 40 ~ "30–40",
      !is.na(work_hours_week_clean) & work_hours_week_clean > 40  ~ ">40",
      TRUE ~ NA_character_
    )
  )

# 4) FINAL unified variable: prioritise exact (binned), else lumped;
#        ONLY for people in paid employment/self-employment
ukb <- ukb %>%
  mutate(
    work_hours_unified_cat = case_when(
      employment_status == "In paid employment or self-employed" &
        !is.na(work_hours_exact_binned) ~ work_hours_exact_binned,

      employment_status == "In paid employment or self-employed" &
        is.na(work_hours_exact_binned) &
        !is.na(work_hours_cat_clean) ~ work_hours_cat_clean,

      TRUE ~ NA_character_
    ),

    work_hours_unified_cat = factor(
      work_hours_unified_cat,
      levels = c("15–<20", "20–<30", "30–40", ">40"),
      ordered = TRUE
    )
  )


#  5) DROP intermediate columns to avoid confusion ----
ukb <- ukb %>% select(-work_hours_cat_clean, -work_hours_exact_binned)



# ============================================================
# IMD cleaning and Assignment via Decision Tree
# INPUTS
#   - Columns (raw IMD): imd_england, imd_scotland, imd_wales
#
# OUTPUTS (new columns added to ukb)
#   - imd_final   : single IMD score per participant using precedence rule
#                  England -> Scotland -> Wales -> NA if all missing
#   - imd_source  : which nation IMD was used to create imd_final
#   - imd_quintile: quintile based on imd_final (computed on non-missing)
#
#==========================================================

# 1) Replace invalid values: set any IMD value < 0 to NA

ukb <- ukb %>%
  mutate(
    across(
      c(imd_england, imd_scotland, imd_wales),
      ~ ifelse(.x < 0, NA_real_, .x)
    )
  )


# 2) Decision-tree assignment England -> Scotland -> Wales -> NA

ukb <- ukb %>%
  mutate(
    imd_final = coalesce(imd_england, imd_scotland, imd_wales),
    # Track which IMD was used
    imd_source = case_when(
      !is.na(imd_england)  ~ "England",
      is.na(imd_england) & !is.na(imd_scotland) ~ "Scotland",
      is.na(imd_england) & is.na(imd_scotland) & !is.na(imd_wales) ~ "Wales",
      TRUE ~ NA_character_
    )
  )

# QC table: counts + percentages of which source was used (including missing)
imd_source_qc <- ukb %>%
  count(imd_source) %>%
  mutate(
    percent = percent(n / sum(n), accuracy = 0.1)
  )

print(imd_source_qc)


# 3) Create quintiles for comparability

ukb <- ukb %>%
  mutate(
    imd_quintile = ifelse(
      is.na(imd_final),
      NA_integer_,
      ntile(imd_final, 5)
    ),
    # Transform to ordered factor
    imd_quintile = factor(imd_quintile, levels = 1:5, ordered = TRUE)
  )

#==========================================================================================
# Pre-processing (Medication categories, oral contraception, HRT, ICD Co-morbidities)
#==========================================================================================
# 1. Number of Medication Categories (0, 1-2, >=3)
ukb <- ukb %>%
  mutate(
    medication_category = case_when(
      total_medications == 0     ~ "0",
      total_medications %in% 1:2 ~ "1-2",
      total_medications >= 3     ~ ">=3",
      TRUE                       ~ NA_character_
    ),
    medication_category = factor(
      medication_category, 
      levels = c("0", "1-2", ">=3")
    )
  )


# 2. Oral Contraception
## Set male results all to "No", then "Prefer not to answer" and "Do not know" to NA
ukb <- ukb %>%
mutate(
    oral_contraception = as.character(oral_contraception),
    oral_contraception = case_when(
      sex == "Male" ~ "No",
      oral_contraception == "Prefer not to answer" ~ NA_character_,
      oral_contraception == "Do not know"           ~ NA_character_,
      TRUE                                          ~ oral_contraception
    ),
        oral_contraception = as.factor(oral_contraception)
  )

table(ukb$sex, ukb$oral_contraception, useNA = "always")


# 3. HRT
## Set male results all to "No", then "Prefer not to answer" and "Do not know" to NA
ukb <- ukb %>%
mutate(
    hrt = as.character(hrt),
    hrt = case_when(
      sex == "Male" ~ "No",
      hrt == "Prefer not to answer" ~ NA_character_,
      hrt == "Do not know"           ~ NA_character_,
      TRUE                                          ~ hrt
    ),
        hrt = as.factor(hrt)
  )

table(ukb$sex, ukb$hrt, useNA = "always")



# 4. Coding Presence of ICD Comorbidities
## 4a. ICD Comorbidity Codes

Diabetes_ICD_10 <- c("E10", "E11","E12", "E13", "E14")
Diabetes_ICD_9 <- c("250")

Lipidemia_ICD_10 <- c("E78")
Lipidemia_ICD_9 <- c("272")

Mental_Health_ICD_10 <- c("F20", "F31", "F32")
Mental_Health_ICD_9 <- c("231", "296", "298", "311")

Migraine_ICD_10 <- c("G43")
Migraine_ICD_9 <- c("346")

Hypertension_ICD_10 <- c("I10")
Hypertension_ICD_9 <- c("401")

AF_ICD_10 <- c("I48")
AF_ICD_9<- c("42731", "42732")

Atopy_ICD_10 <- c("J30", "J45", "L20")
Atopy_ICD_9 <- c("477", "493", "691")

Autoimmune_ICD_10 <- c("L23", "L40", "L93", "L95", "M05", "M06")
Autoimmune_ICD_9 <- c("692", "695", "696", "709", "714")

CKD_ICD_10 <- c("N18")
CKD_ICD_9<- c("585")

## 4b. ICD Comorbidity Code Checking Function

check_icd <- function(df, icd10_list, icd9_list) {
  res <- df %>%
    transmute(
      flag = if_any(
        c(starts_with("primary_diagnoses_icd10.0."), starts_with("secondary_diagnoses_icd10.0.")),
        ~ substr(as.character(.), 1, 3) %in% icd10_list
      ) | 
      if_any(
        c(starts_with("primary_diagnoses_icd9.0."), starts_with("secondary_diagnoses_icd9.0.")),
        ~ substr(as.character(.), 1, 3) %in% icd9_list
      )
    )
  return(as.integer(res$flag))
}


## 4c. Apply ICD Comorbidity Code Checking Function
ukb <- ukb %>%
  mutate(
    ICD_Diabetes      = check_icd(., Diabetes_ICD_10, Diabetes_ICD_9),
    ICD_Lipidemia      = check_icd(., Lipidemia_ICD_10, Lipidemia_ICD_9),
    ICD_Mental_Health  = check_icd(., Mental_Health_ICD_10, Mental_Health_ICD_9),
    ICD_Migraine       = check_icd(., Migraine_ICD_10, Migraine_ICD_9),
    ICD_Hypertension   = check_icd(., Hypertension_ICD_10, Hypertension_ICD_9),
    ICD_AF             = check_icd(., AF_ICD_10, AF_ICD_9),
    ICD_Atopy          = check_icd(., Atopy_ICD_10, Atopy_ICD_9),
    ICD_Autoimmune     = check_icd(., Autoimmune_ICD_10, Autoimmune_ICD_9),
    ICD_CKD            = check_icd(., CKD_ICD_10, CKD_ICD_9)
  )

#=========================================================================================================
# Exclusion (ICD Congenital Disease, Cholesterol Lowering Medication, Cancer within 2 years of recruitment)
#==========================================================================================================
# 1. Removing Specified ICD Congenital Diseases

## 1a. ICD Congenital Codes
exclude_codes_ICD10 <- c(
  "Q99", "Q98", "Q97", "Q96", "Q95", "Q93", "Q92", "Q91", "Q89", "Q87", 
  "Q86", "Q85", "Q81", "Q80", "Q79", "Q78", "Q77", "Q76", "Q67", "Q64", 
  "Q63", "Q62", "Q61", "Q60", "Q56", "Q45", "Q44", "Q43", "Q42", "Q41", 
  "Q34", "Q33", "Q32", "Q28", "Q27", "Q26", "Q25", "Q24", "Q23", "Q22", 
  "Q21", "Q20", "Q07", "Q06", "Q05", "Q04", "Q03", "Q02", "Q01", "Q00"
)


exclude_codes_ICD9 <- c("237", "740", "741", "742", "745", "746", "747", "748", "751", "752", "753", "754", "755", "756", "757", "758", "759", "760")



## 1b. Remove participants with specified ICD10 or ICD9 codes

ukb <- ukb %>%
filter(
  !if_any(
    c(starts_with("primary_diagnoses_icd10.0."), 
      starts_with("secondary_diagnoses_icd10.0.")),
    ~ substr(., 1, 3) %in% exclude_codes_ICD10
  )
)

ukb <- ukb %>%
filter(
  !if_any(
    c(starts_with("primary_diagnoses_icd9.0."), 
      starts_with("secondary_diagnoses_icd9.0.")),
    ~ substr(., 1, 3) %in% exclude_codes_ICD9
  )
)


# 2. Exclusion - Cholesterol lowering medications

## 2a. Remove female only questions for male participants (vice versa)

ukb <- ukb %>%
  mutate(
    # Fix Male columns
    across(
      c(
        chol_bp_diabetes_hormone_medication.0.0,
        chol_bp_diabetes_hormone_medication.0.1,
        chol_bp_diabetes_hormone_medication.0.2,
        chol_bp_diabetes_hormone_medication.0.3
      ),
      ~ ifelse(as.character(sex) == "Male", NA, as.character(.))
    ),
    
    # Fix Female columns
    across(
      c(
        chol_bp_diabetes_medication.0.0,
        chol_bp_diabetes_medication.0.1,
        chol_bp_diabetes_medication.0.2
      ),
      ~ ifelse(as.character(sex) == "Female", NA, as.character(.))
    )
  )

## 2b. Remove all participants on Cholesterol lowering medications

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


# 3. Exclusion - Participants with cancer within 2 years of recruitment


# 3a. Ensure diagnosis ages/years treated as numeric variables (or converted to NA if applicable)
ukb <- ukb %>%
  mutate(
    across(
      c(
        cancer_first_dx.0.0,
        cancer_first_dx.0.1,
        cancer_first_dx.0.2,
        cancer_first_dx.0.3,
        cancer_first_dx.0.4,
        cancer_first_dx.0.5
      ),
      ~ as.numeric(as.character(.))
    )
  )


# 3b. Create year and age recruited variables

ukb <- ukb %>%
  mutate(
    date_recr = as.Date(date_recr),
    year_recr = as.numeric(format(date_recr, "%Y")),
    age_recr = year_recr - yob
  )
# 3c. Compare diagnosis ages/years to year and age recruited variables
ukb <- ukb %>%
  filter(
    !if_any(
      c(
        cancer_first_dx.0.0,
        cancer_first_dx.0.1,
        cancer_first_dx.0.2,
        cancer_first_dx.0.3,
        cancer_first_dx.0.4,
        cancer_first_dx.0.5
      ),
      ~ !is.na(.) & (
          ( . >= (year_recr - 2) & . <= year_recr ) | 
          ( . >= (age_recr - 2)  & . <= age_recr )
      )
    )
  )
#==========================================================================================







# save
saveRDS(ukb, "ukb_G1_preprocessed.rds")



