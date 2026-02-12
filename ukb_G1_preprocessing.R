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

# quick checks
if ("sys_bp_automatic" %in% names(ukb)) print(summary(ukb$sys_bp_automatic))
if ("sbp_mean" %in% names(ukb)) print(summary(ukb$sbp_mean))
if ("dbp_mean" %in% names(ukb)) print(summary(ukb$dbp_mean))

##############################################################################
# recode ethnicity -> 5 groups (White / South Asian / Chinese / Black / Other)
# "Prefer not to answer" and "Do not know" -> missing (NA)
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

# quick checks
print(table(ukb$ethnicity_5cat, useNA = "ifany"))
print(table(ukb$ethnicity, ukb$ethnicity_5cat, useNA = "ifany"))

##############################################################################
# recode urban_rural -> 2 groups (Urban / Rural)
# "Postcode not linkable" -> missing (NA)
##############################################################################
x <- ukb$urban_rural

ur2 <- case_when(
  grepl("Postcode not linkable", x) ~ NA_character_,
  grepl("Urban", x) ~ "Urban",
  grepl("Town and Fringe|Village|Hamlet|Rural|Small Town", x) ~ "Rural",
  TRUE ~ NA_character_
)

ukb$urban_rural_2cat <- factor(ur2, levels = c("Urban", "Rural"))

# quick checks
print(table(ukb$urban_rural_2cat, useNA = "ifany"))
print(table(ukb$urban_rural, ukb$urban_rural_2cat, useNA = "ifany"))

##############################################################################
# Air pollution cleanup (2010): remove NO, keep NO2; drop PM10, keep PM2.5
# - set negative values to NA (no physical meaning)
##############################################################################

# NO2: keep, set negatives to NA
ukb$no2_2010[ukb$no2_2010 < 0] <- NA

# NO: drop (highly collinear with NO2)
if ("no_2010" %in% names(ukb)) ukb$no_2010 <- NULL

# PM2.5: keep
# PM10: drop (moderate correlation with PM2.5; not averaged)
if ("pm10_2010" %in% names(ukb)) ukb$pm10_2010 <- NULL

# quick checks
if ("no2_2010" %in% names(ukb)) print(summary(ukb$no2_2010))
if ("pm2_5_2010" %in% names(ukb)) print(summary(ukb$pm2_5_2010))

##############################################################################
# Greenspace / Water / Noise cleanup (1000m + 24h)
# - greenspace_pct_1000m and water_pct_1000m: negative values -> NA
# - greenspace_pct_1000m: >100 -> NA (percent out of range)
##############################################################################

# greenspace (%)
if ("greenspace_pct_1000m" %in% names(ukb)) {
  ukb$greenspace_pct_1000m[ukb$greenspace_pct_1000m < 0] <- NA
  ukb$greenspace_pct_1000m[ukb$greenspace_pct_1000m > 100] <- NA
}

# water (%)
if ("water_pct_1000m" %in% names(ukb)) {
  ukb$water_pct_1000m[ukb$water_pct_1000m < 0] <- NA
}

# noise_24h (no cleaning needed based on your summary, but keep a quick check)
if ("noise_24h" %in% names(ukb)) {
  print(summary(ukb$noise_24h))
}

# quick checks
if ("greenspace_pct_1000m" %in% names(ukb)) print(summary(ukb$greenspace_pct_1000m))
if ("water_pct_1000m" %in% names(ukb)) print(summary(ukb$water_pct_1000m))

##############################################################################
# Biomarker cleanup
# - keep blood_sample_attempted as factor (metadata)
# - drop urine_device (no variation)
# - convert biomarker strings to numeric
# - set negative values to NA (UKB special codes)
##############################################################################

# metadata
if ("blood_sample_attempted" %in% names(ukb)) {
  ukb$blood_sample_attempted <- factor(ukb$blood_sample_attempted)
  print(table(ukb$blood_sample_attempted, useNA = "ifany"))
}

if ("urine_device" %in% names(ukb)) {
  print(table(ukb$urine_device, useNA = "ifany"))
  if (length(unique(na.omit(ukb$urine_device))) <= 1) ukb$urine_device <- NULL
}

# biomarker columns (numeric)
biomarkers_num <- c(
  "sodium_in_urine","creatinine","microalbumin",
  "wbc_count","rbc_count","hemoglobin","mcv","platelet_count","mpv",
  "lymphocyte_count","monocyte_count","neutrophil_count",
  "eosinophil_count","basophil_count","nrbc_count",
  "albumin","alkaline_phosphatase","alanine_aminotransferase",
  "aspartate_aminotransferase","direct_bilirubin","urea","calcium",
  "cystatin_c","glucose","hba1c","igf_1","phosphate",
  "total_bilirubin","total_blood_protein","total_triglyceride",
  "urate","blood_vitamin_d","crp","lipoprotein_a"
)

biomarkers_num <- biomarkers_num[biomarkers_num %in% names(ukb)]

# convert character -> numeric
for (v in biomarkers_num) {
  if (is.character(ukb[[v]])) ukb[[v]] <- as.numeric(ukb[[v]])
}

# negative values -> NA
for (v in biomarkers_num) {
  if (is.numeric(ukb[[v]])) ukb[[v]][ukb[[v]] < 0] <- NA
}

# quick checks
print(sapply(biomarkers_num, function(v) mean(is.na(ukb[[v]]))))

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
      .data[[alcohol_status]] == "Never" & .data[[alcohol_freq]] != "Never" ~ NA_character_,
      .data[[alcohol_status]] != "Never" & .data[[alcohol_freq]] == "Never" ~ NA_character_,
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

#################################YILIU########################################


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


