# Libraries
library(dplyr)

# path
setwd("/rds/general/project/hda_25-26/live/TDS/TDS_Group1")

# Read data
ukb <- readRDS("ukb_G1_Raw.rds")

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

# save
saveRDS(ukb, "ukb_G1_preprocessed.rds")
