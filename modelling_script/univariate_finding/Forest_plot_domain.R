suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(forcats)
})

# paths
in_dir  <- "/rds/general/project/hda_25-26/live/TDS/TDS_Group1"
out_dir <- file.path(in_dir, "modelling_script", "univariate_finding")
uni_dir <- file.path(out_dir, "univariate_output")

dir.create(uni_dir, showWarnings = FALSE, recursive = TRUE)

# read results
term_df <- read.csv(
  file.path(uni_dir, "univariate_cvd_results_term_level.csv"),
  stringsAsFactors = FALSE
)

pred_df <- read.csv(
  file.path(uni_dir, "univariate_cvd_results_predictor_level.csv"),
  stringsAsFactors = FALSE
)

# predictor -> field mapping
field_map <- data.frame(
  predictor = c(
    # Diet
    "fruit_intake_fresh","tea_intake","coffee_intake",
    "oily_fish_3cat","salt_3cat","alcohol_combined",
    "redwine_group","sf_score",
    
    # Lifestyle
    "tv_group","total_met_group",
    "work_hours_week_clean","smoking_3cat",
    "sleep_quality_score",
    
    # Mental health
    "neuroticism_score","mh_satis_mean_score",
    "mood_disorder","stress_count_2yr","stress_group_2yr",
    
    # Environment
    "no2_2010","pm10_2010","pm2_5_2010","noise_24h",
    "greenspace_pct_1000m","temp_average","urban_rural_2cat",
    
    # Health history
    "diabetes_bin","famhx_cvd",
    "ICD_Diabetes","ICD_Lipidemia","ICD_Hypertension",
    "ICD_AF","ICD_CKD","ICD_Mental_Health","ICD_Migraine",
    "ICD_Atopy","ICD_Autoimmune","medication_category",
    
    # Physical measures
    "sbp_mean","dbp_mean","bmi","creatinine_1","birth_weight_clean",
    
    # Socio-demographics
    "education_2","employment_2cat",
    "living_with_partner","imd_final",
    "imd_source","imd_quintile"
  ),
  field = c(
    rep("Diet", 8),
    rep("Lifestyle", 5),
    rep("Mental health", 5),
    rep("Environment", 7),
    rep("Health and medical history", 12),
    rep("Physical measures", 5),
    rep("Socio-demographics", 6)
  ),
  stringsAsFactors = FALSE
)

# biomarker domains
blood_biomarkers <- c(
  "wbc_count","rbc_count","hemoglobin","mcv","platelet_count","mpv",
  "lymphocyte_count","monocyte_count","neutrophil_count",
  "eosinophil_count","basophil_count","nrbc_count"
)

biochem_biomarkers <- c(
  "creatinine","cystatin_c","urea","urate","microalbumin",
  "albumin","alkaline_phosphatase",
  "alanine_aminotransferase","aspartate_aminotransferase",
  "direct_bilirubin","total_bilirubin","total_blood_protein",
  "cholesterol","hdl_cholesterol","ldl_cholesterol",
  "lipoprotein_a","apolipoprotein_a","apolipoprotein_b",
  "total_triglyceride",
  "glucose","hba1c","igf_1","crp","blood_vitamin_d",
  "calcium","phosphate"
)

biomarker_map <- bind_rows(
  data.frame(
    predictor = blood_biomarkers,
    field = "Blood biomarkers",
    stringsAsFactors = FALSE
  ),
  data.frame(
    predictor = biochem_biomarkers,
    field = "Biochemistry biomarkers",
    stringsAsFactors = FALSE
  )
)

# extra predictors
extra_map <- data.frame(
  predictor = c(
    "shift_work","mixed_shift","job_shift_work_clean","ever_nights","work_hours_unified_cat",
    "household_size_clean","living_alone",
    "self_health_bin","gp_anxdep","hrt"
  ),
  field = c(
    rep("Lifestyle", 5),
    rep("Socio-demographics", 2),
    rep("Health and medical history", 3)
  ),
  stringsAsFactors = FALSE
)

field_map <- bind_rows(field_map, biomarker_map, extra_map) %>%
  distinct(predictor, .keep_all = TRUE)

# automatically add predictors not in field_map
all_predictors <- unique(c(term_df$predictor, pred_df$predictor))
missing_predictors <- setdiff(all_predictors, field_map$predictor)

if (length(missing_predictors) > 0) {
  
  auto_field <- function(x) {
    if (grepl("creatinine|cholesterol|glucose|hba1c|albumin|microalbumin|bilirubin|protein|triglyceride|apolipoprotein|cystatin|urea|urate|crp|vitamin_d|calcium|phosphate|igf", x, ignore.case = TRUE))
      return("Biochemistry biomarkers")
    if (grepl("wbc|rbc|hemoglobin|mcv|platelet|lymphocyte|monocyte|neutrophil|eosinophil|basophil|nrbc", x, ignore.case = TRUE))
      return("Blood biomarkers")
    if (grepl("imd|education|employment|partner|household|living", x, ignore.case = TRUE))
      return("Socio-demographics")
    if (grepl("smoking|alcohol|sleep|tv|met|work|shift|night", x, ignore.case = TRUE))
      return("Lifestyle")
    if (grepl("neuroticism|stress|mood|satis|mental|anxdep", x, ignore.case = TRUE))
      return("Mental health")
    if (grepl("no2|pm|noise|greenspace|temp|urban", x, ignore.case = TRUE))
      return("Environment")
    if (grepl("diabetes|hypertension|af|ckd|migraine|autoimmune|atopy|medication|hrt|self_health|health_bin", x, ignore.case = TRUE))
      return("Health and medical history")
    return("Other")
  }
  
  auto_map <- data.frame(
    predictor = missing_predictors,
    field = sapply(missing_predictors, auto_field),
    stringsAsFactors = FALSE
  )
  
  field_map <- bind_rows(field_map, auto_map) %>%
    distinct(predictor, .keep_all = TRUE)
}

field_order <- c(
  "Blood biomarkers",
  "Biochemistry biomarkers",
  "Diet",
  "Environment",
  "Health and medical history",
  "Lifestyle",
  "Mental health",
  "Physical measures",
  "Socio-demographics",
  "Other"
)

# term-level prep
term_plot_df <- term_df %>%
  left_join(field_map, by = "predictor") %>%
  mutate(
    field = ifelse(is.na(field), "Other", field),
    field = factor(field, levels = field_order),
    label = term,
    OR_plot = OR,
    OR_low = OR_CI_low,
    OR_high = OR_CI_high,
    p_plot = raw_p_value,
    bonf_plot = Bonferroni_p_value,
    bonf_sig = !is.na(bonf_plot) & bonf_plot < 0.05
  ) %>%
  filter(
    is.finite(OR_plot), is.finite(OR_low), is.finite(OR_high),
    OR_plot > 0, OR_low > 0, OR_high > 0
  ) %>%
  mutate(effect = abs(log(OR_plot))) %>%
  group_by(field) %>%
  arrange(desc(bonf_sig), p_plot, desc(effect), .by_group = TRUE) %>%
  slice_head(n = 15) %>%
  ungroup()

term_plot_df <- term_plot_df %>%
  group_by(field) %>%
  mutate(label = fct_reorder(label, OR_plot)) %>%
  ungroup()

# predictor-level prep
pred_plot_df <- pred_df %>%
  left_join(field_map, by = "predictor") %>%
  mutate(
    field = ifelse(is.na(field), "Other", field),
    field = factor(field, levels = field_order),
    label = predictor,
    OR_plot = top_OR,
    OR_low = top_OR_CI_low,
    OR_high = top_OR_CI_high,
    p_plot = predictor_raw_p_value,
    bonf_plot = predictor_Bonferroni_p_value,
    bonf_sig = !is.na(bonf_plot) & bonf_plot < 0.05
  ) %>%
  filter(
    is.finite(OR_plot), is.finite(OR_low), is.finite(OR_high),
    OR_plot > 0, OR_low > 0, OR_high > 0
  ) %>%
  mutate(effect = abs(log(OR_plot))) %>%
  group_by(field) %>%
  arrange(desc(bonf_sig), p_plot, desc(effect), .by_group = TRUE) %>%
  slice_head(n = 15) %>%
  ungroup()

pred_plot_df <- pred_plot_df %>%
  group_by(field) %>%
  mutate(label = fct_reorder(label, OR_plot)) %>%
  ungroup()

# plotting helper
make_forest_plot <- function(plot_df, title_suffix = "") {
  ggplot(plot_df, aes(x = OR_plot, y = label)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_errorbarh(aes(xmin = OR_low, xmax = OR_high), height = 0.2, colour = "grey45") +
    geom_point(aes(colour = bonf_sig), size = 2.5) +
    scale_colour_manual(
      values = c("TRUE" = "red", "FALSE" = "grey40"),
      name = "Bonferroni significant"
    ) +
    scale_x_log10(
      breaks = c(0.5, 1, 3, 5),
      limits = c(0.5, 5)
    ) +
    facet_wrap(~field, scales = "free_y", ncol = 3, drop = TRUE) +
    labs(
      title = "Univariate associations with has_cvd (logistic regression)",
      subtitle = paste0("Odds ratios with 95% CI; red = Bonferroni-adjusted p < 0.05", title_suffix),
      x = "Odds ratio",
      y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
}

# build plots
p_term <- make_forest_plot(term_plot_df, " [term level]")
p_pred <- make_forest_plot(pred_plot_df, " [predictor level]")

# save plots
ggsave(
  file.path(uni_dir, "univariate_forest_term_level_by_field.png"),
  p_term, width = 14, height = 9, dpi = 300
)

ggsave(
  file.path(uni_dir, "univariate_forest_term_level_by_field.pdf"),
  p_term, width = 14, height = 9
)

ggsave(
  file.path(uni_dir, "univariate_forest_predictor_level_by_field.png"),
  p_pred, width = 14, height = 9, dpi = 300
)

ggsave(
  file.path(uni_dir, "univariate_forest_predictor_level_by_field.pdf"),
  p_pred, width = 14, height = 9
)