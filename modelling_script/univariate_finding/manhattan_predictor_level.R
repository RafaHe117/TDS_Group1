suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(ggrepel)
})

# ---- directories ----
in_dir  <- "/rds/general/project/hda_25-26/live/TDS/TDS_Group1"
out_dir <- file.path(in_dir, "modelling_script", "univariate_finding")
uni_dir <- file.path(out_dir, "univariate_output")

dir.create(uni_dir, showWarnings = FALSE, recursive = TRUE)
setwd(in_dir)

# ---- load results ----
final_results_predictor <- read.csv(
  file.path(uni_dir, "univariate_cvd_results_predictor_level.csv"),
  stringsAsFactors = FALSE
)

multiple_testing_summary <- read.csv(
  file.path(uni_dir, "multiple_testing_summary.csv"),
  stringsAsFactors = FALSE
)

# ---- base predictor-domain mapping ----
domain_map_base <- data.frame(
  predictor = c(
    "fruit_intake_fresh","tea_intake","coffee_intake",
    "oily_fish_3cat","salt_3cat","alcohol_combined","redwine_group",
    "sf_score",
    "tv_group","total_met_group",
    "work_hours_week_clean","smoking_3cat",
    "neuroticism_score","sleep_quality_score","mh_satis_mean_score",
    "mood_disorder","stress_count_2yr","stress_group_2yr",
    "no2_2010","pm10_2010","pm2_5_2010","noise_24h",
    "greenspace_pct_1000m","temp_average","urban_rural_2cat",
    "wbc_count","rbc_count","hemoglobin","mcv","platelet_count","mpv",
    "lymphocyte_count","monocyte_count","neutrophil_count",
    "eosinophil_count","basophil_count",
    "creatinine","cystatin_c","urea","urate",
    "albumin","alkaline_phosphatase","alanine_aminotransferase",
    "aspartate_aminotransferase","direct_bilirubin","total_bilirubin",
    "total_blood_protein",
    "cholesterol","hdl_cholesterol","ldl_cholesterol","lipoprotein_a",
    "apolipoprotein_a","apolipoprotein_b","total_triglyceride",
    "glucose","hba1c","igf_1","crp","blood_vitamin_d",
    "calcium","phosphate",
    "sbp_mean","dbp_mean",
    "bmi","birth_weight_clean","diabetes_bin","famhx_cvd",
    "ICD_Diabetes","ICD_Lipidemia","ICD_Hypertension","ICD_AF",
    "ICD_CKD","ICD_Mental_Health","ICD_Migraine","ICD_Atopy",
    "ICD_Autoimmune",
    "education_2","employment_2cat","living_with_partner",
    "imd_final","imd_source","imd_quintile","medication_category"
  ),
  domain = c(
    rep("Diet", 8),
    rep("Lifestyle", 4),
    rep("Psychosocial", 6),
    rep("Environment", 7),
    rep("Haematology", 11),
    rep("Renal", 4),
    rep("Liver_Protein", 7),
    rep("Lipid", 7),
    rep("Metabolic_Inflammation", 7),
    rep("Blood_Pressure", 2),
    rep("Clinical_History", 13),
    rep("Socioeconomic", 7)
  ),
  stringsAsFactors = FALSE
)

# ---- additional predictors (safe if absent in results) ----
domain_map_extra <- data.frame(
  predictor = c(
    "shift_work",
    "mixed_shift",
    "household_size_clean",
    "living_alone",
    "job_shift_work_clean",
    "ever_nights",
    "nrbc_count",
    "microalbumin",
    "hrt",
    "self_health_bin",
    "gp_anxdep",
    "work_hours_unified_cat",
    "creatinine_1"
  ),
  domain = c(
    "Lifestyle",
    "Lifestyle",
    "Socioeconomic",
    "Socioeconomic",
    "Lifestyle",
    "Lifestyle",
    "Haematology",
    "Renal",
    "Clinical_History",
    "Clinical_History",
    "Clinical_History",
    "Lifestyle",
    "Renal"
  ),
  stringsAsFactors = FALSE
)

# combine mappings and remove duplicates
domain_map <- bind_rows(domain_map_base, domain_map_extra) %>%
  distinct(predictor, .keep_all = TRUE)

# ---- prepare data for Manhattan plot ----
manhattan_df <- final_results_predictor %>%
  mutate(
    predictor_raw_p_value = pmax(predictor_raw_p_value, .Machine$double.xmin)
  ) %>%
  left_join(domain_map, by = "predictor") %>%
  mutate(
    domain = ifelse(is.na(domain), "Other", domain),
    term_label = case_when(
      !is.na(top_term) & top_term != "" ~ top_term,
      predictor == "sf_score" ~ "saturated_score",
      TRUE ~ predictor
    ),
    neg_log10_p = -log10(predictor_raw_p_value)
  ) %>%
  filter(is.finite(neg_log10_p))

# ---- domain ordering ----
domain_order <- c(
  "Diet","Lifestyle","Psychosocial","Environment",
  "Haematology","Renal","Liver_Protein","Lipid",
  "Metabolic_Inflammation","Blood_Pressure",
  "Clinical_History","Socioeconomic","Other"
)

manhattan_df <- manhattan_df %>%
  mutate(domain = factor(domain, levels = domain_order)) %>%
  arrange(domain, predictor_raw_p_value, predictor, term_label) %>%
  mutate(index = row_number())

# ---- compute domain bands ----
domain_summary <- manhattan_df %>%
  group_by(domain) %>%
  summarise(
    xmin = min(index) - 0.5,
    xmax = max(index) + 0.5,
    center = mean(index),
    .groups = "drop"
  ) %>%
  mutate(
    band_id = row_number(),
    band_fill = ifelse(band_id %% 2 == 1, "grey98", "grey94")
  )

domain_boundaries <- domain_summary %>%
  arrange(xmin) %>%
  slice(-n()) %>%
  transmute(xintercept = xmax)

# ---- significance thresholds ----
nominal_line <- -log10(0.05)

bonf_threshold <- multiple_testing_summary %>%
  filter(level == "predictor_level") %>%
  pull(bonferroni_threshold)

if (length(bonf_threshold) != 1 || is.na(bonf_threshold)) {
  stop("Could not uniquely read predictor-level Bonferroni threshold from multiple_testing_summary.csv")
}

bonf_line <- -log10(pmax(bonf_threshold, .Machine$double.xmin))

# ---- threshold dataframe for legend ----
threshold_df <- data.frame(
  type = c("Nominal p = 0.05", "Bonferroni"),
  y = c(nominal_line, bonf_line),
  stringsAsFactors = FALSE
)

# ---- cap y-axis for readability ----
y_cap <- 12

manhattan_df <- manhattan_df %>%
  mutate(
    neg_log10_p_plot = pmin(neg_log10_p, y_cap),
    is_truncated = neg_log10_p > y_cap
  )

# ---- label predictors with raw p < 0.05 ----
label_df <- manhattan_df %>%
  filter(predictor_raw_p_value < 0.05) %>%
  distinct(predictor, .keep_all = TRUE) %>%
  mutate(
    label = ifelse(
      is_truncated,
      paste0(predictor, " (", round(neg_log10_p, 1), ")"),
      predictor
    )
  )

# ---- domain colours ----
domain_cols <- c(
  "Diet" = "#E76F51",
  "Lifestyle" = "#F4A261",
  "Psychosocial" = "#C9A227",
  "Environment" = "#7CAE00",
  "Haematology" = "#2FBF71",
  "Renal" = "#2EC4B6",
  "Liver_Protein" = "#1FBFCF",
  "Lipid" = "#27A7E7",
  "Metabolic_Inflammation" = "#5C88FF",
  "Blood_Pressure" = "#B983FF",
  "Clinical_History" = "#F15BB5",
  "Socioeconomic" = "#FF6FAE",
  "Other" = "grey50"
)

threshold_cols <- c(
  "Nominal p = 0.05" = "grey50",
  "Bonferroni" = "#C44E52"
)

threshold_lty <- c(
  "Nominal p = 0.05" = "dashed",
  "Bonferroni" = "dashed"
)

# ---- build Manhattan plot ----
p_manhattan <- ggplot(manhattan_df, aes(x = index, y = neg_log10_p_plot)) +
  geom_rect(
    data = domain_summary,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = band_fill),
    inherit.aes = FALSE,
    alpha = 0.7,
    color = NA
  ) +
  scale_fill_identity() +
  geom_vline(
    data = domain_boundaries,
    aes(xintercept = xintercept),
    inherit.aes = FALSE,
    linewidth = 0.55,
    color = "grey70"
  ) +
  geom_hline(
    data = threshold_df,
    aes(yintercept = y, color = type, linetype = type),
    linewidth = 0.5,
    alpha = 0.85,
    show.legend = TRUE
  ) +
  geom_point(aes(color = domain), size = 2.2, alpha = 0.9, show.legend = TRUE) +
  geom_point(
    data = manhattan_df %>% filter(is_truncated),
    aes(x = index, y = y_cap),
    inherit.aes = FALSE,
    shape = 24,
    size = 3,
    fill = "black"
  ) +
  geom_point(
    data = label_df,
    aes(x = index, y = neg_log10_p_plot),
    inherit.aes = FALSE,
    size = 2.8,
    color = "black"
  ) +
  ggrepel::geom_text_repel(
    data = label_df,
    aes(x = index, y = neg_log10_p_plot, label = label),
    size = 3.2,
    box.padding = 0.45,
    point.padding = 0.3,
    segment.size = 0.28,
    max.overlaps = Inf,
    seed = 42
  ) +
  scale_color_manual(
    values = c(domain_cols, threshold_cols),
    breaks = names(threshold_cols),
    labels = names(threshold_cols)
  ) +
  scale_linetype_manual(
    values = threshold_lty,
    breaks = names(threshold_lty),
    labels = names(threshold_lty)
  ) +
  scale_x_continuous(
    breaks = domain_summary$center,
    labels = domain_summary$domain,
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  coord_cartesian(ylim = c(0, y_cap), clip = "off") +
  labs(
    title = "Univariate associations across predictor domains",
    subtitle = "Predictors with raw p < 0.05 labelled",
    x = NULL,
    y = expression(-log[10](italic(p))),
    color = NULL,
    linetype = NULL
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        shape = NA,
        linewidth = 0.9,
        alpha = 1
      ),
      order = 1
    ),
    linetype = guide_legend(
      override.aes = list(
        color = unname(threshold_cols),
        linewidth = 0.9
      ),
      order = 1
    )
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1, face = "bold"),
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.box.background = element_rect(color = "grey70", fill = "white"),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key.width = unit(1.6, "cm"),
    legend.text = element_text(size = 10),
    plot.margin = margin(10, 24, 12, 10)
  )

# ---- save plots ----
ggsave(
  filename = file.path(uni_dir, "univariate_manhattan_by_domain.png"),
  plot = p_manhattan,
  width = 16,
  height = 8,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(uni_dir, "univariate_manhattan_by_domain.pdf"),
  plot = p_manhattan,
  width = 16,
  height = 8,
  bg = "white"
)