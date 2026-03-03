suppressPackageStartupMessages({
  library(dplyr)
  library(broom)
})

in_dir  <- "/rds/general/project/hda_25-26/live/TDS/TDS_Group1"
out_dir <- file.path(in_dir, "modelling_script")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
setwd(in_dir)

x <- readRDS("ukb_G1_imputed_final.rds")

outcome <- "cvd_incident"
confounders <- c("age", "sex", "ethnicity_5cat")

exclude_vars <- c(
  "eid","cvd_incident","cvd_event","cvd_prevalent",
  "cvd_first_date","dod","age_of_death","date_recr","yob"
)

predictors <- setdiff(names(x), c(exclude_vars, confounders))

factor_vars <- c(
  "sex","ethnicity_5cat","urban_rural_2cat","famhx_cvd",
  "mood_disorder","diabetes_bin","smoking_3cat",
  "education_2","alcohol_combined","redwine_group",
  "oily_fish_3cat","salt_3cat","living_with_partner",
  "employment_2cat","imd_quintile","medication_category",
  "ICD_Diabetes","ICD_Lipidemia","ICD_Mental_Health",
  "ICD_Migraine","ICD_Hypertension","ICD_AF",
  "ICD_Atopy","ICD_Autoimmune","ICD_CKD",
  "tv_group","total_met_group","stress_group_2yr"
)
factor_vars <- intersect(factor_vars, names(x))
x[factor_vars] <- lapply(x[factor_vars], factor)

max_levels <- 50

run_one <- function(var) {
  if (is.factor(x[[var]]) && nlevels(x[[var]]) > max_levels) return(NULL)
  
  fml <- as.formula(paste0(outcome, " ~ ", var, " + ", paste(confounders, collapse = " + ")))
  fit <- try(glm(fml, data = x, family = binomial()), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  
  td <- broom::tidy(fit, conf.int = TRUE)
  td <- td[td$term == var | grepl(paste0("^", var), td$term), , drop = FALSE]
  if (nrow(td) == 0) return(NULL)
  
  td %>%
    mutate(
      predictor = var,
      OR = exp(estimate),
      CI_low = exp(conf.low),
      CI_high = exp(conf.high)
    )
}

res <- lapply(predictors, run_one)

final_results <- bind_rows(res) %>%
  mutate(p_FDR = p.adjust(p.value, method = "BH")) %>%
  arrange(p_FDR)

write.csv(final_results,
          file.path(out_dir, "univariate_cvd_results.csv"),
          row.names = FALSE)

# =========================
# Visualisations 
# =========================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(stringr)
})

fig_dir <- "/rds/general/user/rh725/projects/hda_25-26/live/TDS/TDS_Group1/modelling_script/figure"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

plot_df <- final_results %>%
  mutate(
    term_label = ifelse(is.na(term), predictor, term)
  ) %>%
  filter(is.finite(OR), is.finite(CI_low), is.finite(CI_high)) %>%
  mutate(term_label = str_replace_all(term_label, ":", "×")) %>%
  arrange(p_FDR)

# -----------------------------
# Top 30 forest plot
# -----------------------------
top_n <- 30

top_df <- plot_df[seq_len(min(top_n, nrow(plot_df))), ] %>%
  mutate(term_label = reorder(term_label, OR))

p_top30 <- ggplot(top_df,
                  aes(x = term_label, y = OR, ymin = CI_low, ymax = CI_high)) +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.6) +
  geom_pointrange(linewidth = 0.8) +
  coord_flip() +
  scale_y_log10() +
  labs(
    title = "Top 30 terms by FDR (Adjusted OR, 95% CI)",
    x = NULL,
    y = "Odds Ratio"
  ) +
  theme_bw(base_size = 14)

# PNG
ggsave(
  filename = file.path(fig_dir, "univariate_forest_top30.png"),
  plot = p_top30,
  width = 12,
  height = 10,
  dpi = 300
)

# PDF
ggsave(
  filename = file.path(fig_dir, "univariate_forest_top30.pdf"),
  plot = p_top30,
  width = 12,
  height = 10
)

# -----------------------------
# All terms forest plot
# -----------------------------
all_df <- plot_df %>%
  mutate(term_label = reorder(term_label, OR))

n_terms <- nrow(all_df)
plot_height <- max(8, n_terms * 0.35)

p_all <- ggplot(all_df,
                aes(x = term_label, y = OR, ymin = CI_low, ymax = CI_high)) +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.6) +
  geom_pointrange(linewidth = 0.7) +
  coord_flip() +
  scale_y_log10() +
  labs(
    title = "Univariate adjusted models (All terms, OR 95% CI)",
    x = NULL,
    y = "Odds Ratio"
  ) +
  theme_bw(base_size = 13)

# PNG
ggsave(
  filename = file.path(fig_dir, "univariate_forest_all_terms.png"),
  plot = p_all,
  width = 12,
  height = plot_height,
  dpi = 300,
  limitsize = FALSE
)

# PDF
ggsave(
  filename = file.path(fig_dir, "univariate_forest_all_terms.pdf"),
  plot = p_all,
  width = 12,
  height = plot_height,
  limitsize = FALSE
)

message("Top 30 and full forest plots saved as PNG and PDF to: ", fig_dir)