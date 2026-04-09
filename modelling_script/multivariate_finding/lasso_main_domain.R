library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(forcats)

# =========================
# Set paths
# =========================
setwd("/rds/general/project/hda_25-26/live/TDS/TDS_Group1")

in_file <- file.path(
  "modelling_script",
  "multivariate_finding",
  "lasso_model",
  "forest_plot_data_lambda_min.csv"
)

out_dir <- file.path(
  "modelling_script",
  "multivariate_finding",
  "forest_plot_by_domain"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

plot_file <- file.path(
  out_dir,
  "forest_plot_lambda_min_penalized_by_domain.png"
)

table_file <- file.path(
  out_dir,
  "forest_plot_lambda_min_penalized_by_domain_table.csv"
)

# =========================
# Helper functions
# =========================
clean_term_label <- function(x) {
  x |>
    str_replace_all("_", " ") |>
    str_replace_all("\n", " ") |>
    str_replace_all("([a-z])([A-Z])", "\\1 \\2") |>
    str_replace_all("([0-9])([A-Za-z])", "\\1 \\2") |>
    str_replace_all("([A-Za-z])([0-9])", "\\1 \\2") |>
    str_squish() |>
    str_replace_all("\\bimd\\b", "IMD") |>
    str_replace_all("\\bbmi\\b", "BMI") |>
    str_replace_all("\\bcvd\\b", "CVD") |>
    str_replace_all("\\bpm10\\b", "PM10") |>
    str_replace_all("\\bpm2 5\\b", "PM2.5") |>
    str_replace_all("\\bno2\\b", "NO2") |>
    str_replace_all("\\btv\\b", "TV")
}

pretty_term_label <- function(term, term_plot = NULL) {
  x <- if (!is.null(term_plot)) term_plot else term
  x <- clean_term_label(x)
  
  case_when(
    str_detect(term, "^famhx") ~ "Family history of CVD",
    str_detect(term, "^bmi$") ~ "BMI",
    str_detect(term, "^temp") ~ "Average temperature",
    str_detect(term, "^pm10") ~ "PM10",
    str_detect(term, "^pm2_5") ~ "PM2.5",
    str_detect(term, "^no2") ~ "NO2",
    str_detect(term, "^neurotic") ~ "Neuroticism",
    str_detect(term, "^sleep") ~ "Sleep quality",
    str_detect(term, "^work_hours") ~ "Work hours/week",
    str_detect(term, "^total_met") ~ "Physical activity",
    str_detect(term, "^fruit") ~ "Fresh fruit intake",
    str_detect(term, "^employment") ~ "Employment",
    str_detect(term, "^sf_score$") ~ "Saturated fat score",
    str_detect(term, "^imd_quintile") ~ x |> str_replace("^IMD quintile", "IMD quintile "),
    str_detect(term, "^salt_3cat") ~ x |>
      str_replace("^salt 3 cat", "Salt intake:") |>
      str_replace("^salt 3cat", "Salt intake:"),
    str_detect(term, "^alcohol_combined") ~ x |>
      str_replace("^alcohol combined", "Alcohol:"),
    str_detect(term, "^redwine_group") ~ x |>
      str_replace("^redwine group", "Red wine:") |>
      str_replace(" drinker$", ""),
    str_detect(term, "^oily_fish|^fish") ~ x |>
      str_replace("^oily fish intake", "Oily fish intake:"),
    str_detect(term, "^urban") ~ x |>
      str_replace("^urban rural 2 cat", "Urban/rural:") |>
      str_replace("^urban rural 2cat", "Urban/rural:"),
    TRUE ~ x
  ) |>
    str_squish()
}

assign_domain <- function(term) {
  case_when(
    str_detect(term, "imd|employment|education|income|deprivation") ~ "Socioeconomic",
    str_detect(term, "alcohol|smoking|salt|fruit|fish|tea|coffee|sleep|tv|met|work_hours|physical|activity|redwine|sf_score") ~ "Lifestyle / Behaviour",
    str_detect(term, "bmi|famhx|family|birth_weight|neurotic|stress|satis") ~ "Clinical / Family",
    str_detect(term, "pm10|pm2_5|no2|noise|greenspace|temp|urban|living") ~ "Environmental",
    TRUE ~ "Other"
  )
}

# =========================
# Load and validate data
# =========================
if (!file.exists(in_file)) {
  stop("Input file not found: ", in_file)
}

df <- read_csv(in_file, show_col_types = FALSE)

required_cols <- c("term", "OR", "low", "high")
missing_cols <- setdiff(required_cols, names(df))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

if (!"term_plot" %in% names(df)) {
  df <- df %>% mutate(term_plot = term)
}

df <- df %>%
  filter(
    !is.na(term),
    !is.na(OR),
    !is.na(low),
    !is.na(high),
    is.finite(OR),
    is.finite(low),
    is.finite(high),
    OR > 0,
    low > 0,
    high > 0
  )

if (nrow(df) == 0) {
  stop("No valid rows available after filtering.")
}

# =========================
# Remove unpenalized/confounder terms
# =========================
df <- df %>%
  filter(
    !str_detect(term, "^age$"),
    !str_detect(term, "^sex"),
    !str_detect(term, "^ethnicity_5cat")
  )

if (nrow(df) == 0) {
  stop("No penalized predictors remain after removing confounders.")
}

# =========================
# Clean labels and assign groups
# =========================
df <- df %>%
  mutate(
    term_label = mapply(pretty_term_label, term, term_plot, USE.NAMES = FALSE),
    domain = assign_domain(term),
    direction = case_when(
      OR > 1 ~ "OR > 1",
      OR < 1 ~ "OR < 1",
      TRUE ~ "OR = 1"
    ),
    effect_size = abs(log(OR)),
    ci_label = sprintf("%.2f (%.3f, %.3f)", OR, low, high)
  )

# =========================
# Order terms within domain
# =========================
domain_order <- c(
  "Socioeconomic",
  "Lifestyle / Behaviour",
  "Clinical / Family",
  "Environmental",
  "Other"
)

df <- df %>%
  mutate(domain = factor(domain, levels = domain_order)) %>%
  arrange(domain, desc(effect_size), desc(OR)) %>%
  group_by(domain) %>%
  mutate(term_label_f = factor(term_label, levels = rev(unique(term_label)))) %>%
  ungroup()

# =========================
# Set plot dimensions
# =========================
x_min_raw <- min(df$low, na.rm = TRUE)
x_max_raw <- max(df$high, na.rm = TRUE)

x_min <- max(0.75, x_min_raw * 0.95)
x_max <- x_max_raw * 1.38
text_x <- x_max_raw * 1.16

plot_width <- 12
plot_height <- max(7, 1.2 + 0.46 * nrow(df))

# =========================
# Create forest plot
# =========================
p <- ggplot(df, aes(x = OR, y = term_label_f)) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.7,
    color = "grey45"
  ) +
  geom_errorbarh(
    aes(xmin = low, xmax = high, color = direction),
    height = 0.14,
    linewidth = 0.95,
    alpha = 0.95
  ) +
  geom_point(
    aes(fill = direction, shape = direction),
    size = 3.7,
    stroke = 0.5,
    color = "black"
  ) +
  geom_text(
    aes(x = text_x, label = ci_label),
    hjust = 0,
    size = 3.8,
    color = "grey20"
  ) +
  facet_grid(
    domain ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  scale_x_log10(
    limits = c(x_min, x_max),
    breaks = c(0.8, 0.9, 1.0, 1.1, 1.2, 1.5),
    labels = c("0.8", "0.9", "1.0", "1.1", "1.2", "1.5")
  ) +
  scale_fill_manual(
    values = c(
      "OR > 1" = "#E76F51",
      "OR < 1" = "#2A9D8F",
      "OR = 1" = "#8D99AE"
    ),
    drop = FALSE
  ) +
  scale_color_manual(
    values = c(
      "OR > 1" = "#E76F51",
      "OR < 1" = "#2A9D8F",
      "OR = 1" = "#8D99AE"
    ),
    drop = FALSE
  ) +
  scale_shape_manual(
    values = c(
      "OR > 1" = 21,
      "OR < 1" = 24,
      "OR = 1" = 22
    ),
    drop = FALSE
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "LASSO-selected penalized predictors at lambda.min",
    subtitle = "Refitted logistic regression excluding unpenalized confounders from display",
    x = "Odds ratio",
    y = NULL,
    fill = NULL,
    color = NULL,
    shape = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(size = 12, color = "grey30"),
    axis.text.y = element_text(size = 11.3, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey88"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text = element_text(size = 11),
    plot.margin = margin(15, 125, 15, 20)
  )

ggsave(
  filename = plot_file,
  plot = p,
  width = plot_width,
  height = plot_height,
  dpi = 320
)

# =========================
# Save output table
# =========================
table_out <- df %>%
  select(domain, term, term_label, OR, low, high, ci_label, direction) %>%
  arrange(domain, desc(abs(log(OR))))

write_csv(table_out, table_file)

cat("Done.\n")
cat("Saved plot to:", plot_file, "\n")
cat("Saved table to:", table_file, "\n")