library(dplyr)
library(ggplot2)
library(stringr)
library(forcats)
library(ggh4x)
library(grid)
library(ggtext)
library(htmltools)

# =========================
# Base directory
# =========================
base_out_dir <- "/rds/general/project/hda_25-26/live/TDS/TDS_Group1/modelling_script/stability_analysis/subsample_lasso_sex"

female_csv <- file.path(base_out_dir, "female", "selection_proportions_penalized_female.csv")
male_csv   <- file.path(base_out_dir, "male", "selection_proportions_penalized_male.csv")

female_summary <- file.path(base_out_dir, "female", "summary_female.txt")
male_summary   <- file.path(base_out_dir, "male", "summary_male.txt")

# =========================
# Read calibrated pi threshold
# =========================
read_pi_from_summary <- function(path) {
  lines <- readLines(path, warn = FALSE)
  pi_line <- lines[grepl("^Best pi threshold:", lines)]
  
  if (length(pi_line) == 0) {
    stop(paste("Could not find 'Best pi threshold' in:", path))
  }
  
  as.numeric(sub("^Best pi threshold:\\s*", "", pi_line[1]))
}

female_pi <- read_pi_from_summary(female_summary)
male_pi   <- read_pi_from_summary(male_summary)

# =========================
# Read data
# =========================
female_df <- read.csv(female_csv, stringsAsFactors = FALSE)
male_df   <- read.csv(male_csv, stringsAsFactors = FALSE)

# =========================
# Clean variable names
# =========================
make_pretty_names <- function(x) {
  out <- x
  
  out <- str_replace(out, "^tv_group", "TV time: ")
  out <- str_replace(out, "^total_met_group", "Physical activity: ")
  out <- str_replace(out, "^smoking_3cat", "Smoking status: ")
  out <- str_replace(out, "^alcohol_combined", "Alcohol intake: ")
  out <- str_replace(out, "^redwine_group", "Redwine intake: ")
  out <- str_replace(out, " drinker$", "")
  out <- str_replace(out, "^sf_score$", "Saturated fat score")
  out <- str_replace(out, "^salt_3cat", "Salt intake: ")
  out <- str_replace(out, "^education_2", "Education: ")
  out <- str_replace(out, "^employment_2cat", "Employment: ")
  out <- str_replace(out, "^imd_quintile", "IMD quintile: ")
  out <- str_replace(out, "^oily_fish_3cat", "Oily fish intake: ")
  out <- str_replace(out, "^urban_rural_2cat", "Urban/rural: ")
  out <- str_replace(out, "^living_with_partner", "Living with partner: ")
  out <- str_replace(out, "^stress_group_2yr", "Stress group: ")
  
  out <- str_replace(out, "^fruit_intake_fresh$", "Fresh fruit intake")
  out <- str_replace(out, "^tea_intake$", "Tea intake")
  out <- str_replace(out, "^coffee_intake$", "Coffee intake")
  out <- str_replace(out, "^sleep_quality_score$", "Sleep quality score")
  out <- str_replace(out, "^stress_count_2yr$", "Stress count (2 years)")
  out <- str_replace(out, "^mh_satis_mean_score$", "satisfaction score")
  out <- str_replace(out, "^neuroticism_score$", "Neuroticism score")
  out <- str_replace(out, "^work_hours_week_clean$", "Work hours per week")
  out <- str_replace(out, "^greenspace_pct_1000m$", "Greenspace within 1000 m")
  out <- str_replace(out, "^temp_average$", "Average temperature")
  out <- str_replace(out, "^noise_24h$", "24 h noise")
  out <- str_replace(out, "^pm10_2010$", "PM10 (2010)")
  out <- str_replace(out, "^pm2_5_2010$", "PM2.5 (2010)")
  out <- str_replace(out, "^no2_2010$", "NO2 (2010)")
  out <- str_replace(out, "^bmi$", "BMI")
  
  out |>
    str_replace_all("_", " ") |>
    str_replace_all("pm2 5", "PM2.5") |>
    str_replace_all("pm10", "PM10") |>
    str_replace_all("no2", "NO2") |>
    str_replace_all("imd", "IMD") |>
    str_replace_all("bmi", "BMI") |>
    str_replace_all("tv", "TV") |>
    str_replace_all("cvd", "CVD") |>
    str_replace_all("3cat", "") |>
    str_replace_all("2cat", "") |>
    str_replace_all("Current-", "") |>
    str_replace_all(">", "> ") |>
    str_replace_all("鈥", "–") |>
    str_squish()
}

# =========================
# Assign domain
# =========================
assign_domain <- function(term) {
  case_when(
    str_detect(term, "imd|employment|education|income|deprivation") ~ "Socioeconomic",
    str_detect(term, "alcohol|smoking|salt|fruit|fish|tea|coffee|sleep|tv|met|work_hours|physical|activity|redwine|sf_score") ~ "Lifestyle / Behaviour",
    str_detect(term, "bmi|famhx|family|birth_weight|neurotic|stress|satis") ~ "Clinical / Family",
    str_detect(term, "pm10|pm2_5|no2|noise|greenspace|temp|urban|living") ~ "Environmental",
    TRUE ~ "Other"
  )
}

domain_levels <- c(
  "Socioeconomic",
  "Lifestyle / Behaviour",
  "Clinical / Family",
  "Environmental",
  "Other"
)

# =========================
# Domain strip colours
# =========================
domain_fill <- c(
  "Socioeconomic" = "#F3E9DC",
  "Lifestyle / Behaviour" = "#E6F0EB",
  "Clinical / Family" = "#ECE6F2",
  "Environmental" = "#E5EEF7",
  "Other" = "#F2F2F2"
)

# =========================
# Infer reference groups from terms
# =========================
all_terms <- union(female_df$term, male_df$term)

normalize_level <- function(x) {
  x |>
    str_replace_all("–", "-") |>
    str_replace_all("—", "-") |>
    str_squish() |>
    str_to_lower()
}

infer_single_reference <- function(present_levels, all_levels) {
  present_norm <- normalize_level(present_levels)
  all_norm <- normalize_level(all_levels)
  
  ref_norm <- setdiff(all_norm, present_norm)
  
  if (length(ref_norm) == 1) {
    all_levels[match(ref_norm, all_norm)]
  } else {
    NA_character_
  }
}

infer_reference_groups <- function(terms) {
  ref_list <- list()
  
  # TV time
  tv_present <- terms[str_detect(terms, "^tv_group")] |>
    str_remove("^tv_group") |>
    unique()
  tv_ref <- infer_single_reference(tv_present, c("<1h", "1-3h", "3-5h", ">5h"))
  if (!is.na(tv_ref)) {
    ref_list[["TV time"]] <- ifelse(tv_ref == "<1h", "<1 h/day", tv_ref)
  }
  
  # Physical activity
  met_present <- terms[str_detect(terms, "^total_met_group")] |>
    str_remove("^total_met_group") |>
    unique()
  met_ref <- infer_single_reference(met_present, c("Low", "Moderate", "High"))
  if (!is.na(met_ref)) {
    ref_list[["Physical activity"]] <- str_to_lower(met_ref)
  }
  
  # Smoking
  smoking_present <- terms[str_detect(terms, "^smoking_3cat")] |>
    str_remove("^smoking_3cat") |>
    unique()
  smoking_ref <- infer_single_reference(smoking_present, c("Never", "Previous", "Current"))
  if (!is.na(smoking_ref)) {
    ref_list[["Smoking"]] <- str_to_lower(smoking_ref)
  }
  
  # Alcohol
  alcohol_present <- terms[str_detect(terms, "^alcohol_combined")] |>
    str_remove("^alcohol_combined") |>
    unique()
  alcohol_ref <- infer_single_reference(
    alcohol_present,
    c("Current-High","Current-Moderate","Current-Low","Prefer not to answer","Never","Former")
  )
  if (!is.na(alcohol_ref)) {
    ref_list[["Alcohol"]] <- alcohol_ref
  }
  
  # Salt
  salt_present <- terms[str_detect(terms, "^salt_3cat")] |>
    str_remove("^salt_3cat") |>
    unique()
  salt_ref <- infer_single_reference(salt_present, c("Low", "Medium", "High"))
  if (!is.na(salt_ref)) {
    ref_list[["Salt"]] <- str_to_lower(salt_ref)
  }
  
  # Education
  edu_present <- terms[str_detect(terms, "^education_2")] |>
    str_remove("^education_2") |>
    unique()
  edu_ref <- infer_single_reference(edu_present, c("<=12Yrs", "12-14Yrs", ">15Yrs"))
  if (!is.na(edu_ref)) {
    ref_list[["Education"]] <- edu_ref
  }
  
  # Employment
  emp_present <- terms[str_detect(terms, "^employment_2cat")] |>
    str_remove("^employment_2cat") |>
    unique()
  emp_ref <- infer_single_reference(emp_present, c("Not in paid employment", "In paid employment"))
  if (!is.na(emp_ref)) {
    ref_list[["Employment"]] <- emp_ref
  }
  
  # IMD quintile
  imd_present <- terms[str_detect(terms, "^imd_quintile")] |>
    str_remove("^imd_quintile") |>
    unique()
  imd_ref <- infer_single_reference(imd_present, c("1", "2", "3", "4", "5"))
  if (!is.na(imd_ref)) {
    ref_list[["IMD quintile"]] <- imd_ref
  }
  
  # Oily fish
  fish_present <- terms[str_detect(terms, "^oily_fish_3cat")] |>
    str_remove("^oily_fish_3cat") |>
    unique()
  fish_ref <- infer_single_reference(fish_present, c("Low", "Medium", "High"))
  if (!is.na(fish_ref)) {
    ref_list[["Oily fish"]] <- str_to_lower(fish_ref)
  }
  
  # Urban / rural
  urban_present <- terms[str_detect(terms, "^urban_rural_2cat")] |>
    str_remove("^urban_rural_2cat") |>
    unique()
  urban_ref <- infer_single_reference(urban_present, c("Urban", "Rural"))
  if (!is.na(urban_ref)) {
    ref_list[["Urban/rural"]] <- urban_ref
  }
  
  # Living with partner
  partner_present <- terms[str_detect(terms, "^living_with_partner")] |>
    str_remove("^living_with_partner") |>
    unique()
  partner_ref <- infer_single_reference(partner_present, c("No", "Yes"))
  if (!is.na(partner_ref)) {
    ref_list[["Living with partner"]] <- partner_ref
  }
  
  # Red wine
  redwine_present <- terms[str_detect(terms, "^redwine_group")] |>
    str_remove("^redwine_group") |>
    unique()
  redwine_ref <- infer_single_reference(
    redwine_present,
    c("Non-drinker", "Light drinker", "Moderate drinker", "Heavy drinker")
  )
  if (!is.na(redwine_ref)) {
    ref_list[["Red wine intake"]] <- redwine_ref
  }
  
  # Stress group
  stress_present <- terms[str_detect(terms, "^stress_group_2yr")] |>
    str_remove("^stress_group_2yr") |>
    unique()
  stress_ref <- infer_single_reference(stress_present, c("0-2", "3-4", "5-6"))
  if (!is.na(stress_ref)) {
    ref_list[["Stress group"]] <- stress_ref
  }
  
  ref_list
}

reference_groups <- infer_reference_groups(all_terms)

reference_text <- if (length(reference_groups) > 0) {
  ref_items <- paste(names(reference_groups), unlist(reference_groups), sep = " = ")
  
  n_items <- length(ref_items)
  split_idx <- ceiling(n_items / 2)
  
  line1 <- paste(ref_items[1:split_idx], collapse = "; ")
  line2 <- paste(ref_items[(split_idx + 1):n_items], collapse = "; ")
  
  paste0("Reference groups: ", line1, ";\n", line2, ".")
} else {
  "Reference groups could not be inferred from the available CSV terms."
}
# =========================
# Merge + reorder by domain
# =========================
plot_compare <- full_join(
  female_df |>
    transmute(
      term,
      term_plot = make_pretty_names(term),
      selection_female = selection_proportion
    ),
  male_df |>
    transmute(
      term,
      term_plot = make_pretty_names(term),
      selection_male = selection_proportion
    ),
  by = c("term", "term_plot")
) |>
  mutate(
    selection_female = ifelse(is.na(selection_female), 0, selection_female),
    selection_male   = ifelse(is.na(selection_male), 0, selection_male)
  ) |>
  filter(!(selection_female == 0 & selection_male == 0)) |>
  mutate(
    domain = assign_domain(term),
    domain = factor(domain, levels = domain_levels),
    max_selection = pmax(selection_female, selection_male),
    female_selected = selection_female >= female_pi,
    male_selected   = selection_male >= male_pi,
    any_selected    = female_selected | male_selected,
    female_left = -selection_female,
    male_right = selection_male
  ) |>
  arrange(domain, desc(max_selection), term_plot) |>
  mutate(
    term_plot_wrapped = str_wrap(term_plot, width = 26),
    term_order = factor(term_plot_wrapped, levels = rev(unique(term_plot_wrapped)))
  )

label_df <- plot_compare |>
  distinct(term_order, term_plot_wrapped, any_selected) |>
  mutate(
    term_plot_safe = htmlEscape(term_plot_wrapped),
    term_plot_safe = str_replace_all(term_plot_safe, "\n", "<br>"),
    label_markdown = ifelse(
      any_selected,
      paste0("<span style='color:#FF0000;font-weight:800;'>", term_plot_safe, "</span>"),
      term_plot_safe
    )
  )

label_map <- setNames(label_df$label_markdown, label_df$term_order)

x_pad <- 0.10
x_min <- -max(c(plot_compare$selection_female, female_pi), na.rm = TRUE) - x_pad
x_max <-  max(c(plot_compare$selection_male, male_pi), na.rm = TRUE) + x_pad

# =========================
# Plot
# =========================
p <- ggplot(plot_compare, aes(y = term_order)) +
  
  geom_rect(
    aes(
      xmin = x_min,
      xmax = x_max,
      ymin = -Inf,
      ymax = Inf,
      fill = domain
    ),
    alpha = 0.35,
    inherit.aes = FALSE
  ) +
  
  geom_vline(
    xintercept = 0,
    color = "grey35",
    linewidth = 0.7
  ) +
  
  geom_vline(
    xintercept = -female_pi,
    linetype = "dashed",
    linewidth = 0.9,
    color = "#D95F02"
  ) +
  
  geom_vline(
    xintercept = male_pi,
    linetype = "dashed",
    linewidth = 0.9,
    color = "#1B9E77"
  ) +
  
  # female: not selected -> grey
  geom_segment(
    data = plot_compare |> filter(!female_selected),
    aes(
      x = 0,
      xend = female_left,
      yend = term_order
    ),
    color = "grey80",
    linewidth = 0.9,
    lineend = "round"
  ) +
  
  # female: selected -> full orange
  geom_segment(
    data = plot_compare |> filter(female_selected),
    aes(
      x = 0,
      xend = female_left,
      yend = term_order
    ),
    color = "#D95F02",
    linewidth = 1.00,
    lineend = "round"
  ) +
  
  # male: not selected -> grey
  geom_segment(
    data = plot_compare |> filter(!male_selected),
    aes(
      x = 0,
      xend = male_right,
      yend = term_order
    ),
    color = "grey80",
    linewidth = 0.9,
    lineend = "round"
  ) +
  
  # male: selected -> full green
  geom_segment(
    data = plot_compare |> filter(male_selected),
    aes(
      x = 0,
      xend = male_right,
      yend = term_order
    ),
    color = "#1B9E77",
    linewidth = 1.00,
    lineend = "round"
  ) +
  
  # female points
  geom_point(
    aes(x = female_left, y = term_order),
    color = ifelse(plot_compare$female_selected, "#D95F02", "grey70"),
    size = 3.1
  ) +
  
  # male points
  geom_point(
    aes(x = male_right, y = term_order),
    color = ifelse(plot_compare$male_selected, "#1B9E77", "grey70"),
    size = 3.1
  ) +
  
  facet_grid(
    domain ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  
  scale_fill_manual(
    values = domain_fill,
    guide = "none"
  ) +
  
  scale_y_discrete(labels = label_map) +
  
  scale_x_continuous(
    labels = function(x) abs(x),
    limits = c(x_min, x_max),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  
  coord_cartesian(clip = "off") +
  
  labs(
    x = "Selection proportion",
    y = NULL,
    title = "Subsample lasso stability selection: Female vs Male",
    subtitle = paste0(
      "Variables grouped by domain and ordered by selection stability\n",
      "Female π = ", sprintf("%.2f", female_pi),
      "   |   Male π = ", sprintf("%.2f", male_pi)
    ),
    caption = reference_text
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.35),
    
    axis.text.y = ggtext::element_markdown(
      size = 8.6,
      color = "grey15",
      hjust = 1,
      margin = margin(r = 12)
    ),
    axis.text.y.left = ggtext::element_markdown(
      size = 8.6,
      color = "grey15",
      hjust = 1,
      margin = margin(r = 12)
    ),
    axis.text.x = element_text(size = 10, color = "grey20"),
    axis.title.x = element_text(size = 12, face = "bold"),
    
    legend.position = "none",
    
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(
      size = 11,
      hjust = 0.5,
      color = "grey25",
      lineheight = 1.2
    ),
    plot.caption = element_text(
      size = 9,
      hjust = 0,
      color = "grey25",
      lineheight = 1.15
    ),
    
    strip.placement = "outside",
    strip.background = element_rect(
      fill = "grey88",
      colour = NA
    ),
    strip.text.y.left = element_text(
      angle = 0,
      face = "bold",
      size = 11,
      margin = margin(r = 10, l = 10)
    ),
    
    plot.margin = margin(t = 40, r = 35, b = 45, l = 120)
  )

print(p)

# =========================
# Save files
# =========================
png_file <- file.path(base_out_dir, "butterfly_plot_female_vs_male_by_domain.png")
pdf_file <- file.path(base_out_dir, "butterfly_plot_female_vs_male_by_domain.pdf")

plot_height <- max(10.5, nrow(plot_compare) * 0.24 + 3.2)
plot_width  <- 15.8

ggsave(
  filename = png_file,
  plot = p,
  width = plot_width,
  height = plot_height,
  dpi = 320,
  bg = "white"
)

ggsave(
  filename = pdf_file,
  plot = p,
  width = plot_width,
  height = plot_height,
  device = cairo_pdf,
  bg = "white"
)