library(dplyr)

in_dir <- "/rds/general/project/hda_25-26/live/TDS/TDS_Group1"

out_dir <- file.path(in_dir, "sanity_check_output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

train_df <- readRDS(file.path(in_dir, "split_imputed_data", "ukb_G1_train_imputed.rds"))
test_df  <- readRDS(file.path(in_dir, "split_imputed_data", "ukb_G1_test_imputed.rds"))

outcome <- "cvd_incident"

# Missingness check
missing_summary <- function(df, dataset_name) {
  data.frame(
    dataset = dataset_name,
    variable = names(df),
    n = nrow(df),
    missing_n = sapply(df, function(x) sum(is.na(x))),
    missing_pct = sapply(df, function(x) mean(is.na(x)) * 100),
    stringsAsFactors = FALSE
  ) |>
    arrange(desc(missing_pct), variable)
}

missing_all <- bind_rows(
  missing_summary(train_df, "train"),
  missing_summary(test_df, "test")
)

write.csv(
  missing_all,
  file.path(out_dir, "missingness_train_test_all_variables.csv"),
  row.names = FALSE
)

# Outcome proportion check
outcome_summary <- function(df, dataset_name, outcome) {
  y <- df[[outcome]]
  
  data.frame(
    dataset = dataset_name,
    outcome = outcome,
    n = length(y),
    missing_n = sum(is.na(y)),
    non_missing_n = sum(!is.na(y)),
    cases_n = sum(y == 1, na.rm = TRUE),
    controls_n = sum(y == 0, na.rm = TRUE),
    case_proportion = mean(y == 1, na.rm = TRUE),
    control_proportion = mean(y == 0, na.rm = TRUE),
    case_pct = round(mean(y == 1, na.rm = TRUE) * 100, 2),
    control_pct = round(mean(y == 0, na.rm = TRUE) * 100, 2),
    stringsAsFactors = FALSE
  )
}

outcome_prop <- bind_rows(
  outcome_summary(train_df, "train", outcome),
  outcome_summary(test_df, "test", outcome)
)

write.csv(
  outcome_prop,
  file.path(out_dir, "cvd_incident_proportion_train_test.csv"),
  row.names = FALSE
)

# Train variable type and counts
var_info <- data.frame(
  variable = names(train_df),
  class = sapply(train_df, function(x) class(x)[1]),
  n_total = nrow(train_df),
  n_non_missing = sapply(train_df, function(x) sum(!is.na(x))),
  n_missing = sapply(train_df, function(x) sum(is.na(x))),
  stringsAsFactors = FALSE
)

write.csv(
  var_info,
  file.path(out_dir, "variable_type_and_counts_train.csv"),
  row.names = FALSE
)

# ukb_G1_cleaned.rds variable type and counts
cleaned_df <- readRDS(file.path(in_dir, "ukb_G1_cleaned.rds"))

cleaned_var_info <- data.frame(
  variable = names(cleaned_df),
  class = sapply(cleaned_df, function(x) class(x)[1]),
  n_total = nrow(cleaned_df),
  n_non_missing = sapply(cleaned_df, function(x) sum(!is.na(x))),
  n_missing = sapply(cleaned_df, function(x) sum(is.na(x))),
  stringsAsFactors = FALSE
)

write.csv(
  cleaned_var_info,
  file.path(out_dir, "variable_type_and_counts_ukb_G1_cleaned.csv"),
  row.names = FALSE
)

# RDS variable names and counts summary
rds_info <- data.frame(
  rds_file = c(
    "ukb_G1_raw.rds",
    "ukb_G1_preprocessed.rds",
    "ukb_G1_imputed.rds",
    "ukb_G1_cleaned.rds"
  ),
  rds_path = c(
    file.path(in_dir, "ukb_G1_raw.rds"),
    file.path(in_dir, "ukb_G1_preprocessed.rds"),
    file.path(in_dir, "imputation", "ukb_G1_imputed.rds"),
    file.path(in_dir, "ukb_G1_cleaned.rds")
  ),
  stringsAsFactors = FALSE
)

rds_var_df <- bind_rows(lapply(seq_len(nrow(rds_info)), function(i) {
  obj <- readRDS(rds_info$rds_path[i])
  
  data.frame(
    rds_file = rds_info$rds_file[i],
    variable_name = names(obj),
    variable_type = sapply(obj, function(x) class(x)[1]),
    stringsAsFactors = FALSE
  )
}))

rds_counts <- rds_var_df |>
  group_by(rds_file) |>
  summarise(variable_count = n(), .groups = "drop")

rds_var_summary <- rds_var_df |>
  left_join(rds_counts, by = "rds_file") |>
  select(rds_file, variable_count, variable_name, variable_type) |>
  arrange(rds_file, variable_name)

write.csv(
  rds_var_summary,
  file.path(out_dir, "rds_variable_names_summary.csv"),
  row.names = FALSE
)

cat("Sanity check outputs saved to:\n")
cat(out_dir, "\n\n")

cat("Generated files:\n")
cat("- missingness_train_test_all_variables.csv\n")
cat("- cvd_incident_proportion_train_test.csv\n")
cat("- variable_type_and_counts_train.csv\n")
cat("- variable_type_and_counts_ukb_G1_cleaned.csv\n")
cat("- rds_variable_names_summary.csv\n\n")

print(outcome_prop)