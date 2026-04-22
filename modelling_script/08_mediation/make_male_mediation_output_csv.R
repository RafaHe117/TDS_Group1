suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

BASE_DIR <- "/rds/general/project/hda_25-26/live/TDS/anw16/TDS_Group1"
MED_DIR <- file.path(BASE_DIR, "modelling_script", "08_mediation_outdated")

REFIT_SPLIT_DIR <- file.path(MED_DIR, "outputs", "formal_mediation_refit_split")
OUT_DIR <- file.path(MED_DIR, "outputs", "final_mediation_outputs")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

safe_read_csv <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  read_csv(path, show_col_types = FALSE)
}

find_result_files <- function(refit_split_dir, analysis_name) {
  pattern <- paste0("^formal_mediation_results_", analysis_name, "(_chunk[0-9]+)?\\.csv$")
  files <- list.files(
    refit_split_dir,
    pattern = pattern,
    recursive = TRUE,
    full.names = TRUE
  )
  sort(files)
}

combine_analysis_results <- function(refit_split_dir, analysis_name) {
  files <- find_result_files(refit_split_dir, analysis_name)

  if (length(files) == 0) {
    stop("No result files found for analysis: ", analysis_name)
  }

  df <- bind_rows(lapply(files, safe_read_csv)) %>%
    mutate(
      analysis = analysis_name,
      pathway = paste(exposure_term, mediator, sep = " -> ")
    ) %>%
    distinct()

  req_cols <- c(
    "analysis", "subgroup_value", "exposure_var", "exposure_term", "mediator",
    "n_complete",
    "a_estimate", "a_se", "a_p",
    "b_estimate", "b_se", "b_p",
    "total_effect_c", "total_se", "total_p",
    "direct_effect_cprime", "direct_se", "direct_p",
    "indirect_effect_ab", "indirect_direction", "note",
    "indirect_boot_mean", "indirect_ci_low", "indirect_ci_high",
    "direct_boot_mean", "direct_ci_low", "direct_ci_high",
    "total_boot_mean", "total_ci_low", "total_ci_high",
    "pathway"
  )

  missing_cols <- setdiff(req_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in ", analysis_name, ": ",
      paste(missing_cols, collapse = ", ")
    )
  }

  df %>%
    arrange(exposure_term, mediator)
}

male_df <- combine_analysis_results(REFIT_SPLIT_DIR, "male")

if (is.null(male_df) || nrow(male_df) == 0) {
  stop("No combined formal mediation results available for male analysis.")
}

write_csv(male_df, file.path(OUT_DIR, "final_mediation_results_male.csv"))

cat("Male mediation csv saved to:\n", file.path(OUT_DIR, "final_mediation_results_male.csv"), "\n")
