suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

BASE_DIR <- "/rds/general/project/hda_25-26/live/TDS/anw16/TDS_Group1"
MED_DIR  <- file.path(BASE_DIR, "modelling_script", "08_mediation_outdated")
INPUT_DIR <- file.path(MED_DIR, "inputs")
OUT_DIR   <- file.path(MED_DIR, "outputs")
SHORT_DIR <- file.path(OUT_DIR, "final_shortlists")

dir.create(INPUT_DIR, recursive = TRUE, showWarnings = FALSE)

split_into_two_chunks <- function(df) {
  n <- nrow(df)
  if (n == 0) {
    return(list(df[0, , drop = FALSE], df[0, , drop = FALSE]))
  }
  cut_point <- ceiling(n / 2)
  list(
    df[seq_len(n) <= cut_point, , drop = FALSE],
    df[seq_len(n) >  cut_point, , drop = FALSE]
  )
}

build_pair_df <- function(shortlist_df, analysis_label, subgroup_label) {
  req_cols <- c("base_exposure", "term", "biomarker")
  miss <- setdiff(req_cols, names(shortlist_df))
  if (length(miss) > 0) {
    stop(
      "Shortlist file for ", analysis_label,
      " is missing required columns: ",
      paste(miss, collapse = ", ")
    )
  }

  shortlist_df %>%
    transmute(
      analysis = analysis_label,
      subgroup_label = subgroup_label,
      exposure_var = base_exposure,
      exposure_term = term,
      mediator = biomarker
    ) %>%
    filter(
      !is.na(exposure_var), exposure_var != "",
      !is.na(exposure_term), exposure_term != "",
      !is.na(mediator), mediator != ""
    ) %>%
    distinct() %>%
    arrange(exposure_var, exposure_term, mediator)
}

strict_male_path <- file.path(SHORT_DIR, "strict_shortlist_male.csv")

if (!file.exists(strict_male_path)) {
  stop("Missing file: ", strict_male_path)
}

strict_male <- read_csv(strict_male_path, show_col_types = FALSE)
pairs_male <- build_pair_df(strict_male, "male", "male")

if (nrow(pairs_male) == 0) {
  stop("No male mediation pairs available after strict shortlist filtering.")
}

male_chunks <- split_into_two_chunks(pairs_male)

write_csv(male_chunks[[1]], file.path(INPUT_DIR, "pairs_male_chunk1.csv"))
write_csv(male_chunks[[2]], file.path(INPUT_DIR, "pairs_male_chunk2.csv"))

summary_tbl <- tibble(
  file = c("pairs_male_chunk1.csv", "pairs_male_chunk2.csv"),
  n_rows = c(nrow(male_chunks[[1]]), nrow(male_chunks[[2]]))
)

write_csv(summary_tbl, file.path(INPUT_DIR, "formal_mediation_input_summary_male.csv"))

cat("Done.\n")
cat("Male pair files written to:\n")
cat(" - ", file.path(INPUT_DIR, "pairs_male_chunk1.csv"), "\n", sep = "")
cat(" - ", file.path(INPUT_DIR, "pairs_male_chunk2.csv"), "\n", sep = "")
print(summary_tbl)
