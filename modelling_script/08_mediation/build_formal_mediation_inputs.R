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

# -----------------------------------
# helper: split into 2 chunks
# -----------------------------------
split_into_two_chunks <- function(df) {
  n <- nrow(df)
  
  if (n == 0) {
    return(list(df[0, , drop = FALSE], df[0, , drop = FALSE]))
  }
  
  idx <- seq_len(n)
  cut_point <- ceiling(n / 2)
  
  list(
    df[idx <= cut_point, , drop = FALSE],
    df[idx >  cut_point, , drop = FALSE]
  )
}

# -----------------------------------
# helper: build pair file from shortlist
# -----------------------------------
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

# -----------------------------------
# paths
# -----------------------------------
strict_main_path   <- file.path(SHORT_DIR, "strict_shortlist_main.csv")
strict_female_path <- file.path(SHORT_DIR, "strict_shortlist_female.csv")

if (!file.exists(strict_main_path)) {
  stop("Missing file: ", strict_main_path)
}
if (!file.exists(strict_female_path)) {
  stop("Missing file: ", strict_female_path)
}

strict_main   <- read_csv(strict_main_path, show_col_types = FALSE)
strict_female <- read_csv(strict_female_path, show_col_types = FALSE)

pairs_main   <- build_pair_df(strict_main, "main", "all")
pairs_female <- build_pair_df(strict_female, "female", "female")

if (nrow(pairs_main) == 0) {
  stop("No main mediation pairs available after strict shortlist filtering.")
}

if (nrow(pairs_female) == 0) {
  stop("No female mediation pairs available after strict shortlist filtering.")
}

main_chunks   <- split_into_two_chunks(pairs_main)
female_chunks <- split_into_two_chunks(pairs_female)

# -----------------------------------
# write config
# -----------------------------------
config_tbl <- tribble(
  ~key, ~value,
  "outcome", "cvd_incident",
  "sex_var", "sex",
  "confounders", "age|sex|ethnicity_5cat",
  "n_boot", "1000",
  "min_complete_n", "100"
)

write_csv(config_tbl, file.path(INPUT_DIR, "formal_mediation_config.csv"))

# -----------------------------------
# write pair csvs
# -----------------------------------
write_csv(main_chunks[[1]],   file.path(INPUT_DIR, "pairs_main_chunk1.csv"))
write_csv(main_chunks[[2]],   file.path(INPUT_DIR, "pairs_main_chunk2.csv"))
write_csv(female_chunks[[1]], file.path(INPUT_DIR, "pairs_female_chunk1.csv"))
write_csv(female_chunks[[2]], file.path(INPUT_DIR, "pairs_female_chunk2.csv"))

# optional summary
summary_tbl <- tibble(
  file = c(
    "pairs_main_chunk1.csv",
    "pairs_main_chunk2.csv",
    "pairs_female_chunk1.csv",
    "pairs_female_chunk2.csv"
  ),
  n_rows = c(
    nrow(main_chunks[[1]]),
    nrow(main_chunks[[2]]),
    nrow(female_chunks[[1]]),
    nrow(female_chunks[[2]])
  )
)

write_csv(summary_tbl, file.path(INPUT_DIR, "formal_mediation_input_summary.csv"))

cat("Done.\n")
cat("Config written to:\n", file.path(INPUT_DIR, "formal_mediation_config.csv"), "\n")
cat("Pair files written to:\n")
cat(" - ", file.path(INPUT_DIR, "pairs_main_chunk1.csv"), "\n", sep = "")
cat(" - ", file.path(INPUT_DIR, "pairs_main_chunk2.csv"), "\n", sep = "")
cat(" - ", file.path(INPUT_DIR, "pairs_female_chunk1.csv"), "\n", sep = "")
cat(" - ", file.path(INPUT_DIR, "pairs_female_chunk2.csv"), "\n", sep = "")
cat("Counts:\n")
print(summary_tbl)