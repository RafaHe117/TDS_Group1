suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

BASE_DIR <- "/rds/general/project/hda_25-26/live/TDS/TDS_Group1"
MED_DIR  <- file.path(BASE_DIR, "analysis", "mediation")
IN_DIR   <- file.path(MED_DIR, "inputs")
OUT_DIR  <- file.path(MED_DIR, "outputs")
SHORT_DIR <- file.path(OUT_DIR, "final_shortlists")
dir.create(SHORT_DIR, recursive = TRUE, showWarnings = FALSE)

RAW_EXPOSURE_PATH <- file.path(IN_DIR, "raw_exposure_universe.csv")

raw_exposure_universe <- read_csv(RAW_EXPOSURE_PATH, show_col_types = FALSE) %>%
  pull(exposure) %>%
  unique() %>%
  sort()

term_to_base_var <- function(term,
                             exposure_universe,
                             confounders = c("age", "sex", "ethnicity_5cat")) {
  if (is.na(term) || term == "") return(NA_character_)
  
  valid_prefixes <- unique(c(exposure_universe, confounders))
  valid_prefixes <- valid_prefixes[!is.na(valid_prefixes) & valid_prefixes != ""]
  
  hits <- valid_prefixes[startsWith(term, valid_prefixes)]
  
  if (length(hits) > 0) {
    return(hits[which.max(nchar(hits))])
  }
  
  NA_character_
}

extract_sig_sets <- function(coef_df, p_cutoff = 0.05) {
  coef_df2 <- coef_df %>%
    mutate(
      base_var = vapply(
        term,
        term_to_base_var,
        character(1),
        exposure_universe = raw_exposure_universe
      )
    )
  
  sig_exposures <- coef_df2 %>%
    filter(
      term != "(Intercept)",
      p_value < p_cutoff,
      term_group == "exposure_or_expanded_term",
      !is.na(base_var)
    ) %>%
    pull(base_var) %>%
    unique() %>%
    sort()
  
  sig_biomarkers <- coef_df2 %>%
    filter(
      term != "(Intercept)",
      p_value < p_cutoff,
      term_group == "biomarker_main_effect"
    ) %>%
    pull(term) %>%
    unique() %>%
    sort()
  
  list(
    sig_exposures = sig_exposures,
    sig_biomarkers = sig_biomarkers
  )
}

build_shortlist_one <- function(label, p_cutoff = 0.05) {
  stable_links_path <- file.path(OUT_DIR, label, paste0("stable_links_", label, ".csv"))
  final_coef_path   <- file.path(OUT_DIR, "final_outcome_models", label,
                                 paste0("final_model_coefficients_", label, ".csv"))
  
  stable_links <- read_csv(stable_links_path, show_col_types = FALSE)
  final_coef   <- read_csv(final_coef_path, show_col_types = FALSE)
  
  sig_sets <- extract_sig_sets(final_coef, p_cutoff = p_cutoff)
  
  sig_exposures  <- sig_sets$sig_exposures
  sig_biomarkers <- sig_sets$sig_biomarkers
  
  links2 <- stable_links %>%
    mutate(
      base_exposure = vapply(
        term,
        term_to_base_var,
        character(1),
        exposure_universe = raw_exposure_universe
      ),
      exposure_sig_in_final  = !is.na(base_exposure) & base_exposure %in% sig_exposures,
      biomarker_sig_in_final = biomarker %in% sig_biomarkers
    ) %>%
    filter(selected %in% TRUE)
  
  # 1) strict: both exposure and biomarker significant in final model
  strict_shortlist <- links2 %>%
    filter(exposure_sig_in_final, biomarker_sig_in_final) %>%
    arrange(biomarker, base_exposure, desc(selection_proportion), term) %>%
    distinct(biomarker, base_exposure, .keep_all = TRUE)
  
  # 2) medium: biomarker significant in final model, exposure was selected upstream
  medium_shortlist <- links2 %>%
    filter(biomarker_sig_in_final) %>%
    arrange(biomarker, base_exposure, desc(selection_proportion), term) %>%
    distinct(biomarker, base_exposure, .keep_all = TRUE)
  
  # 3) all selected links: keep everything selected for record
  all_selected_links <- links2 %>%
    arrange(biomarker, base_exposure, desc(selection_proportion), term)
  
  # save
  write_csv(
    strict_shortlist,
    file.path(SHORT_DIR, paste0("strict_shortlist_", label, ".csv"))
  )
  
  write_csv(
    medium_shortlist,
    file.path(SHORT_DIR, paste0("medium_shortlist_", label, ".csv"))
  )
  
  write_csv(
    all_selected_links,
    file.path(SHORT_DIR, paste0("all_selected_links_", label, ".csv"))
  )
  
  summary_tbl <- tibble(
    analysis = label,
    p_cutoff = p_cutoff,
    n_sig_exposures_final = length(sig_exposures),
    n_sig_biomarkers_final = length(sig_biomarkers),
    n_selected_links = nrow(all_selected_links),
    n_medium_shortlist = nrow(medium_shortlist),
    n_strict_shortlist = nrow(strict_shortlist)
  )
  
  write_csv(
    summary_tbl,
    file.path(SHORT_DIR, paste0("summary_shortlist_", label, ".csv"))
  )
  
  cat("Done:", label, "\n")
  cat("  sig exposures in final model:", length(sig_exposures), "\n")
  cat("  sig biomarkers in final model:", length(sig_biomarkers), "\n")
  cat("  all selected links:", nrow(all_selected_links), "\n")
  cat("  medium shortlist:", nrow(medium_shortlist), "\n")
  cat("  strict shortlist:", nrow(strict_shortlist), "\n\n")
}

build_shortlist_one("main", p_cutoff = 0.05)
build_shortlist_one("male", p_cutoff = 0.05)
build_shortlist_one("female", p_cutoff = 0.05)

cat("All shortlist files saved to:\n", SHORT_DIR, "\n")
