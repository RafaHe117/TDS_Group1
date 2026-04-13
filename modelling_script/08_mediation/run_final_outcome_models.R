suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

BASE_DIR <- "/rds/general/project/hda_25-26/live/TDS/TDS_Group1"

DATA_PATH <- file.path(BASE_DIR, "split_imputed_data/ukb_G1_train_imputed.rds")

MED_DIR <- file.path(BASE_DIR, "analysis", "mediation")
OUT_BASE <- file.path(MED_DIR, "outputs", "final_outcome_models")
dir.create(OUT_BASE, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# helper: clean term -> base variable
# -----------------------------
INPUT_DIR <- file.path(MED_DIR, "inputs")
RAW_EXPOSURE_PATH <- file.path(INPUT_DIR, "raw_exposure_universe.csv")

raw_exposure_universe <- read_csv(RAW_EXPOSURE_PATH, show_col_types = FALSE) %>%
  pull(exposure) %>%
  unique() %>%
  sort()

term_to_base_var <- function(term, exposure_universe, confounders = c("age", "sex", "ethnicity_5cat")) {
  if (is.na(term) || term == "") return(NA_character_)
  
  valid_prefixes <- unique(c(exposure_universe, confounders))
  valid_prefixes <- valid_prefixes[!is.na(valid_prefixes) & valid_prefixes != ""]
  
  # find all prefixes that match the beginning of the term
  hits <- valid_prefixes[startsWith(term, valid_prefixes)]
  
  if (length(hits) > 0) {
    # choose the longest match
    return(hits[which.max(nchar(hits))])
  }
  
  NA_character_
}

# -----------------------------
# helper: choose selected exposures/biomarkers from stable links
# -----------------------------
extract_selected_sets <- function(stable_links_df) {
  req_cols <- c("biomarker", "term", "selected")
  miss <- setdiff(req_cols, names(stable_links_df))
  if (length(miss) > 0) {
    stop("stable_links file missing columns: ", paste(miss, collapse = ", "))
  }
  
  stable_links_df <- stable_links_df %>%
    mutate(
      base_var = vapply(
        term,
        term_to_base_var,
        character(1),
        exposure_universe = raw_exposure_universe
      ),
      term_role = case_when(
        !is.na(base_var) & base_var == term ~ "plain_term",
        TRUE ~ "expanded_term"
      )
    )
  
  selected_rows <- stable_links_df %>%
    filter(selected %in% TRUE)
  
  selected_exposures <- selected_rows %>%
    filter(!is.na(base_var)) %>%
    pull(base_var) %>%
    unique() %>%
    sort()
  
  selected_biomarkers <- selected_rows %>%
    pull(biomarker) %>%
    unique() %>%
    sort()
  
  list(
    selected_exposures = selected_exposures,
    selected_biomarkers = selected_biomarkers
  )
}

# -----------------------------
# helper: fit one final model
# -----------------------------
run_final_model <- function(df, analysis_name, stable_links_path, sex_subset = NULL) {
  out_dir <- file.path(OUT_BASE, analysis_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  stable_links <- read_csv(stable_links_path, show_col_types = FALSE)
  sets <- extract_selected_sets(stable_links)
  
  selected_exposures <- sets$selected_exposures
  selected_biomarkers <- sets$selected_biomarkers
  
  confounders <- c("age", "sex", "ethnicity_5cat")
  outcome <- "cvd_incident"
  
  dat <- df
  
  # sex-stratified subset
  if (!is.null(sex_subset)) {
    if (!("sex" %in% names(dat))) {
      stop("sex variable not found in data.")
    }
    
    # data shows levels Female / Male
    dat <- dat %>% filter(as.character(sex) == sex_subset)
    
    # in sex-stratified models, drop sex from confounders
    confounders <- setdiff(confounders, "sex")
  }
  
  needed_vars <- unique(c(outcome, confounders, selected_exposures, selected_biomarkers))
  needed_vars <- needed_vars[needed_vars %in% names(dat)]
  
  dat_model <- dat[, needed_vars, drop = FALSE]
  dat_model <- dat_model[complete.cases(dat_model), , drop = FALSE]
  
  if (nrow(dat_model) < 100) {
    stop("Too few complete-case rows for ", analysis_name, ": ", nrow(dat_model))
  }
  
  rhs_terms <- unique(c(selected_exposures, selected_biomarkers, confounders))
  rhs_terms <- rhs_terms[rhs_terms %in% names(dat_model)]
  
  form <- as.formula(
    paste(outcome, "~", paste(rhs_terms, collapse = " + "))
  )
  
  fit <- glm(form, data = dat_model, family = binomial())
  
  sm <- summary(fit)$coefficients
  coef_df <- data.frame(
    term = rownames(sm),
    estimate = sm[, 1],
    std_error = sm[, 2],
    statistic = sm[, 3],
    p_value = sm[, 4],
    stringsAsFactors = FALSE
  ) %>%
    as_tibble() %>%
    mutate(
      odds_ratio = exp(estimate),
      ci_low_or = exp(estimate - 1.96 * std_error),
      ci_high_or = exp(estimate + 1.96 * std_error),
      term_group = case_when(
        term == "(Intercept)" ~ "intercept",
        term %in% confounders ~ "confounder",
        term %in% selected_biomarkers ~ "biomarker_main_effect",
        TRUE ~ "exposure_or_expanded_term"
      )
    )
  
  # mark exposure family for expanded terms
  coef_df <- coef_df %>%
    mutate(
      exposure_family = vapply(
        term,
        term_to_base_var,
        character(1),
        exposure_universe = raw_exposure_universe
      ),
      exposure_family = ifelse(
        exposure_family %in% selected_exposures,
        exposure_family,
        NA_character_
      )
    )
  
  summary_df <- tibble(
    analysis = analysis_name,
    sex_subset = ifelse(is.null(sex_subset), "all", sex_subset),
    n_rows_input = nrow(dat),
    n_rows_model = nrow(dat_model),
    n_selected_exposures = length(selected_exposures),
    n_selected_biomarkers = length(selected_biomarkers),
    n_terms_in_model = length(rhs_terms)
  )
  
  selected_exposure_df <- tibble(exposure = selected_exposures)
  selected_biomarker_df <- tibble(biomarker = selected_biomarkers)
  
  write_csv(summary_df, file.path(out_dir, paste0("summary_", analysis_name, ".csv")))
  write_csv(selected_exposure_df, file.path(out_dir, paste0("selected_exposures_", analysis_name, ".csv")))
  write_csv(selected_biomarker_df, file.path(out_dir, paste0("selected_biomarkers_", analysis_name, ".csv")))
  write_csv(coef_df, file.path(out_dir, paste0("final_model_coefficients_", analysis_name, ".csv")))
  
  saveRDS(fit, file.path(out_dir, paste0("final_model_", analysis_name, ".rds")))
  
  cat("Done:", analysis_name, "\n")
  cat("Rows in model:", nrow(dat_model), "\n")
  cat("Selected exposures:", length(selected_exposures), "\n")
  cat("Selected biomarkers:", length(selected_biomarkers), "\n")
  cat("Output dir:", out_dir, "\n\n")
}

# -----------------------------
# main
# -----------------------------
df <- readRDS(DATA_PATH)

run_final_model(
  df = df,
  analysis_name = "main",
  stable_links_path = file.path(MED_DIR, "outputs", "main", "stable_links_main.csv"),
  sex_subset = NULL
)

run_final_model(
  df = df,
  analysis_name = "male",
  stable_links_path = file.path(MED_DIR, "outputs", "male", "stable_links_male.csv"),
  sex_subset = "Male"
)

run_final_model(
  df = df,
  analysis_name = "female",
  stable_links_path = file.path(MED_DIR, "outputs", "female", "stable_links_female.csv"),
  sex_subset = "Female"
)

cat("All final outcome models completed.\n")