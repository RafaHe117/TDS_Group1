suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(Matrix)
  library(sharp)
  library(parallelly)
})

find_project_root <- function(start = getwd()) {
  cur <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (
      dir.exists(file.path(cur, "split_imputed_data")) &&
      dir.exists(file.path(cur, "modelling_script")) &&
      dir.exists(file.path(cur, "analysis"))
    ) {
      return(cur)
    }
    parent <- dirname(cur)
    if (parent == cur) {
      stop("Could not locate project root from: ", start)
    }
    cur <- parent
  }
}

get_requested_cores <- function(default = 1L) {
  candidates <- c(
    Sys.getenv("PBS_NCPUS", unset = ""),
    Sys.getenv("PBS_NP", unset = ""),
    Sys.getenv("NCPUS", unset = ""),
    Sys.getenv("MY_NCORES", unset = "")
  )
  candidates <- suppressWarnings(as.integer(candidates))
  candidates <- candidates[!is.na(candidates) & candidates > 0]
  if (length(candidates) == 0) return(as.integer(default))
  candidates[1]
}

coerce_numeric_like <- function(x) {
  if (is.numeric(x)) return(x)
  if (is.logical(x)) return(factor(x))
  if (is.factor(x) || is.character(x)) {
    x_chr <- as.character(x)
    non_missing <- x_chr[!is.na(x_chr) & x_chr != ""]
    if (length(non_missing) > 0) {
      suppress_num <- suppressWarnings(as.numeric(non_missing))
      if (all(!is.na(suppress_num))) {
        return(as.numeric(x_chr))
      }
    }
    return(factor(x))
  }
  x
}

build_predictor_matrix <- function(dat, predictors) {
  x_df <- dat[, predictors, drop = FALSE]
  for (nm in names(x_df)) {
    x_df[[nm]] <- coerce_numeric_like(x_df[[nm]])
  }
  Matrix::sparse.model.matrix(~ . - 1, data = x_df)
}

make_confounder_term_index <- function(colnms, confounders) {
  idx <- rep(FALSE, length(colnms))
  names(idx) <- colnms
  for (cf in confounders) {
    idx <- idx | grepl(paste0("^", cf), colnms)
  }
  idx
}

detect_biomarkers <- function(df, biomarker_file) {
  if (!file.exists(biomarker_file)) {
    stop("Missing biomarker candidate file: ", biomarker_file)
  }
  biom <- read_csv(biomarker_file, show_col_types = FALSE)
  if (!"biomarker" %in% names(biom)) {
    stop("Column 'biomarker' not found in: ", biomarker_file)
  }
  vars <- biom$biomarker
  vars <- vars[vars %in% names(df)]
  vars <- unique(vars[!is.na(vars) & vars != ""])
  sort(vars)
}

run_one_biomarker <- function(df, biomarker, predictors, selected_terms, confounders,
                              k_subsamples, tau, n_cores) {
  vars_needed <- unique(c(biomarker, predictors))
  dat <- df[, vars_needed, drop = FALSE]
  dat <- dat[complete.cases(dat), , drop = FALSE]
  
  if (nrow(dat) < 200) {
    return(NULL)
  }
  
  x_mat <- build_predictor_matrix(dat, predictors)
  y_vec <- dat[[biomarker]]
  
  available_selected_terms <- intersect(selected_terms, colnames(x_mat))
  confounder_idx_full <- make_confounder_term_index(colnames(x_mat), confounders)
  confounder_terms <- colnames(x_mat)[confounder_idx_full]
  
  keep_terms <- unique(c(confounder_terms, available_selected_terms))
  keep_terms <- keep_terms[keep_terms %in% colnames(x_mat)]
  
  if (length(available_selected_terms) == 0) {
    return(NULL)
  }
  
  x_use <- as.matrix(x_mat[, keep_terms, drop = FALSE])
  
  pf <- rep(1, ncol(x_use))
  names(pf) <- colnames(x_use)
  pf[colnames(x_use) %in% confounder_terms] <- 0
  
  stab <- sharp::VariableSelection(
    xdata = x_use,
    ydata = y_vec,
    family = "gaussian",
    implementation = PenalisedRegression,
    penalty.factor = pf,
    alpha = 1,
    standardize = TRUE,
    verbose = FALSE,
    K = k_subsamples,
    resampling = "subsampling",
    tau = tau,
    n_cores = n_cores,
    n_cat = 3
  )
  
  best_par <- sharp::Argmax(stab)
  best_lambda <- best_par[1, 1]
  best_pi <- best_par[1, 2]
  
  sel_prop <- sharp::SelectionProportions(stab)
  
  tibble(
    biomarker = biomarker,
    n_complete = nrow(dat),
    term = names(sel_prop),
    selection_proportion = as.numeric(sel_prop),
    best_pi = as.numeric(best_pi),
    best_lambda = as.numeric(best_lambda),
    selected = as.numeric(sel_prop) >= as.numeric(best_pi),
    term_role = ifelse(names(sel_prop) %in% confounder_terms, "confounder", "selected_exposure_term")
  ) %>%
    arrange(desc(selection_proportion), term)
}

set.seed(2026)

project_root <- find_project_root()

med_dir <- file.path(project_root, "analysis", "mediation")
inputs_dir <- file.path(med_dir, "inputs")
out_dir <- file.path(med_dir, "outputs", "main")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

data_path <- file.path(project_root, "split_imputed_data", "ukb_G1_train_imputed.rds")
raw_exposure_path <- file.path(inputs_dir, "raw_exposure_universe.csv")
selected_terms_path <- file.path(inputs_dir, "selected_terms_main.csv")
biomarker_file <- file.path(inputs_dir, "biomarker_candidates.csv")

if (!file.exists(data_path)) stop("Missing data file: ", data_path)
if (!file.exists(raw_exposure_path)) stop("Missing raw exposure universe file: ", raw_exposure_path)
if (!file.exists(selected_terms_path)) stop("Missing selected terms file: ", selected_terms_path)

df <- readRDS(data_path)

outcome <- "cvd_incident"
confounders <- c("age", "sex", "ethnicity_5cat")

raw_exposure_vars <- read_csv(raw_exposure_path, show_col_types = FALSE)$exposure
raw_exposure_vars <- raw_exposure_vars[raw_exposure_vars %in% names(df)]

selected_terms_df <- read_csv(selected_terms_path, show_col_types = FALSE)
selected_terms <- unique(selected_terms_df$term)

predictors <- unique(c(confounders, raw_exposure_vars))
predictors <- predictors[predictors %in% names(df)]

missing_core <- setdiff(c(outcome, confounders), names(df))
if (length(missing_core) > 0) {
  stop("Missing required columns in data: ", paste(missing_core, collapse = ", "))
}

biomarkers <- detect_biomarkers(df, biomarker_file)

write_csv(
  tibble(biomarker = biomarkers),
  file.path(out_dir, "detected_biomarkers_main.csv")
)

if (length(biomarkers) == 0) {
  stop("No biomarker candidates available in data.")
}

requested_cores <- get_requested_cores(default = 2L)
available_cores <- parallelly::availableCores()
n_cores <- min(requested_cores, available_cores)

k_subsamples <- 100L
tau_subsample <- 0.5

message("Detected biomarkers: ", length(biomarkers))
message("Selected design terms from prior stability stage: ", length(selected_terms))
message("n_cores = ", n_cores)
message("n_cat = 3")

results <- vector("list", length(biomarkers))
names(results) <- biomarkers

for (i in seq_along(biomarkers)) {
  biomarker_i <- biomarkers[i]
  message("[", i, "/", length(biomarkers), "] Running biomarker: ", biomarker_i)
  
  results[[i]] <- run_one_biomarker(
    df = df,
    biomarker = biomarker_i,
    predictors = predictors,
    selected_terms = selected_terms,
    confounders = confounders,
    k_subsamples = k_subsamples,
    tau = tau_subsample,
    n_cores = n_cores
  )
}

final_res <- bind_rows(results)

write_csv(
  final_res,
  file.path(out_dir, "biomarker_stability_main.csv")
)

stable_links <- final_res %>%
  filter(term_role == "selected_exposure_term", selected %in% TRUE) %>%
  arrange(biomarker, desc(selection_proportion), term)

write_csv(
  stable_links,
  file.path(out_dir, "stable_links_main.csv")
)

summary_df <- tibble(
  analysis = "main",
  n_biomarkers_tested = length(biomarkers),
  n_rows_output = nrow(final_res),
  n_stable_links = nrow(stable_links),
  n_selected_terms_input = length(selected_terms),
  n_raw_exposure_vars = length(raw_exposure_vars),
  n_cores = n_cores,
  k_subsamples = k_subsamples,
  tau = tau_subsample,
  n_cat = 3
)

write_csv(
  summary_df,
  file.path(out_dir, "summary_main.csv")
)

message("Done.")
message("Outputs saved to: ", out_dir)