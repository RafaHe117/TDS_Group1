suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
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

project_root <- find_project_root()

med_dir <- file.path(project_root, "analysis", "mediation")
inputs_dir <- file.path(med_dir, "inputs")
dir.create(inputs_dir, recursive = TRUE, showWarnings = FALSE)

main_raw_exposure_path <- file.path(
  project_root,
  "modelling_script", "stability_analysis", "subsample_lasso", "exposure_list.csv"
)

main_selected_path <- file.path(
  project_root,
  "modelling_script", "stability_analysis", "subsample_lasso",
  "selection_proportions_penalized.csv"
)

male_selected_path <- file.path(
  project_root,
  "modelling_script", "stability_analysis", "subsample_lasso_sex", "male",
  "selection_proportions_penalized_male.csv"
)

female_selected_path <- file.path(
  project_root,
  "modelling_script", "stability_analysis", "subsample_lasso_sex", "female",
  "selection_proportions_penalized_female.csv"
)

if (!file.exists(main_raw_exposure_path)) {
  stop("Missing raw exposure universe file: ", main_raw_exposure_path)
}

raw_exposure_df <- read_csv(main_raw_exposure_path, show_col_types = FALSE)

if (!"exposure" %in% names(raw_exposure_df)) {
  stop("Column 'exposure' not found in raw exposure input file.")
}

raw_exposure_df <- raw_exposure_df %>%
  distinct(exposure) %>%
  filter(!is.na(exposure), exposure != "") %>%
  arrange(exposure)

write_csv(
  raw_exposure_df,
  file.path(inputs_dir, "raw_exposure_universe.csv")
)

build_selected_terms <- function(path_in, label, out_dir) {
  if (!file.exists(path_in)) {
    stop("Missing selected-terms file for ", label, ": ", path_in)
  }
  
  df <- read_csv(path_in, show_col_types = FALSE)
  
  required_cols <- c("term", "selection_proportion", "selected")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in ", label, " file: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  out <- df %>%
    filter(selected %in% TRUE) %>%
    distinct(term, .keep_all = TRUE) %>%
    arrange(desc(selection_proportion), term)
  
  write_csv(
    out,
    file.path(out_dir, paste0("selected_terms_", label, ".csv"))
  )
  
  message(label, ": ", nrow(out), " selected terms")
}

build_selected_terms(main_selected_path, "main", inputs_dir)
build_selected_terms(male_selected_path, "male", inputs_dir)
build_selected_terms(female_selected_path, "female", inputs_dir)

message("Done.")
message("Inputs saved to: ", inputs_dir)