# analysis/table1/table1_main_bysex_before_after_run.R
# ------------------------------------------------------------
# PURPOSE
# Create FOUR new main-paper Table 1 tables only:
#   1) Female only, stratified by cvd_incident, BEFORE imputation
#   2) Female only, stratified by cvd_incident, AFTER imputation
#   3) Male only, stratified by cvd_incident, BEFORE imputation
#   4) Male only, stratified by cvd_incident, AFTER imputation
# NOTE
# - It does NOT modify any previous scripts
# - It does NOT overwrite old table1 outputs
# - sex is removed from row variables because tables are split by sex
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(grid)
  library(gridExtra)
})

source("analysis/table1/table1_config.R")
source("analysis/table1/table1_utils.R")

# -----------------------------
# Options
# -----------------------------
RENDER_PNG <- TRUE

# Put all new outputs in a separate folder
OUTDIR_BYSEX <- file.path("analysis", "table1", "output_main_bysex_before_after")
dir.create(OUTDIR_BYSEX, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Read data
# -----------------------------
message("Reading clean dataset ...")
df_clean <- readRDS(PATH_CLEAN)

message("Reading imputed dataset ...")
df_imputed <- readRDS(PATH_IMPUTED)

# -----------------------------
# Main variables for paper
# Remove sex because sex is used to split the tables
# -----------------------------
vars_main_use <- setdiff(vars_main, vars_exclude_always)
vars_main_use <- setdiff(vars_main_use, "sex")

# keep only variables present in BOTH datasets for stability
vars_main_use <- intersect(vars_main_use, names(df_clean))
vars_main_use <- intersect(vars_main_use, names(df_imputed))

if (length(vars_main_use) == 0) {
  stop("No variables available for sex-specific main Table 1 after filtering.", call. = FALSE)
}

# -----------------------------
# Check required columns
# -----------------------------
required_cols <- c("sex", "cvd_incident")
missing_clean <- setdiff(required_cols, names(df_clean))
missing_imp   <- setdiff(required_cols, names(df_imputed))

if (length(missing_clean) > 0) {
  stop("Clean dataset missing required columns: ",
       paste(missing_clean, collapse = ", "), call. = FALSE)
}
if (length(missing_imp) > 0) {
  stop("Imputed dataset missing required columns: ",
       paste(missing_imp, collapse = ", "), call. = FALSE)
}

# -----------------------------
# Helper: sex subset
# -----------------------------
get_sex_subset <- function(data, sex_value, data_name = "dataset") {
  idx <- !is.na(data$sex) & as.character(data$sex) == sex_value
  out <- data[idx, , drop = FALSE]
  
  if (nrow(out) == 0) {
    stop(
      paste0(
        "No rows found for sex = '", sex_value, "' in ", data_name, ". ",
        "Observed non-missing sex values are: ",
        paste(sort(unique(as.character(data$sex[!is.na(data$sex)]))), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  out
}

# -----------------------------
# Helper: write CSV
# -----------------------------
write_tbl <- function(tbl, filename) {
  outpath <- file.path(OUTDIR_BYSEX, filename)
  write.csv(tbl, outpath, row.names = FALSE, na = "")
  message("Wrote CSV: ", outpath)
  outpath
}

# -----------------------------
# Helper: clean for PNG
# -----------------------------
prepare_for_png <- function(df) {
  df2 <- df
  df2[is.na(df2)] <- ""
  
  if ("p_value" %in% names(df2)) {
    df2$p_value <- as.character(df2$p_value)
    df2$p_value[is.na(df2$p_value)] <- ""
    df2$p_value[df2$p_value == "NA"] <- ""
  }
  
  df2
}

# -----------------------------
# Helper: render CSV to PNG
# -----------------------------
render_csv_to_png <- function(csv_path, png_path,
                              width = 15, height = 10, res = 220) {
  df_png <- read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)
  df_png <- prepare_for_png(df_png)
  
  tbl_theme <- ttheme_minimal(
    base_size = 9,
    core = list(
      fg_params = list(cex = 0.8),
      padding = unit(c(2.5, 2.5), "mm")
    ),
    colhead = list(
      fg_params = list(fontface = "bold", cex = 0.85)
    )
  )
  
  g <- tableGrob(df_png, rows = NULL, theme = tbl_theme)
  
  png(filename = png_path, width = width, height = height, units = "in", res = res)
  grid.newpage()
  grid.text(
    label = tools::file_path_sans_ext(basename(csv_path)),
    x = 0.02, y = 0.98,
    just = c("left", "top"),
    gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "Times")
  )
  pushViewport(viewport(x = 0.5, y = 0.46, width = 0.98, height = 0.90))
  grid.draw(g)
  popViewport()
  dev.off()
  
  message("Rendered PNG: ", png_path)
}

# -----------------------------
# Core builder: one table
# -----------------------------
build_one_table <- function(data,
                            sex_value,
                            sex_label_for_file,
                            stage_label,              # "before" or "after"
                            add_missing_column,
                            data_name = "dataset") {
  
  message("Building ", stage_label, " table for sex = ", sex_value, " from ", data_name, " ...")
  
  data_sex <- get_sex_subset(data, sex_value = sex_value, data_name = data_name)
  
  tbl <- make_table1_two_group(
    df = data_sex,
    vars = vars_main_use,
    strata_var = "cvd_incident",
    strata_levels = c(0, 1),
    strata_labels = c(label_outcome_0, label_outcome_1),
    add_missing_overall = add_missing_column,
    collapse_topk = NULL
  )
  
  csv_name <- paste0("table1_main_", sex_label_for_file, "_outcome_", stage_label, ".csv")
  csv_path <- write_tbl(tbl, csv_name)
  
  if (RENDER_PNG) {
    png_name <- paste0("table1_main_", sex_label_for_file, "_outcome_", stage_label, ".png")
    png_path <- file.path(OUTDIR_BYSEX, png_name)
    render_csv_to_png(csv_path, png_path)
  }
  
  invisible(tbl)
}

# -----------------------------
# Run all four
# Rule:
# - BEFORE = clean data + missingness column
# - AFTER  = imputed data + no missingness column
# -----------------------------

# Female - before
build_one_table(
  data = df_clean,
  sex_value = "Female",
  sex_label_for_file = "female",
  stage_label = "before",
  add_missing_column = TRUE,
  data_name = "clean dataset"
)

# Female - after
build_one_table(
  data = df_imputed,
  sex_value = "Female",
  sex_label_for_file = "female",
  stage_label = "after",
  add_missing_column = FALSE,
  data_name = "imputed dataset"
)

# Male - before
build_one_table(
  data = df_clean,
  sex_value = "Male",
  sex_label_for_file = "male",
  stage_label = "before",
  add_missing_column = TRUE,
  data_name = "clean dataset"
)

# Male - after
build_one_table(
  data = df_imputed,
  sex_value = "Male",
  sex_label_for_file = "male",
  stage_label = "after",
  add_missing_column = FALSE,
  data_name = "imputed dataset"
)

message("Done. All new sex-specific main tables saved to: ", OUTDIR_BYSEX)