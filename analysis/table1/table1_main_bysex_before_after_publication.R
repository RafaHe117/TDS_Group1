# analysis/table1/table1_main_bysex_before_after_publication.R
# ------------------------------------------------------------
# PURPOSE
# Publication-ready sex-specific Table 1 (main paper)
# Generates FOUR tables:
#   1) Female only, outcome-stratified, BEFORE imputation
#   2) Female only, outcome-stratified, AFTER imputation
#   3) Male only, outcome-stratified, BEFORE imputation
#   4) Male only, outcome-stratified, AFTER imputation
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

OUTDIR_BYSEX <- file.path("analysis", "table1", "output_main_bysex_before_after_publication")
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
# Remove sex because sex is used to split tables
# -----------------------------
vars_main_use <- setdiff(vars_main, vars_exclude_always)
vars_main_use <- setdiff(vars_main_use, "sex")
vars_main_use <- intersect(vars_main_use, names(df_clean))
vars_main_use <- intersect(vars_main_use, names(df_imputed))

if (length(vars_main_use) == 0) {
  stop("No variables available for publication table after filtering.", call. = FALSE)
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

# ============================================================
# Publication label dictionaries
# ============================================================

# Variable labels for main table
var_label_map <- c(
  age                  = "Age, years",
  ethnicity_5cat       = "Ethnicity",
  imd_quintile         = "Index of Multiple Deprivation quintile",
  urban_rural_2cat     = "Urban/rural residence",
  education_2          = "Educational attainment",
  employment_2cat      = "Employment status",
  living_with_partner  = "Living with partner",
  work_hours_week_clean= "Working hours per week",
  bmi                  = "Body mass index, kg/m²",
  sbp_mean             = "Systolic blood pressure, mmHg",
  dbp_mean             = "Diastolic blood pressure, mmHg",
  diabetes_bin         = "Diabetes",
  famhx_cvd            = "Family history of cardiovascular disease",
  birth_weight_clean   = "Birth weight, kg",
  smoking_3cat         = "Smoking status",
  alcohol_combined     = "Alcohol intake"
)

# Level labels for categorical variables
level_label_map <- list(
  ethnicity_5cat = c(
    "White" = "White",
    "Mixed" = "Mixed",
    "Asian" = "Asian or Asian British",
    "Black" = "Black or Black British",
    "Chinese" = "Chinese",
    "Other" = "Other ethnic group"
  ),
  
  imd_quintile = c(
    "1" = "1 (most deprived)",
    "2" = "2",
    "3" = "3",
    "4" = "4",
    "5" = "5 (least deprived)",
    "Q1" = "1 (most deprived)",
    "Q2" = "2",
    "Q3" = "3",
    "Q4" = "4",
    "Q5" = "5 (least deprived)"
  ),
  
  urban_rural_2cat = c(
    "Urban" = "Urban",
    "Rural" = "Rural"
  ),
  
  education_2 = c(
    "College/University degree" = "College/university degree",
    "No College/University degree" = "No college/university degree",
    "Degree" = "College/university degree",
    "No degree" = "No college/university degree"
  ),
  
  employment_2cat = c(
    "Employed" = "Employed",
    "Not employed" = "Not employed"
  ),
  
  living_with_partner = c(
    "No" = "No",
    "Yes" = "Yes",
    "FALSE" = "No",
    "TRUE" = "Yes"
  ),
  
  diabetes_bin = c(
    "No" = "No",
    "Yes" = "Yes",
    "0" = "No",
    "1" = "Yes",
    "FALSE" = "No",
    "TRUE" = "Yes"
  ),
  
  famhx_cvd = c(
    "No" = "No",
    "Yes" = "Yes",
    "0" = "No",
    "1" = "Yes",
    "FALSE" = "No",
    "TRUE" = "Yes"
  ),
  
  smoking_3cat = c(
    "Never" = "Never",
    "Previous" = "Former",
    "Former" = "Former",
    "Current" = "Current"
  ),
  
  alcohol_combined = c(
    "Never" = "Never",
    "Special occasions only" = "Special occasions only",
    "1-3/month" = "1-3 times/month",
    "1-2/week" = "1-2 times/week",
    "3-4/week" = "3-4 times/week",
    "Daily/almost daily" = "Daily or almost daily"
  )
)

# Optional variable order in final table
final_var_order <- c(
  "Age, years",
  "Ethnicity",
  "Index of Multiple Deprivation quintile",
  "Urban/rural residence",
  "Educational attainment",
  "Employment status",
  "Living with partner",
  "Working hours per week",
  "Body mass index, kg/m²",
  "Systolic blood pressure, mmHg",
  "Diastolic blood pressure, mmHg",
  "Diabetes",
  "Family history of cardiovascular disease",
  "Birth weight, kg",
  "Smoking status",
  "Alcohol intake"
)

# ============================================================
# Helper functions
# ============================================================

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

write_tbl <- function(tbl, filename) {
  outpath <- file.path(OUTDIR_BYSEX, filename)
  write.csv(tbl, outpath, row.names = FALSE, na = "")
  message("Wrote CSV: ", outpath)
  outpath
}

prepare_for_png <- function(df) {
  df2 <- df
  df2[is.na(df2)] <- ""
  
  if ("P value" %in% names(df2)) {
    df2[["P value"]] <- as.character(df2[["P value"]])
    df2[["P value"]][is.na(df2[["P value"]])] <- ""
    df2[["P value"]][df2[["P value"]] == "NA"] <- ""
  }
  
  df2
}

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

# ------------------------------------------------------------
# Relabel rows to publication style
# ------------------------------------------------------------
relabel_table_publication <- function(tbl, add_missing_column = FALSE) {
  tbl2 <- tbl
  
  # Keep original variable id for mapping row levels
  tbl2$variable_raw <- tbl2$variable
  
  # Replace variable labels
  variable_rows <- !is.na(tbl2$variable_raw) & tbl2$variable_raw != ""
  tbl2$variable[variable_rows] <- ifelse(
    tbl2$variable_raw[variable_rows] %in% names(var_label_map),
    unname(var_label_map[tbl2$variable_raw[variable_rows]]),
    tbl2$variable_raw[variable_rows]
  )
  
  # Relabel categorical levels according to original variable
  current_var <- NA_character_
  for (i in seq_len(nrow(tbl2))) {
    if (!is.na(tbl2$variable_raw[i]) && tbl2$variable_raw[i] != "") {
      current_var <- tbl2$variable_raw[i]
    }
    
    if (!is.na(tbl2$level[i]) && tbl2$level[i] != "" && !is.na(current_var)) {
      lvl_trim <- trimws(tbl2$level[i])
      
      if (current_var %in% names(level_label_map)) {
        mapper <- level_label_map[[current_var]]
        if (lvl_trim %in% names(mapper)) {
          tbl2$level[i] <- paste0("  ", unname(mapper[lvl_trim]))
        } else {
          tbl2$level[i] <- paste0("  ", lvl_trim)
        }
      } else {
        tbl2$level[i] <- paste0("  ", lvl_trim)
      }
    }
  }
  
  # Rename columns to publication style
  names(tbl2)[names(tbl2) == "variable"] <- "Characteristic"
  names(tbl2)[names(tbl2) == "level"]    <- ""
  names(tbl2)[names(tbl2) == "p_value"]  <- "P value"
  
  if (add_missing_column) {
    names(tbl2)[names(tbl2) == "Missing n (%) (Overall)"] <- "Missing n (%)"
  }
  
  # Drop helper column
  tbl2$variable_raw <- NULL
  
  # Reorder rows according to publication order
  if ("Characteristic" %in% names(tbl2)) {
    row_group <- tbl2$Characteristic
    row_group[row_group == ""] <- NA
    row_group <- zoo::na.locf(row_group, na.rm = FALSE)
    
    ord_match <- match(row_group, final_var_order)
    ord_match[is.na(ord_match)] <- 999
    
    tbl2$.ord1 <- ord_match
    tbl2$.ord2 <- seq_len(nrow(tbl2))
    tbl2 <- tbl2[order(tbl2$.ord1, tbl2$.ord2), , drop = FALSE]
    tbl2$.ord1 <- NULL
    tbl2$.ord2 <- NULL
  }
  
  rownames(tbl2) <- NULL
  tbl2
}

# ------------------------------------------------------------
# Build one publication table
# ------------------------------------------------------------
build_one_table <- function(data,
                            sex_value,
                            sex_label_for_file,
                            stage_label,
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
  
  tbl_pub <- relabel_table_publication(tbl, add_missing_column = add_missing_column)
  
  csv_name <- paste0("table1_main_", sex_label_for_file, "_outcome_", stage_label, "_publication.csv")
  csv_path <- write_tbl(tbl_pub, csv_name)
  
  if (RENDER_PNG) {
    png_name <- paste0("table1_main_", sex_label_for_file, "_outcome_", stage_label, "_publication.png")
    png_path <- file.path(OUTDIR_BYSEX, png_name)
    render_csv_to_png(csv_path, png_path)
  }
  
  invisible(tbl_pub)
}

# ============================================================
# Run all four
# ============================================================

# Need zoo only for na.locf in row ordering
if (!requireNamespace("zoo", quietly = TRUE)) {
  stop("Package 'zoo' is required for publication row ordering. Please install it first.", call. = FALSE)
}

build_one_table(
  data = df_clean,
  sex_value = "Female",
  sex_label_for_file = "female",
  stage_label = "before",
  add_missing_column = TRUE,
  data_name = "clean dataset"
)

build_one_table(
  data = df_imputed,
  sex_value = "Female",
  sex_label_for_file = "female",
  stage_label = "after",
  add_missing_column = FALSE,
  data_name = "imputed dataset"
)

build_one_table(
  data = df_clean,
  sex_value = "Male",
  sex_label_for_file = "male",
  stage_label = "before",
  add_missing_column = TRUE,
  data_name = "clean dataset"
)

build_one_table(
  data = df_imputed,
  sex_value = "Male",
  sex_label_for_file = "male",
  stage_label = "after",
  add_missing_column = FALSE,
  data_name = "imputed dataset"
)

message("Done. Publication-ready sex-specific main tables saved to: ", OUTDIR_BYSEX)