# analysis/table1/table1_icd_run.R
# ------------------------------------------------------------
# PURPOSE
# Standalone ICD Table 1-style tables
# - outcome stratified by cvd_incident
# - before imputation: with Missing n (%) (Overall)
# - after imputation: without missing column
# - ICD flags forced to categorical Yes/No
# - final PNG rendering hides NA in p_value cells
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(gridExtra)
  library(grid)
})

source("analysis/table1/table1_config.R")
source("analysis/table1/table1_utils.R")

# -----------------------------
# Paths
# -----------------------------
PATH_CLEAN   <- "ukb_G1_cleaned.rds"
PATH_IMPUTED <- "ukb_G1_imputed_final.rds"

OUTDIR_ICD <- file.path("analysis", "table1", "icd_cat_output")
dir.create(OUTDIR_ICD, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# ICD variables
# -----------------------------
vars_icd <- c(
  "ICD_Diabetes",
  "ICD_Lipidemia",
  "ICD_Hypertension",
  "ICD_AF",
  "ICD_CKD",
  "ICD_Mental_Health",
  "ICD_Migraine",
  "ICD_Atopy",
  "ICD_Autoimmune"
)

# -----------------------------
# Helpers
# -----------------------------
keep_present <- function(df, vars) {
  intersect(vars, names(df))
}

coerce_icd_yesno <- function(df, vars) {
  vars <- intersect(vars, names(df))
  
  for (v in vars) {
    x <- df[[v]]
    
    if (is.logical(x)) {
      df[[v]] <- factor(x, levels = c(FALSE, TRUE), labels = c("No", "Yes"))
      
    } else if (is.numeric(x) || is.integer(x)) {
      x2 <- ifelse(
        is.na(x), NA,
        ifelse(x == 1, "Yes",
               ifelse(x == 0, "No", as.character(x)))
      )
      df[[v]] <- factor(x2)
      
    } else {
      x_chr <- trimws(as.character(x))
      x_std <- ifelse(
        is.na(x_chr), NA,
        ifelse(x_chr %in% c("1", "Yes", "YES", "yes", "Y", "True", "TRUE", "true"),
               "Yes",
               ifelse(x_chr %in% c("0", "No", "NO", "no", "N", "False", "FALSE", "false"),
                      "No",
                      x_chr))
      )
      df[[v]] <- factor(x_std)
    }
  }
  
  df
}

write_tbl <- function(tbl, filename) {
  outpath <- file.path(OUTDIR_ICD, filename)
  write.csv(tbl, outpath, row.names = FALSE, na = "")
  message("Wrote CSV: ", outpath)
}

# For PNG display only:
# hide NA / "NA" in p_value column
prepare_for_png <- function(df) {
  df2 <- df
  
  if ("p_value" %in% names(df2)) {
    df2$p_value <- as.character(df2$p_value)
    df2$p_value[is.na(df2$p_value)] <- ""
    df2$p_value[df2$p_value == "NA"] <- ""
  }
  
  df2
}

render_csv_to_png <- function(csv_path, png_path,
                              width = 16, height = 8, res = 220) {
  df <- read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)
  df <- prepare_for_png(df)
  
  tbl_theme <- ttheme_minimal(
    core = list(fg_params = list(cex = 0.8)),
    colhead = list(fg_params = list(fontface = "bold", cex = 0.85))
  )
  
  g <- tableGrob(df, rows = NULL, theme = tbl_theme)
  
  png(filename = png_path, width = width, height = height, units = "in", res = res)
  grid.newpage()
  pushViewport(viewport(y = 0.48, height = 0.96))
  grid.draw(g)
  dev.off()
  
  message("Rendered PNG: ", png_path)
}

# -----------------------------
# Read data
# -----------------------------
message("Reading clean data ...")
clean <- readRDS(PATH_CLEAN)

message("Reading imputed data ...")
imputed <- readRDS(PATH_IMPUTED)

vars_icd_clean   <- keep_present(clean, vars_icd)
vars_icd_imputed <- keep_present(imputed, vars_icd)

if (length(vars_icd_clean) == 0) {
  stop("No ICD variables found in clean dataset.", call. = FALSE)
}
if (length(vars_icd_imputed) == 0) {
  stop("No ICD variables found in imputed dataset.", call. = FALSE)
}

clean   <- coerce_icd_yesno(clean, vars_icd_clean)
imputed <- coerce_icd_yesno(imputed, vars_icd_imputed)

# -----------------------------
# Before imputation table
# -----------------------------
message("Generating ICD before-imputation table ...")
tbl_icd_before <- make_table1_two_group(
  df = clean,
  vars = vars_icd_clean,
  strata_var = "cvd_incident",
  strata_levels = c(0, 1),
  strata_labels = c(label_outcome_0, label_outcome_1),
  add_missing_overall = TRUE,
  collapse_topk = NULL
)

write_tbl(tbl_icd_before, "table1_icd_outcome_before_missing.csv")

# -----------------------------
# After imputation table
# -----------------------------
message("Generating ICD after-imputation table ...")
tbl_icd_after <- make_table1_two_group(
  df = imputed,
  vars = vars_icd_imputed,
  strata_var = "cvd_incident",
  strata_levels = c(0, 1),
  strata_labels = c(label_outcome_0, label_outcome_1),
  add_missing_overall = FALSE,
  collapse_topk = NULL
)

write_tbl(tbl_icd_after, "table1_icd_outcome_after.csv")

# -----------------------------
# Render PNGs
# -----------------------------
render_csv_to_png(
  csv_path = file.path(OUTDIR_ICD, "table1_icd_outcome_before_missing.csv"),
  png_path = file.path(OUTDIR_ICD, "table1_icd_outcome_before_missing.png")
)

render_csv_to_png(
  csv_path = file.path(OUTDIR_ICD, "table1_icd_outcome_after.csv"),
  png_path = file.path(OUTDIR_ICD, "table1_icd_outcome_after.png")
)

message("Done. ICD tables saved to: ", OUTDIR_ICD)