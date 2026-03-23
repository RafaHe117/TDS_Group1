# analysis/table1/table1_run.R
# ------------------------------------------------------------
# Runner: generates 6 tables per spec
# Output: CSVs under analysis/table1/output/
# ------------------------------------------------------------

source("analysis/table1/table1_config.R")
source("analysis/table1/table1_utils.R")

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

message("Reading data ...")
clean   <- readRDS(PATH_CLEAN)
imputed <- readRDS(PATH_IMPUTED)

# Remove always-excluded vars from each set
vars_main_use    <- setdiff(vars_main, vars_exclude_always)
vars_appendix_use<- setdiff(vars_appendix, vars_exclude_always)
vars_bio_use     <- setdiff(vars_bio, vars_exclude_always)

# Helper: write CSV with consistent naming
write_tbl <- function(tbl, filename) {
  outpath <- file.path(OUTDIR, filename)
  write.csv(tbl, outpath, row.names = FALSE)
  message("Wrote: ", outpath)
}

# -------------------------
# 1) MAIN (paper): outcome-stratified (before/after)
#    REQUIREMENT: all BEFORE tables include Missing n (%) (Overall)
# -------------------------
message("Main table by outcome (before imputation) ...")
tbl_main_outcome_before <- make_table1_two_group(
  df = clean,
  vars = vars_main_use,
  strata_var = "cvd_incident",
  strata_levels = c(0, 1),
  strata_labels = c(label_outcome_0, label_outcome_1),
  add_missing_overall = TRUE,     # <-- REQUIRED
  collapse_topk = NULL
)
write_tbl(tbl_main_outcome_before, "table1_main_outcome_before.csv")

message("Main table by outcome (after imputation) ...")
tbl_main_outcome_after <- make_table1_two_group(
  df = imputed,
  vars = vars_main_use,
  strata_var = "cvd_incident",
  strata_levels = c(0, 1),
  strata_labels = c(label_outcome_0, label_outcome_1),
  add_missing_overall = FALSE,
  collapse_topk = NULL
)
write_tbl(tbl_main_outcome_after, "table1_main_outcome_after.csv")

# -------------------------
# 2) APPENDIX: outcome-stratified (before includes missingness + TopK+Other)
# -------------------------
message("Appendix table by outcome (before imputation) ...")
tbl_app_outcome_before <- make_table1_two_group(
  df = clean,
  vars = vars_appendix_use,
  strata_var = "cvd_incident",
  strata_levels = c(0, 1),
  strata_labels = c(label_outcome_0, label_outcome_1),
  add_missing_overall = TRUE,                 # <-- REQUIRED
  collapse_topk = vars_appendix_use,          # collapse categoricals in appendix only
  topk = TOPK_LEVELS
)
write_tbl(tbl_app_outcome_before, "table1_appendix_outcome_before_missing.csv")

message("Appendix table by outcome (after imputation) ...")
tbl_app_outcome_after <- make_table1_two_group(
  df = imputed,
  vars = vars_appendix_use,
  strata_var = "cvd_incident",
  strata_levels = c(0, 1),
  strata_labels = c(label_outcome_0, label_outcome_1),
  add_missing_overall = FALSE,
  collapse_topk = vars_appendix_use,
  topk = TOPK_LEVELS
)
write_tbl(tbl_app_outcome_after, "table1_appendix_outcome_after.csv")

# -------------------------
# 3) BIO / POPULATION: sex-stratified (before/after)
#    Note: sex is the stratifier, do NOT include it as a row variable.
#    REQUIREMENT: all BEFORE tables include Missing n (%) (Overall)
# -------------------------
message("Biological table by sex (before imputation) ...")
tbl_bio_sex_before <- make_table1_two_group(
  df = clean,
  vars = vars_bio_use,
  strata_var = "sex",
  strata_levels = c("Female", "Male"),
  strata_labels = c(label_sex_f, label_sex_m),
  add_missing_overall = TRUE,    # <-- REQUIRED
  collapse_topk = NULL
)
write_tbl(tbl_bio_sex_before, "table1_bio_bysex_before.csv")

message("Biological table by sex (after imputation) ...")
tbl_bio_sex_after <- make_table1_two_group(
  df = imputed,
  vars = vars_bio_use,
  strata_var = "sex",
  strata_levels = c("Female", "Male"),
  strata_labels = c(label_sex_f, label_sex_m),
  add_missing_overall = FALSE,
  collapse_topk = NULL
)
write_tbl(tbl_bio_sex_after, "table1_bio_bysex_after.csv")

message("All done.")