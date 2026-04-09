suppressPackageStartupMessages(library(miceRanger))

setwd("/rds/general/project/hda_25-26/live/TDS/TDS_Group1")
in_file  <- "ukb_G1_cleaned.rds"
out_dir  <- "imputation"
out_file <- file.path(out_dir, "ukb_G1_imputed.rds")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(20260225)

ukb <- readRDS(in_file)
ukb <- as.data.frame(ukb)

# remove derived vars if present (we regenerate later)
if ("tv_group" %in% names(ukb)) ukb$tv_group <- NULL
if ("total_met_group" %in% names(ukb)) ukb$total_met_group <- NULL

# convert ALL ordered factors -> plain factors (important for miceRanger stability)
is_ord <- sapply(ukb, is.ordered)
if (any(is_ord)) {
  ord_names <- names(ukb)[is_ord]
  for (v in ord_names) ukb[[v]] <- factor(ukb[[v]], ordered = FALSE)
  cat("Converted ordered factors:", paste(ord_names, collapse = ", "), "\n")
  flush.console()
}

get_cores <- function() {
  for (k in c("MY_NCORES", "NCPUS", "PBS_NCPUS", "OMP_NUM_THREADS")) {
    x <- Sys.getenv(k, "")
    if (nzchar(x)) return(as.integer(x))
  }
  1L
}
n_cores <- get_cores()
if (is.na(n_cores) || n_cores < 1) n_cores <- 1L
options(ranger.num.threads = n_cores)

cat(
  "MY_NCORES =", Sys.getenv("MY_NCORES"),
  "| NCPUS =", Sys.getenv("NCPUS"),
  "| PBS_NCPUS =", Sys.getenv("PBS_NCPUS"),
  "| OMP_NUM_THREADS =", Sys.getenv("OMP_NUM_THREADS"),
  "| using n_cores =", n_cores, "\n"
)
flush.console()

imp_vars <- c(
  "smoking_3cat","crp","diabetes_bin","basophil_count",
  "oily_fish_3cat","lipoprotein_a",
  "total_met_min_wk","sf_score","tv_hours_day",
  "salt_3cat","eosinophil_count","sleep_quality_score","mh_satis_mean_score",
  "creatinine","coffee_intake","tea_intake","fruit_intake_fresh",
  "alanine_aminotransferase","lymphocyte_count","total_triglyceride",
  "monocyte_count","redwine_group","greenspace_pct_1000m",
  "total_bilirubin","direct_bilirubin","blood_vitamin_d",
  "aspartate_aminotransferase","imd_final","imd_source","imd_quintile",
  "neutrophil_count","alkaline_phosphatase","wbc_count","no2_2010",
  "urea","hdl_cholesterol","igf_1","urate","employment_2cat",
  "work_hours_week_clean","platelet_count","glucose","ldl_cholesterol",
  "birth_weight_clean","living_with_partner"
)
imp_vars <- intersect(imp_vars, names(ukb))
if (length(imp_vars) == 0) stop("None of imp_vars exist in dataset.")

# force key categoricals to plain factor (not ordered)
cat_vars <- intersect(
  c("smoking_3cat","diabetes_bin","oily_fish_3cat","salt_3cat","redwine_group",
    "imd_quintile","employment_2cat","imd_source","living_with_partner"),
  names(ukb)
)
for (v in cat_vars) ukb[[v]] <- factor(ukb[[v]], ordered = FALSE)

# -------------------------
# DROP LEAKAGE / ID / FUTURE-INFO columns from predictors
# (these will NOT be used as predictors; they can still exist in ukb output)
# -------------------------
drop_predictors <- c(
  "eid",
  "cvd_event", "cvd_prevalent", "cvd_incident", "cvd_first_date",
  "dod", "age_of_death"
)
drop_predictors <- intersect(drop_predictors, names(ukb))
if (length(drop_predictors) > 0) {
  cat("Dropping predictors:", paste(drop_predictors, collapse = ", "), "\n")
  flush.console()
}

# build imputation frame: only numeric/factor/logical + targets
keep_cols <- names(ukb)[sapply(ukb, function(z) is.numeric(z) || is.factor(z) || is.logical(z))]
keep_cols <- unique(c(keep_cols, imp_vars))

# remove leakage predictors from keep_cols (key change)
keep_cols <- setdiff(keep_cols, drop_predictors)

ukb_imp <- ukb[, keep_cols, drop = FALSE]

# only targets that actually have missing
imp_vars <- imp_vars[sapply(imp_vars, function(v) any(is.na(ukb_imp[[v]])))]
if (length(imp_vars) == 0) {
  cat("No missing in target vars. Saving...\n")
  saveRDS(ukb, out_file)
  quit(save = "no")
}

miss_before <- sapply(imp_vars, function(v) sum(is.na(ukb_imp[[v]])))
cat("Missing before (top):\n")
print(head(sort(miss_before, decreasing = TRUE), 10))
flush.console()

MAXITER <- 5
NTREES  <- 75

# -------------------------
# batch 1
# -------------------------
batch1_vars <- intersect(
  c("smoking_3cat","crp","diabetes_bin","total_met_min_wk","tv_hours_day",
    "sleep_quality_score","mh_satis_mean_score","coffee_intake","tea_intake",
    "fruit_intake_fresh"),
  imp_vars
)
cat("Batch1:", paste(batch1_vars, collapse = ", "), "\n")
flush.console()

if (length(batch1_vars) > 0) {
  missing1 <- setdiff(batch1_vars, names(ukb_imp))
  if (length(missing1)) stop("Batch1 vars missing from ukb_imp: ", paste(missing1, collapse = ", "))
  
  imp1 <- miceRanger(
    data = ukb_imp,
    m = 1,
    maxiter = MAXITER,
    num.trees = NTREES,
    verbose = TRUE,
    vars = batch1_vars,
    valueSelector = "value",
    num.threads = n_cores
  )
  
  cd1 <- completeData(imp1)
  if (!is.list(cd1) || length(cd1) < 1) stop("completeData(imp1) returned unexpected structure.")
  comp1 <- as.data.frame(cd1[[1]])
  
  missing1b <- setdiff(batch1_vars, names(comp1))
  if (length(missing1b)) stop("Batch1 vars missing from completed data: ", paste(missing1b, collapse = ", "))
  
  ukb[, batch1_vars]     <- comp1[, batch1_vars, drop = FALSE]
  ukb_imp[, batch1_vars] <- comp1[, batch1_vars, drop = FALSE]
}

# rebuild ukb_imp safely from current ukb (keep same predictor set)
keep_cols2 <- intersect(keep_cols, names(ukb))
ukb_imp <- ukb[, keep_cols2, drop = FALSE]

# -------------------------
# batch 2
# (run everything EXCEPT oily_fish_3cat first; then run oily_fish_3cat alone)
# -------------------------
batch2_vars_all  <- setdiff(imp_vars, batch1_vars)
batch2_vars_main <- setdiff(batch2_vars_all, "oily_fish_3cat")
batch2_vars_oily <- intersect(batch2_vars_all, "oily_fish_3cat")

cat("Batch2(main):", paste(batch2_vars_main, collapse = ", "), "\n")
flush.console()

if (length(batch2_vars_main) > 0) {
  missing2 <- setdiff(batch2_vars_main, names(ukb_imp))
  if (length(missing2)) {
    cat("DEBUG: missing batch2(main) vars:\n")
    print(missing2)
    stop("Batch2(main) vars missing from ukb_imp.")
  }
  
  imp2 <- miceRanger(
    data = ukb_imp,
    m = 1,
    maxiter = MAXITER,
    num.trees = NTREES,
    verbose = TRUE,
    vars = batch2_vars_main,
    valueSelector = "value",
    num.threads = n_cores
  )
  
  cd2 <- completeData(imp2)
  if (!is.list(cd2) || length(cd2) < 1) stop("completeData(imp2) returned unexpected structure.")
  comp2 <- as.data.frame(cd2[[1]])
  
  missing2b <- setdiff(batch2_vars_main, names(comp2))
  if (length(missing2b)) stop("Batch2(main) vars missing from completed data: ", paste(missing2b, collapse = ", "))
  
  ukb[, batch2_vars_main] <- comp2[, batch2_vars_main, drop = FALSE]
  ukb_imp[, batch2_vars_main] <- comp2[, batch2_vars_main, drop = FALSE]
}

cat("Batch2(oily):", paste(batch2_vars_oily, collapse = ", "), "\n")
flush.console()

if (length(batch2_vars_oily) > 0) {
  if ("oily_fish_3cat" %in% names(ukb_imp)) {
    ukb_imp$oily_fish_3cat <- factor(ukb_imp$oily_fish_3cat, ordered = FALSE)
  }
  
  cat("PRE-OILY: has oily_fish_3cat =", "oily_fish_3cat" %in% names(ukb_imp),
      "| class =", paste(class(ukb_imp$oily_fish_3cat), collapse=","), "\n")
  flush.console()
  
  imp_oily <- miceRanger(
    data = ukb_imp,
    m = 1,
    maxiter = MAXITER,
    num.trees = NTREES,
    verbose = TRUE,
    vars = batch2_vars_oily,
    valueSelector = "value",
    num.threads = n_cores
  )
  
  cd_oily <- completeData(imp_oily)
  if (!is.list(cd_oily) || length(cd_oily) < 1) stop("completeData(imp_oily) returned unexpected structure.")
  comp_oily <- as.data.frame(cd_oily[[1]])
  
  if (!("oily_fish_3cat" %in% names(comp_oily))) stop("oily_fish_3cat missing from completed data (oily run).")
  ukb[, "oily_fish_3cat"] <- comp_oily[, "oily_fish_3cat", drop = FALSE]
}

# -------------------------
# regenerate derived group vars
# -------------------------
if ("tv_hours_day" %in% names(ukb)) {
  tv_num <- ukb$tv_hours_day
  tv_num[tv_num < 0 | tv_num > 24] <- NA_real_
  ukb$tv_group <- factor(
    ifelse(is.na(tv_num), NA,
           ifelse(tv_num < 1, "<1h",
                  ifelse(tv_num <= 3, "1–3h",
                         ifelse(tv_num <= 5, "3–5h", ">5h")))),
    levels = c("<1h", "1–3h", "3–5h", ">5h")
  )
}

if ("total_met_min_wk" %in% names(ukb)) {
  met_raw <- ukb$total_met_min_wk
  met_raw[met_raw < 0] <- NA_real_
  ukb$total_met_group <- factor(
    ifelse(is.na(met_raw), NA,
           ifelse(met_raw < 600, "Low",
                  ifelse(met_raw < 3000, "Moderate", "High"))),
    levels = c("Low", "Moderate", "High")
  )
}

miss_after <- sapply(imp_vars, function(v) sum(is.na(ukb[[v]])))
print(data.frame(variable = imp_vars, missing_before = miss_before, missing_after = miss_after))

saveRDS(ukb, out_file)
cat("Saved:", out_file, "\n")
flush.console()