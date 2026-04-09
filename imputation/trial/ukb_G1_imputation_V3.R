suppressPackageStartupMessages(library(miceRanger))

# paths
setwd("/rds/general/project/hda_25-26/live/TDS/TDS_Group1")
in_file  <- "ukb_G1_cleaned.rds"
out_dir  <- "imputation"
out_file <- file.path(out_dir, "ukb_G1_imputed_final.rds")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(20260225)

# load
ukb <- readRDS(in_file)
ukb <- as.data.frame(ukb)

# drop derived vars if present (will regenerate later)
if ("tv_group" %in% names(ukb)) ukb$tv_group <- NULL
if ("total_met_group" %in% names(ukb)) ukb$total_met_group <- NULL

# convert ALL ordered factors -> plain factors
is_ord <- sapply(ukb, is.ordered)
if (any(is_ord)) {
  ord_names <- names(ukb)[is_ord]
  for (v in ord_names) ukb[[v]] <- factor(ukb[[v]], ordered = FALSE)
  cat("Converted ordered factors:", paste(ord_names, collapse = ", "), "\n")
  flush.console()
}

# ensure unique column names
stopifnot(!anyDuplicated(names(ukb)))

# threads
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

# vars to impute
imp_vars <- c(
  "smoking_3cat","crp","diabetes_bin","basophil_count",
  "oily_fish_3cat","lipoprotein_a","total_met_min_wk",
  "sf_score","tv_hours_day","salt_3cat","eosinophil_count",
  "sleep_quality_score","mh_satis_mean_score","creatinine","coffee_intake",
  "tea_intake","fruit_intake_fresh","alanine_aminotransferase","lymphocyte_count",
  "total_triglyceride","monocyte_count","redwine_group","greenspace_pct_1000m",
  "total_bilirubin","direct_bilirubin","blood_vitamin_d","aspartate_aminotransferase",
  "imd_final","imd_source","imd_quintile","neutrophil_count","alkaline_phosphatase",
  "wbc_count","no2_2010","urea","hdl_cholesterol","igf_1","urate",
  "employment_2cat","work_hours_week_clean","platelet_count","glucose",
  "ldl_cholesterol","birth_weight_clean","living_with_partner"
)

# keep only variables that exist
imp_vars <- intersect(imp_vars, names(ukb))
if (length(imp_vars) == 0) stop("None of imp_vars exist in dataset.")

# factors (force to plain factor, not ordered)
cat_vars <- intersect(
  c("smoking_3cat","diabetes_bin","oily_fish_3cat","salt_3cat","redwine_group",
    "imd_quintile","employment_2cat","imd_source","living_with_partner"),
  names(ukb)
)
for (v in cat_vars) ukb[[v]] <- factor(ukb[[v]], ordered = FALSE)

# keep usable columns only (numeric/factor/logical)
keep_cols <- names(ukb)[sapply(ukb, function(z) is.numeric(z) || is.factor(z) || is.logical(z))]
keep_cols <- union(keep_cols, imp_vars)
ukb_imp <- ukb[, keep_cols, drop = FALSE]

# ensure target vars are in ukb_imp (defensive)
imp_vars <- intersect(imp_vars, names(ukb_imp))
if (length(imp_vars) == 0) stop("Target vars missing from ukb_imp after filtering.")

# only vars with missing
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

# settings
MAXITER <- 10
NTREES  <- 100

# run imputation
cat("Imputing vars:", paste(imp_vars, collapse = ", "), "\n")
flush.console()

impute <- miceRanger(
  data = ukb_imp,
  m = 1,
  maxiter = MAXITER,
  num.trees = NTREES,
  verbose = TRUE,
  vars = imp_vars,
  valueSelector = "value",
  num.threads = n_cores
)

# complete + coerce to data.frame to avoid data.table column-selection issues
cd <- completeData(impute)
stopifnot(is.list(cd), length(cd) >= 1)
comp <- as.data.frame(cd[[1]])

# defensive: only update columns that exist in comp
imp_vars2 <- intersect(imp_vars, names(comp))
if (length(imp_vars2) == 0) stop("No imputed columns found in completed data.")

ukb[, imp_vars2] <- comp[, imp_vars2, drop = FALSE]

# derived groups
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

# save
miss_after <- sapply(imp_vars2, function(v) sum(is.na(ukb[[v]])))
print(data.frame(variable = imp_vars2,
                 missing_before = miss_before[imp_vars2],
                 missing_after = miss_after))

saveRDS(ukb, out_file)
cat("Saved:", out_file, "\n")
flush.console()