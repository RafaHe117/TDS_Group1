suppressPackageStartupMessages(library(miceRanger))

setwd("/rds/general/project/hda_25-26/live/TDS/TDS_Group1")

in_file  <- "ukb_G1_cleaned.rds"
out_dir  <- "imputation"
out_file <- file.path(out_dir, "ukb_G1_imputed.rds")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(20260225)

ukb <- readRDS(in_file)
ukb <- as.data.frame(ukb)

# ordered -> factor
ord <- names(ukb)[sapply(ukb, is.ordered)]
if (length(ord) > 0) {
  for (v in ord) {
    ukb[[v]] <- factor(ukb[[v]], ordered = FALSE)
  }
}

# cores
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

# protect ID / outcome / date-like
id_like      <- grep("^eid$|id$", names(ukb), ignore.case = TRUE, value = TRUE)
date_like    <- grep("date|death|time", names(ukb), ignore.case = TRUE, value = TRUE)
outcome_like <- grep("^cvd", names(ukb), ignore.case = TRUE, value = TRUE)

non_impute <- unique(c(id_like, date_like, outcome_like))

# convert character -> factor (excluding protected vars)
char_vars <- names(ukb)[sapply(ukb, is.character)]
char_vars <- setdiff(char_vars, non_impute)
for (v in char_vars) {
  ukb[[v]] <- factor(ukb[[v]])
}

# eligible vars
eligible <- names(ukb)[
  sapply(ukb, function(x) is.numeric(x) || is.factor(x) || is.logical(x))
]

# targets = any missing (1 missing counts)
imp_vars <- eligible[
  sapply(eligible, function(v) any(is.na(ukb[[v]])))
]

# remove protected vars
imp_vars <- setdiff(imp_vars, non_impute)

if (length(imp_vars) == 0) {
  saveRDS(ukb, out_file)
  quit(save = "no")
}

cat("Total imputation targets:", length(imp_vars), "\n")

# predictor frame (exclude protected)
predictor_cols <- setdiff(eligible, non_impute)
ukb_imp <- ukb[, predictor_cols, drop = FALSE]

imp_vars <- intersect(imp_vars, names(ukb_imp))

# hyperparameters
MAXITER <- 10
NTREES  <- 100

# detect sparse categorical vars
sparse_cat <- imp_vars[
  sapply(imp_vars, function(v)
    is.factor(ukb_imp[[v]]) &&
      length(unique(na.omit(ukb_imp[[v]]))) < 3
  )
]

main_vars <- setdiff(imp_vars, sparse_cat)

# main batch
if (length(main_vars) > 0) {
  imp_main <- miceRanger(
    data = ukb_imp,
    m = 1,
    maxiter = MAXITER,
    num.trees = NTREES,
    vars = main_vars,
    valueSelector = "meanMatch",
    num.threads = n_cores,
    verbose = TRUE
  )
  
  comp_main <- as.data.frame(completeData(imp_main)[[1]])
  ukb[, main_vars]     <- comp_main[, main_vars, drop = FALSE]
  ukb_imp[, main_vars] <- comp_main[, main_vars, drop = FALSE]
}

# sparse factors (single runs)
if (length(sparse_cat) > 0) {
  for (v in sparse_cat) {
    imp_single <- miceRanger(
      data = ukb_imp,
      m = 1,
      maxiter = MAXITER,
      num.trees = NTREES,
      vars = v,
      valueSelector = "meanMatch",
      num.threads = n_cores,
      verbose = TRUE
    )
    
    comp_single <- as.data.frame(completeData(imp_single)[[1]])
    ukb[[v]]     <- comp_single[[v]]
    ukb_imp[[v]] <- comp_single[[v]]
  }
}

# QC
miss_after <- sapply(imp_vars, function(v) sum(is.na(ukb[[v]])))
cat("Remaining missing (should be 0):\n")
print(miss_after[miss_after > 0])

stopifnot(all(sapply(imp_vars, function(v) sum(is.na(ukb[[v]]))) == 0))

saveRDS(ukb, out_file)
cat("Saved:", out_file, "\n")

gc()