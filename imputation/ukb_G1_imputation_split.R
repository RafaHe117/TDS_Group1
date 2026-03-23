suppressPackageStartupMessages({
  library(miceRanger)
})

setwd("/rds/general/project/hda_25-26/live/TDS/TDS_Group1")

in_file <- "ukb_G1_cleaned.rds"
out_dir <- "imputation"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

train_file <- file.path(out_dir, "ukb_G1_train_imputed.rds")
test_file  <- file.path(out_dir, "ukb_G1_test_imputed.rds")

set.seed(20260225)

outcome <- "cvd_incident"
train_p <- 0.7
MAXITER <- 10
NTREES  <- 100

ukb <- readRDS(in_file)
ukb <- as.data.frame(ukb)

if (!outcome %in% names(ukb)) {
  stop("Outcome variable not found: ", outcome)
}

ord <- names(ukb)[sapply(ukb, is.ordered)]
if (length(ord) > 0) {
  for (v in ord) ukb[[v]] <- factor(ukb[[v]], ordered = FALSE)
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

id_like      <- grep("^eid$|id$", names(ukb), ignore.case = TRUE, value = TRUE)
date_like    <- grep("date|death|time", names(ukb), ignore.case = TRUE, value = TRUE)
outcome_like <- outcome
non_impute   <- unique(c(id_like, date_like, outcome_like))

char_vars <- names(ukb)[sapply(ukb, is.character)]
char_vars <- setdiff(char_vars, non_impute)
for (v in char_vars) ukb[[v]] <- factor(ukb[[v]])

y <- ukb[[outcome]]
if (any(is.na(y))) {
  stop("Outcome contains NA. Please handle missing outcome before split.")
}
y_fac <- if (is.factor(y)) y else factor(y)

train_idx <- unlist(
  lapply(split(seq_len(nrow(ukb)), y_fac), function(idx) {
    n_cls <- length(idx)
    n_tr  <- round(train_p * n_cls)
    sample(idx, size = n_tr, replace = FALSE)
  }),
  use.names = FALSE
)

train_idx <- sort(train_idx)
test_idx  <- setdiff(seq_len(nrow(ukb)), train_idx)

train_data <- ukb[train_idx, , drop = FALSE]
test_data  <- ukb[test_idx,  , drop = FALSE]

cat("Train rows:", nrow(train_data), "\n")
cat("Test rows :", nrow(test_data), "\n")
cat("Outcome proportion in full data:\n")
print(prop.table(table(ukb[[outcome]])))
cat("Outcome proportion in train:\n")
print(prop.table(table(train_data[[outcome]])))
cat("Outcome proportion in test:\n")
print(prop.table(table(test_data[[outcome]])))

get_imp_setup <- function(dat, non_impute) {
  eligible <- names(dat)[
    sapply(dat, function(x) is.numeric(x) || is.factor(x) || is.logical(x))
  ]
  
  imp_vars <- eligible[
    sapply(eligible, function(v) any(is.na(dat[[v]])))
  ]
  
  imp_vars <- setdiff(imp_vars, non_impute)
  predictor_cols <- setdiff(eligible, non_impute)
  
  list(
    eligible = eligible,
    predictor_cols = predictor_cols,
    imp_vars = imp_vars
  )
}

setup_train    <- get_imp_setup(train_data, non_impute)
predictor_cols <- setup_train$predictor_cols
imp_vars_train <- setup_train$imp_vars

cat("Number of predictor columns:", length(predictor_cols), "\n")
cat("Outcome in predictors:", outcome %in% predictor_cols, "\n")
cat("EID in predictors:", any(grepl("^eid$", predictor_cols, ignore.case = TRUE)), "\n")
cat("Any ID-like var in predictors:", any(grepl("^eid$|id$", predictor_cols, ignore.case = TRUE)), "\n")

if (length(imp_vars_train) == 0) {
  saveRDS(train_data, train_file)
  saveRDS(test_data, test_file)
  
  cat("No train variables need imputation.\n")
  cat("Saved:\n")
  cat(" -", train_file, "\n")
  cat(" -", test_file, "\n")
  
  gc()
  quit(save = "no")
}

cat("Total train imputation targets:", length(imp_vars_train), "\n")

train_imp <- train_data[, predictor_cols, drop = FALSE]
test_imp  <- test_data[, predictor_cols, drop = FALSE]

fac_cols <- predictor_cols[sapply(train_imp, is.factor)]
for (v in fac_cols) {
  test_imp[[v]] <- factor(test_imp[[v]], levels = levels(train_imp[[v]]))
}

imp_vars <- intersect(imp_vars_train, names(train_imp))

sparse_cat <- imp_vars[
  sapply(imp_vars, function(v) {
    is.factor(train_imp[[v]]) &&
      length(unique(na.omit(train_imp[[v]]))) < 3
  })
]

main_vars <- setdiff(imp_vars, sparse_cat)

if (length(main_vars) > 0) {
  cat("Main vars:", length(main_vars), "\n")
  
  imp_main_train <- miceRanger(
    data = train_imp,
    m = 1,
    maxiter = MAXITER,
    num.trees = NTREES,
    vars = main_vars,
    valueSelector = "meanMatch",
    num.threads = n_cores,
    verbose = TRUE,
    saveModels = TRUE,
    returnModels = TRUE,
    returnRF = TRUE
  )
  
  comp_train_main <- as.data.frame(completeData(imp_main_train)[[1]])
  train_data[, main_vars] <- comp_train_main[, main_vars, drop = FALSE]
  train_imp[, main_vars]  <- comp_train_main[, main_vars, drop = FALSE]
  
  imp_main_test <- miceRanger(
    data = test_imp,
    m = 1,
    maxiter = MAXITER,
    num.trees = NTREES,
    vars = main_vars,
    valueSelector = "meanMatch",
    num.threads = n_cores,
    verbose = TRUE,
    saveModels = FALSE,
    returnModels = FALSE,
    returnRF = FALSE,
    forest = imp_main_train$rfImpute
  )
  
  comp_test_main <- as.data.frame(completeData(imp_main_test)[[1]])
  test_data[, main_vars] <- comp_test_main[, main_vars, drop = FALSE]
  test_imp[, main_vars]  <- comp_test_main[, main_vars, drop = FALSE]
}

if (length(sparse_cat) > 0) {
  cat("Sparse categorical vars:", length(sparse_cat), "\n")
  
  for (v in sparse_cat) {
    cat("Imputing sparse variable:", v, "\n")
    
    imp_single_train <- miceRanger(
      data = train_imp,
      m = 1,
      maxiter = MAXITER,
      num.trees = NTREES,
      vars = v,
      valueSelector = "meanMatch",
      num.threads = n_cores,
      verbose = TRUE,
      saveModels = TRUE,
      returnModels = TRUE,
      returnRF = TRUE
    )
    
    comp_single_train <- as.data.frame(completeData(imp_single_train)[[1]])
    train_data[[v]] <- comp_single_train[[v]]
    train_imp[[v]]  <- comp_single_train[[v]]
    
    test_imp[[v]] <- factor(test_imp[[v]], levels = levels(train_imp[[v]]))
    
    imp_single_test <- miceRanger(
      data = test_imp,
      m = 1,
      maxiter = MAXITER,
      num.trees = NTREES,
      vars = v,
      valueSelector = "meanMatch",
      num.threads = n_cores,
      verbose = TRUE,
      saveModels = FALSE,
      returnModels = FALSE,
      returnRF = FALSE,
      forest = imp_single_train$rfImpute
    )
    
    comp_single_test <- as.data.frame(completeData(imp_single_test)[[1]])
    test_data[[v]] <- comp_single_test[[v]]
    test_imp[[v]]  <- comp_single_test[[v]]
  }
}

miss_after_train <- sapply(imp_vars, function(v) sum(is.na(train_data[[v]])))
miss_after_test  <- sapply(imp_vars, function(v) sum(is.na(test_data[[v]])))

cat("Remaining missing in train:\n")
print(miss_after_train[miss_after_train > 0])

cat("Remaining missing in test:\n")
print(miss_after_test[miss_after_test > 0])

stopifnot(all(sapply(imp_vars, function(v) sum(is.na(train_data[[v]]))) == 0))

saveRDS(train_data, train_file)
saveRDS(test_data, test_file)

cat("Saved:\n")
cat(" -", train_file, "\n")
cat(" -", test_file, "\n")

gc()