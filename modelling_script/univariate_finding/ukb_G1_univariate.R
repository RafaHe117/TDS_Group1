suppressPackageStartupMessages({
  library(dplyr)
  library(broom)
  library(parallel)
})

in_dir  <- "/rds/general/project/hda_25-26/live/TDS/TDS_Group1"
out_dir <- file.path(in_dir, "modelling_script", "univariate_finding")
uni_dir <- file.path(out_dir, "univariate_output")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(uni_dir, showWarnings = FALSE, recursive = TRUE)

setwd(in_dir)

x <- readRDS("ukb_G1_imputed_final.rds")

# convert character variables to factor (safe guard for real data)
char_vars <- names(x)[sapply(x, is.character)]
x[char_vars] <- lapply(x[char_vars], factor)

outcome <- "cvd_incident"
confounders <- c("age", "sex", "ethnicity_5cat")

exclude_vars <- c(
  "eid", "cvd_incident", "cvd_event", "cvd_prevalent",
  "cvd_first_date", "dod", "age_of_death", "date_recr", "yob"
)

predictors <- setdiff(names(x), c(exclude_vars, confounders))

factor_vars <- c(
  "sex", "ethnicity_5cat", "urban_rural_2cat", "famhx_cvd",
  "mood_disorder", "diabetes_bin", "smoking_3cat",
  "education_2", "alcohol_combined", "redwine_group",
  "oily_fish_3cat", "salt_3cat", "living_with_partner",
  "employment_2cat", "imd_quintile", "medication_category",
  "ICD_Diabetes", "ICD_Lipidemia", "ICD_Mental_Health",
  "ICD_Migraine", "ICD_Hypertension", "ICD_AF",
  "ICD_Atopy", "ICD_Autoimmune", "ICD_CKD",
  "tv_group", "total_met_group", "stress_group_2yr"
)

# factor convert
factor_vars <- intersect(factor_vars, names(x))
x[factor_vars] <- lapply(x[factor_vars], factor)

max_levels <- 50
alpha_level <- 0.05

# -------- helper functions --------
is_binary_numeric <- function(v) {
  if (!is.numeric(v)) return(FALSE)
  uv <- unique(v[!is.na(v)])
  length(uv) <= 2
}

is_continuous_numeric <- function(v) {
  if (!is.numeric(v)) return(FALSE)
  uv <- unique(v[!is.na(v)])
  if (length(uv) <= 2) return(FALSE)
  s <- sd(v, na.rm = TRUE)
  if (is.na(s) || s == 0) return(FALSE)
  TRUE
}

has_variation <- function(v) {
  uv <- unique(v[!is.na(v)])
  length(uv) > 1
}

make_bh_detail_table <- function(df, p_col, id_cols, alpha = 0.05) {
  out <- df %>%
    select(all_of(id_cols), raw_p_value = all_of(p_col)) %>%
    mutate(raw_p_value = pmax(raw_p_value, .Machine$double.xmin)) %>%
    arrange(raw_p_value) %>%
    mutate(
      rank = row_number(),
      n_tests = n(),
      BH_critical_value = (rank / n_tests) * alpha,
      is_BH_significant = raw_p_value <= BH_critical_value,
      BH_adjusted_p_value = p.adjust(raw_p_value, method = "BH"),
      Bonferroni_adjusted_p_value = p.adjust(raw_p_value, method = "bonferroni")
    )
  out
}

# -------- auto-detect variable types --------
continuous_numeric_predictors <- predictors[
  sapply(x[predictors], is_continuous_numeric)
]

binary_numeric_predictors <- predictors[
  sapply(x[predictors], is_binary_numeric)
]

factor_predictors <- predictors[
  sapply(x[predictors], is.factor)
]

other_predictors <- setdiff(
  predictors,
  c(continuous_numeric_predictors, binary_numeric_predictors, factor_predictors)
)

message("========================================")
message("AUTO DETECTION SUMMARY")
message("========================================")

message("Continuous numeric predictors (will be standardized):")
print(continuous_numeric_predictors)

message(" ")

message("Binary numeric predictors (NOT standardized):")
print(binary_numeric_predictors)

message(" ")

message("Factor predictors:")
print(factor_predictors)

message(" ")

message("Other predictors:")
print(other_predictors)

message(" ")

message("Counts:")
message("Continuous numeric: ", length(continuous_numeric_predictors))
message("Binary numeric: ", length(binary_numeric_predictors))
message("Factor: ", length(factor_predictors))
message("Other: ", length(other_predictors))
message("========================================")

if (length(continuous_numeric_predictors) > 0) {
  message("Preview of continuous numeric predictors before standardization:")
  print(head(continuous_numeric_predictors, 20))
  
  precheck_df <- data.frame(
    predictor = continuous_numeric_predictors,
    mean_before = sapply(x[continuous_numeric_predictors], function(v) mean(v, na.rm = TRUE)),
    sd_before   = sapply(x[continuous_numeric_predictors], function(v) sd(v, na.rm = TRUE))
  )
  print(head(precheck_df, 20))
}

# -------- standardize continuous numeric predictors only --------
if (length(continuous_numeric_predictors) > 0) {
  x[continuous_numeric_predictors] <- lapply(
    x[continuous_numeric_predictors],
    function(v) as.numeric(scale(v))
  )
}

# -------- post-standardization quick check --------
if (length(continuous_numeric_predictors) > 0) {
  message("========================================")
  message("POST-STANDARDIZATION CHECK")
  message("========================================")
  
  postcheck_df <- data.frame(
    predictor = continuous_numeric_predictors,
    mean_after = sapply(x[continuous_numeric_predictors], function(v) mean(v, na.rm = TRUE)),
    sd_after   = sapply(x[continuous_numeric_predictors], function(v) sd(v, na.rm = TRUE))
  )
  print(head(postcheck_df, 20))
  
  message("Note: mean_after should be close to 0, sd_after close to 1.")
  message("========================================")
}

# optional: drop predictors with no variation
predictors <- predictors[sapply(x[predictors], has_variation)]

# cores
n_cores <- as.integer(Sys.getenv("MY_NCORES", unset = "1"))
if (is.na(n_cores) || n_cores < 1) n_cores <- 1

message("========================================")
message("Start: ", Sys.time())
message("Cores: ", n_cores)
message("Predictors after filtering: ", length(predictors))
message("========================================")

run_one <- function(i) {
  
  var <- predictors[i]
  total <- length(predictors)
  pid <- Sys.getpid()
  
  message(sprintf("[%s] START %d/%d : %s (pid=%s)",
                  Sys.time(), i, total, var, pid))
  flush.console()
  
  if (!var %in% names(x)) {
    message(sprintf("[%s] SKIP  %d/%d : %s (missing)",
                    Sys.time(), i, total, var))
    flush.console()
    return(NULL)
  }
  
  if (is.factor(x[[var]]) && nlevels(x[[var]]) > max_levels) {
    message(sprintf("[%s] SKIP  %d/%d : %s (levels>%d)",
                    Sys.time(), i, total, var, max_levels))
    flush.console()
    return(NULL)
  }
  
  fml <- as.formula(
    paste0(outcome, " ~ ", var, " + ", paste(confounders, collapse = " + "))
  )
  
  fit <- try(glm(fml, data = x, family = binomial()), silent = TRUE)
  if (inherits(fit, "try-error")) {
    message(sprintf("[%s] FAIL  %d/%d : %s (glm)",
                    Sys.time(), i, total, var))
    flush.console()
    return(NULL)
  }
  
  td <- try(broom::tidy(fit, conf.int = TRUE), silent = TRUE)
  if (inherits(td, "try-error")) {
    message(sprintf("[%s] FAIL  %d/%d : %s (tidy)",
                    Sys.time(), i, total, var))
    flush.console()
    return(NULL)
  }
  
  td <- td[td$term == var | grepl(paste0("^", var), td$term), , drop = FALSE]
  if (nrow(td) == 0) {
    message(sprintf("[%s] SKIP  %d/%d : %s (no term)",
                    Sys.time(), i, total, var))
    flush.console()
    return(NULL)
  }
  
  out <- td %>%
    mutate(
      predictor = var,
      predictor_type = case_when(
        is.factor(x[[var]]) ~ "factor",
        var %in% continuous_numeric_predictors ~ "continuous_numeric_standardized",
        var %in% binary_numeric_predictors ~ "binary_numeric_not_standardized",
        is.numeric(x[[var]]) ~ "numeric_not_standardized",
        TRUE ~ "other"
      ),
      beta = estimate,
      beta_CI_low = conf.low,
      beta_CI_high = conf.high,
      OR = exp(estimate),
      OR_CI_low = exp(conf.low),
      OR_CI_high = exp(conf.high)
    ) %>%
    select(
      predictor, predictor_type, term,
      beta, beta_CI_low, beta_CI_high,
      OR, OR_CI_low, OR_CI_high,
      p.value
    )
  
  message(sprintf("[%s] DONE  %d/%d : %s",
                  Sys.time(), i, total, var))
  flush.console()
  
  out
}

message("[", Sys.time(), "] parallel start")

res <- mclapply(
  seq_along(predictors),
  run_one,
  mc.cores = n_cores,
  mc.preschedule = FALSE
)

message("[", Sys.time(), "] parallel end")

res <- Filter(Negate(is.null), res)

if (length(res) == 0) {
  stop("No valid results returned.")
}

final_results_term <- bind_rows(res) %>%
  mutate(
    raw_p_value = pmax(p.value, .Machine$double.xmin),
    FDR_p_value = p.adjust(raw_p_value, method = "BH"),
    Bonferroni_p_value = p.adjust(raw_p_value, method = "bonferroni")
  ) %>%
  select(
    predictor, predictor_type, term,
    beta, beta_CI_low, beta_CI_high,
    OR, OR_CI_low, OR_CI_high,
    raw_p_value, FDR_p_value, Bonferroni_p_value
  ) %>%
  arrange(FDR_p_value, raw_p_value)

write.csv(
  final_results_term,
  file.path(uni_dir, "univariate_cvd_results_term_level.csv"),
  row.names = FALSE
)

final_results_predictor <- final_results_term %>%
  group_by(predictor, predictor_type) %>%
  summarise(
    n_terms = n(),
    top_term = term[which.min(raw_p_value)],
    top_beta = beta[which.min(raw_p_value)],
    top_beta_CI_low = beta_CI_low[which.min(raw_p_value)],
    top_beta_CI_high = beta_CI_high[which.min(raw_p_value)],
    top_OR = OR[which.min(raw_p_value)],
    top_OR_CI_low = OR_CI_low[which.min(raw_p_value)],
    top_OR_CI_high = OR_CI_high[which.min(raw_p_value)],
    predictor_raw_p_value = min(raw_p_value),
    .groups = "drop"
  ) %>%
  mutate(
    predictor_FDR_p_value = p.adjust(predictor_raw_p_value, method = "BH"),
    predictor_Bonferroni_p_value = p.adjust(predictor_raw_p_value, method = "bonferroni")
  ) %>%
  arrange(predictor_FDR_p_value, predictor_raw_p_value)

write.csv(
  final_results_predictor,
  file.path(uni_dir, "univariate_cvd_results_predictor_level.csv"),
  row.names = FALSE
)

# -------- additional outputs: multiple testing correction values --------

alpha_level <- 0.05

n_term_tests <- nrow(final_results_term)
n_predictor_tests <- nrow(final_results_predictor)

multiple_testing_summary <- data.frame(
  level = c("term_level", "predictor_level"),
  n_tests = c(n_term_tests, n_predictor_tests),
  alpha = c(alpha_level, alpha_level),
  bonferroni_threshold = c(
    alpha_level / n_term_tests,
    alpha_level / n_predictor_tests
  ),
  n_BH_significant = c(
    sum(final_results_term$FDR_p_value < alpha_level, na.rm = TRUE),
    sum(final_results_predictor$predictor_FDR_p_value < alpha_level, na.rm = TRUE)
  ),
  n_Bonferroni_significant = c(
    sum(final_results_term$Bonferroni_p_value < alpha_level, na.rm = TRUE),
    sum(final_results_predictor$predictor_Bonferroni_p_value < alpha_level, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  multiple_testing_summary,
  file.path(uni_dir, "multiple_testing_summary.csv"),
  row.names = FALSE
)

message("Saved:")
message(" - ", file.path(uni_dir, "univariate_cvd_results_term_level.csv"))
message(" - ", file.path(uni_dir, "univariate_cvd_results_predictor_level.csv"))
message(" - ", file.path(uni_dir, "multiple_testing_summary.csv"))

message("Finished: ", Sys.time())