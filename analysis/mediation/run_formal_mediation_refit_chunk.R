suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
})

BASE_DIR <- "/rds/general/project/hda_25-26/live/TDS/TDS_Group1"
MED_DIR  <- file.path(BASE_DIR, "analysis", "mediation")
DATA_PATH <- file.path(BASE_DIR, "split_imputed_data", "ukb_G1_train_imputed.rds")

INPUT_DIR <- file.path(MED_DIR, "inputs")
OUT_DIR   <- file.path(MED_DIR, "outputs", "formal_mediation_refit_split")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

CONFIG_PATH <- file.path(INPUT_DIR, "formal_mediation_config.csv")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript run_formal_mediation_refit_chunk.R <pair_csv> <job_label> <subgroup_label>")
}

PAIR_CSV <- args[1]
JOB_LABEL <- args[2]
SUBGROUP_LABEL <- args[3]   # use: all / female / male

# --------------------------------------------------
# config
# --------------------------------------------------
cfg <- read_csv(CONFIG_PATH, show_col_types = FALSE)

get_cfg <- function(k) {
  val <- cfg %>% filter(key == k) %>% pull(value)
  if (length(val) != 1 || is.na(val)) stop("Missing config key: ", k)
  val
}

OUTCOME <- get_cfg("outcome")
SEX_VAR <- get_cfg("sex_var")
CONFOUNDERS <- strsplit(get_cfg("confounders"), "\\|")[[1]]
CONFOUNDERS <- CONFOUNDERS[CONFOUNDERS != ""]
N_BOOT <- as.integer(get_cfg("n_boot"))
MIN_COMPLETE_N <- as.integer(get_cfg("min_complete_n"))

# --------------------------------------------------
# helpers
# --------------------------------------------------
norm_txt <- function(x) {
  x <- tolower(as.character(x))
  gsub("[^a-z0-9]+", "", x)
}

infer_subset_value <- function(label, data_vec) {
  if (norm_txt(label) %in% c("all", "main", "overall")) return(NULL)
  
  vals <- unique(as.character(data_vec))
  vals <- vals[!is.na(vals)]
  vals_norm <- norm_txt(vals)
  label_norm <- norm_txt(label)
  
  alias_map <- list(
    male   = c("male", "m", "men", "man"),
    female = c("female", "f", "women", "woman")
  )
  
  target_alias <- alias_map[[label_norm]]
  if (is.null(target_alias)) {
    hit <- vals[vals_norm == label_norm]
    if (length(hit) == 1) return(hit)
    stop("Could not infer subgroup value for label=", label,
         " from values: ", paste(vals, collapse = ", "))
  }
  
  hit <- vals[vals_norm %in% target_alias]
  if (length(hit) == 1) return(hit)
  
  hit <- vals[sapply(vals_norm, function(z) any(stringr::str_detect(z, target_alias)))]
  if (length(hit) == 1) return(hit)
  
  stop("Could not infer subgroup value for label=", label,
       " from values: ", paste(vals, collapse = ", "))
}

safe_coef <- function(fit, coef_name) {
  cf <- summary(fit)$coefficients
  if (!(coef_name %in% rownames(cf))) {
    return(c(estimate = NA_real_, se = NA_real_, p = NA_real_))
  }
  c(
    estimate = unname(cf[coef_name, 1]),
    se       = unname(cf[coef_name, 2]),
    p        = unname(cf[coef_name, 4])
  )
}

get_complete_data <- function(df, vars) {
  vars <- unique(vars)
  vars <- vars[vars %in% names(df)]
  dat <- df[, vars, drop = FALSE]
  dat <- dat[complete.cases(dat), , drop = FALSE]
  dat
}

bootstrap_one_pair <- function(dat, x_var, x_term, m_var, y_var, confounders, n_boot = 1000) {
  n <- nrow(dat)
  indirect_vals <- rep(NA_real_, n_boot)
  direct_vals   <- rep(NA_real_, n_boot)
  total_vals    <- rep(NA_real_, n_boot)
  
  for (b in seq_len(n_boot)) {
    idx <- sample.int(n, size = n, replace = TRUE)
    d <- dat[idx, , drop = FALSE]
    
    fit_m <- try(
      glm(
        as.formula(paste(m_var, "~", paste(c(x_var, confounders), collapse = " + "))),
        data = d,
        family = gaussian()
      ),
      silent = TRUE
    )
    
    fit_y_total <- try(
      glm(
        as.formula(paste(y_var, "~", paste(c(x_var, confounders), collapse = " + "))),
        data = d,
        family = binomial()
      ),
      silent = TRUE
    )
    
    fit_y_direct <- try(
      glm(
        as.formula(paste(y_var, "~", paste(c(x_var, m_var, confounders), collapse = " + "))),
        data = d,
        family = binomial()
      ),
      silent = TRUE
    )
    
    if (inherits(fit_m, "try-error") ||
        inherits(fit_y_total, "try-error") ||
        inherits(fit_y_direct, "try-error")) {
      next
    }
    
    a_est  <- safe_coef(fit_m, x_term)["estimate"]
    b_est  <- safe_coef(fit_y_direct, m_var)["estimate"]
    c_est  <- safe_coef(fit_y_total, x_term)["estimate"]
    cp_est <- safe_coef(fit_y_direct, x_term)["estimate"]
    
    indirect_vals[b] <- a_est * b_est
    direct_vals[b]   <- cp_est
    total_vals[b]    <- c_est
  }
  
  tibble(
    indirect_boot_mean = mean(indirect_vals, na.rm = TRUE),
    indirect_ci_low    = quantile(indirect_vals, 0.025, na.rm = TRUE, names = FALSE),
    indirect_ci_high   = quantile(indirect_vals, 0.975, na.rm = TRUE, names = FALSE),
    direct_boot_mean   = mean(direct_vals, na.rm = TRUE),
    direct_ci_low      = quantile(direct_vals, 0.025, na.rm = TRUE, names = FALSE),
    direct_ci_high     = quantile(direct_vals, 0.975, na.rm = TRUE, names = FALSE),
    total_boot_mean    = mean(total_vals, na.rm = TRUE),
    total_ci_low       = quantile(total_vals, 0.025, na.rm = TRUE, names = FALSE),
    total_ci_high      = quantile(total_vals, 0.975, na.rm = TRUE, names = FALSE)
  )
}

run_one_pair <- function(df, analysis_label, subgroup_value, x_var, x_term, m_var) {
  dat <- df
  confounders <- CONFOUNDERS
  
  if (!is.null(subgroup_value)) {
    dat <- dat %>% filter(as.character(.data[[SEX_VAR]]) == subgroup_value)
    confounders <- setdiff(confounders, SEX_VAR)
  }
  
  required_vars <- unique(c(OUTCOME, x_var, m_var, confounders))
  missing_vars <- setdiff(required_vars, names(dat))
  if (length(missing_vars) > 0) {
    stop("Missing variables for pair ", x_term, " -> ", m_var, ": ",
         paste(missing_vars, collapse = ", "))
  }
  
  dat_cc <- get_complete_data(dat, c(OUTCOME, x_var, m_var, confounders))
  
  if (nrow(dat_cc) < MIN_COMPLETE_N) {
    return(tibble(
      analysis = analysis_label,
      subgroup_value = ifelse(is.null(subgroup_value), "all", subgroup_value),
      exposure_var = x_var,
      exposure_term = x_term,
      mediator = m_var,
      n_complete = nrow(dat_cc),
      a_estimate = NA_real_,
      a_se = NA_real_,
      a_p = NA_real_,
      b_estimate = NA_real_,
      b_se = NA_real_,
      b_p = NA_real_,
      total_effect_c = NA_real_,
      total_se = NA_real_,
      total_p = NA_real_,
      direct_effect_cprime = NA_real_,
      direct_se = NA_real_,
      direct_p = NA_real_,
      indirect_effect_ab = NA_real_,
      indirect_direction = NA_character_,
      note = paste0("Too few complete-case rows: ", nrow(dat_cc))
    ))
  }
  
  if (!is.numeric(dat_cc[[m_var]])) {
    stop("Mediator is not numeric: ", m_var)
  }
  
  fit_m <- glm(
    as.formula(paste(m_var, "~", paste(c(x_var, confounders), collapse = " + "))),
    data = dat_cc,
    family = gaussian()
  )
  
  fit_y_total <- glm(
    as.formula(paste(OUTCOME, "~", paste(c(x_var, confounders), collapse = " + "))),
    data = dat_cc,
    family = binomial()
  )
  
  fit_y_direct <- glm(
    as.formula(paste(OUTCOME, "~", paste(c(x_var, m_var, confounders), collapse = " + "))),
    data = dat_cc,
    family = binomial()
  )
  
  a_info  <- safe_coef(fit_m, x_term)
  b_info  <- safe_coef(fit_y_direct, m_var)
  c_info  <- safe_coef(fit_y_total, x_term)
  cp_info <- safe_coef(fit_y_direct, x_term)
  
  indirect_ab <- unname(a_info["estimate"] * b_info["estimate"])
  
  boot_tbl <- bootstrap_one_pair(
    dat = dat_cc,
    x_var = x_var,
    x_term = x_term,
    m_var = m_var,
    y_var = OUTCOME,
    confounders = confounders,
    n_boot = N_BOOT
  )
  
  bind_cols(
    tibble(
      analysis = analysis_label,
      subgroup_value = ifelse(is.null(subgroup_value), "all", subgroup_value),
      exposure_var = x_var,
      exposure_term = x_term,
      mediator = m_var,
      n_complete = nrow(dat_cc),
      a_estimate = unname(a_info["estimate"]),
      a_se = unname(a_info["se"]),
      a_p = unname(a_info["p"]),
      b_estimate = unname(b_info["estimate"]),
      b_se = unname(b_info["se"]),
      b_p = unname(b_info["p"]),
      total_effect_c = unname(c_info["estimate"]),
      total_se = unname(c_info["se"]),
      total_p = unname(c_info["p"]),
      direct_effect_cprime = unname(cp_info["estimate"]),
      direct_se = unname(cp_info["se"]),
      direct_p = unname(cp_info["p"]),
      indirect_effect_ab = indirect_ab,
      indirect_direction = ifelse(is.na(indirect_ab), NA_character_,
                                  ifelse(indirect_ab > 0, "positive", "negative")),
      note = NA_character_
    ),
    boot_tbl
  )
}

# --------------------------------------------------
# run
# --------------------------------------------------
df <- readRDS(DATA_PATH)

global_missing <- setdiff(unique(c(OUTCOME, CONFOUNDERS)), names(df))
if (length(global_missing) > 0) {
  stop("Global variables missing in data: ", paste(global_missing, collapse = ", "))
}

pairs <- read_csv(PAIR_CSV, show_col_types = FALSE)

req_cols <- c("analysis", "subgroup_label", "exposure_var", "exposure_term", "mediator")
miss <- setdiff(req_cols, names(pairs))
if (length(miss) > 0) {
  stop("Pair csv missing columns: ", paste(miss, collapse = ", "))
}

subgroup_value <- infer_subset_value(SUBGROUP_LABEL, df[[SEX_VAR]])

out_subdir <- file.path(OUT_DIR, JOB_LABEL)
dir.create(out_subdir, recursive = TRUE, showWarnings = FALSE)

res_list <- vector("list", nrow(pairs))

for (i in seq_len(nrow(pairs))) {
  cat("Running ", JOB_LABEL, " pair ", i, "/", nrow(pairs), ": ",
      pairs$exposure_term[i], " -> ", pairs$mediator[i], "\n", sep = "")
  
  res_list[[i]] <- run_one_pair(
    df = df,
    analysis_label = pairs$analysis[i],
    subgroup_value = subgroup_value,
    x_var = pairs$exposure_var[i],
    x_term = pairs$exposure_term[i],
    m_var = pairs$mediator[i]
  )
}

results <- bind_rows(res_list)

summary_tbl <- tibble(
  job_label = JOB_LABEL,
  subgroup_label = SUBGROUP_LABEL,
  subgroup_value = ifelse(is.null(subgroup_value), "all", subgroup_value),
  n_pairs = nrow(pairs),
  n_complete_results = sum(!is.na(results$indirect_effect_ab))
)

write_csv(pairs, file.path(out_subdir, paste0("pairs_", JOB_LABEL, ".csv")))
write_csv(results, file.path(out_subdir, paste0("formal_mediation_results_", JOB_LABEL, ".csv")))
write_csv(summary_tbl, file.path(out_subdir, paste0("summary_", JOB_LABEL, ".csv")))

cat("Done: ", JOB_LABEL, "\n", sep = "")
cat("Output: ", out_subdir, "\n", sep = "")