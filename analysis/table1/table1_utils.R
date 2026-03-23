# analysis/table1/table1_utils.R
# ------------------------------------------------------------
# Utility functions to build Table 1 as CSV (no gtsummary)
# Includes:
# - robust "numeric-like" coercion for factor/character
# - missingness column (Missing n (%) Overall) when requested
# - optional TopK+Other collapse for true categoricals only
# ------------------------------------------------------------

# -------------------------
# Formatting helpers
# -------------------------
fmt_mean_sd <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("")
  sprintf("%.2f (%.2f)", mean(x), sd(x))
}

fmt_n_pct <- function(logical_vec, denom_n) {
  n <- sum(logical_vec, na.rm = TRUE)
  if (is.na(denom_n) || denom_n == 0) return("")
  sprintf("%d (%.1f%%)", n, 100 * n / denom_n)
}

fmt_missing_n_pct <- function(x, denom_n) {
  n_miss <- sum(is.na(x))
  if (is.na(denom_n) || denom_n == 0) return("")
  sprintf("%d (%.1f%%)", n_miss, 100 * n_miss / denom_n)
}

fmt_p <- function(p) {
  if (length(p) == 0 || is.null(p) || is.na(p)) return("")
  p <- as.numeric(p)
  if (is.na(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

# -------------------------
# Tests
# -------------------------
p_continuous_ttest <- function(x, g) {
  ok <- !is.na(x) & !is.na(g)
  x <- x[ok]; g <- g[ok]
  if (length(unique(g)) < 2) return(NA_real_)
  tryCatch(t.test(x ~ g)$p.value, error = function(e) NA_real_)
}

p_categorical_chi_fisher <- function(x, g) {
  ok <- !is.na(x) & !is.na(g)
  x <- x[ok]; g <- g[ok]
  if (length(x) == 0 || length(unique(g)) < 2) return(NA_real_)
  x <- as.factor(x)
  tab <- table(x, g)
  if (nrow(tab) < 2) return(NA_real_)
  
  p_chi <- tryCatch(chisq.test(tab)$p.value, error = function(e) NA_real_)
  expct <- tryCatch(chisq.test(tab)$expected, error = function(e) NULL)
  
  # Fisher only for 2x2 with expected < 5
  if (!is.null(expct) && all(dim(tab) == c(2, 2)) && any(expct < 5)) {
    return(tryCatch(fisher.test(tab)$p.value, error = function(e) p_chi))
  }
  p_chi
}

# -------------------------
# Robust numeric coercion
# -------------------------
is_numeric_like <- function(x_chr, min_success_rate = 0.95) {
  # x_chr: character vector
  x_chr_trim <- trimws(x_chr)
  x_chr_trim[x_chr_trim %in% c("", "NA", "NaN", "nan")] <- NA_character_
  
  # keep non-missing only
  v <- x_chr_trim[!is.na(x_chr_trim)]
  if (length(v) == 0) return(FALSE)
  
  suppressWarnings(num <- as.numeric(v))
  success_rate <- mean(!is.na(num))
  isTRUE(success_rate >= min_success_rate)
}

coerce_numeric_if_possible <- function(x, min_success_rate = 0.95) {
  # If factor/character but essentially numeric -> convert to numeric
  if (is.factor(x)) {
    x_chr <- as.character(x)
    if (is_numeric_like(x_chr, min_success_rate)) {
      suppressWarnings(x_num <- as.numeric(trimws(x_chr)))
      return(x_num)
    }
    return(x) # keep as factor
  }
  if (is.character(x)) {
    if (is_numeric_like(x, min_success_rate)) {
      suppressWarnings(x_num <- as.numeric(trimws(x)))
      return(x_num)
    }
    return(as.factor(x))
  }
  x
}

# -------------------------
# Collapse high-cardinality categoricals (TopK + Other)
# -------------------------
collapse_topk_other <- function(x, topk = 10) {
  # preserve NA; collapse rare levels; expects categorical-like input
  x0 <- x
  x <- as.character(x)
  x_non_na <- x[!is.na(x)]
  if (length(x_non_na) == 0) return(as.factor(x0))
  
  tb <- sort(table(x_non_na), decreasing = TRUE)
  if (length(tb) <= topk) return(as.factor(x0))
  
  keep <- names(tb)[1:topk]
  x[!is.na(x) & !(x %in% keep)] <- "Other"
  factor(x, levels = c(keep, "Other"))
}

collapse_df_categoricals <- function(df, vars, topk = 10) {
  for (v in vars) {
    if (!v %in% names(df)) next
    x <- df[[v]]
    
    # IMPORTANT: try numeric coercion first; only collapse if still categorical
    x2 <- coerce_numeric_if_possible(x)
    if (is.numeric(x2) || is.integer(x2)) {
      df[[v]] <- x2
      next
    }
    
    if (is.character(x2) || is.factor(x2) || is.logical(x2)) {
      df[[v]] <- collapse_topk_other(x2, topk = topk)
    } else {
      df[[v]] <- x2
    }
  }
  df
}

# -------------------------
# Core table builder: 2-group stratification + overall
# -------------------------
make_table1_two_group <- function(df,
                                  vars,
                                  strata_var,
                                  strata_levels,
                                  strata_labels,
                                  add_missing_overall = FALSE,
                                  collapse_topk = NULL,
                                  topk = 10,
                                  numeric_like_success = 0.95) {
  
  stopifnot(length(strata_levels) == 2, length(strata_labels) == 2)
  
  # Keep vars present
  vars <- intersect(vars, names(df))
  
  # Prepare df with needed columns only
  keep_cols <- intersect(c(strata_var, vars), names(df))
  df2 <- df[, keep_cols, drop = FALSE]
  
  # Robust mapping of strata into exactly 2 labels
  g_raw <- df2[[strata_var]]
  if (is.logical(g_raw)) {
    g_key <- as.character(as.integer(g_raw)) # FALSE/TRUE -> "0"/"1"
  } else {
    g_key <- as.character(g_raw)
  }
  
  level_keys <- as.character(strata_levels)
  g <- ifelse(g_key == level_keys[1], strata_labels[1],
              ifelse(g_key == level_keys[2], strata_labels[2], NA_character_))
  
  df2[[strata_var]] <- g
  groups <- strata_labels
  
  # Convert character -> factor; logical -> factor (Yes/No)
  chr_vars <- names(df2)[sapply(df2, is.character)]
  if (length(chr_vars) > 0) df2[chr_vars] <- lapply(df2[chr_vars], as.factor)
  
  logi_vars <- names(df2)[sapply(df2, is.logical)]
  if (length(logi_vars) > 0) {
    df2[logi_vars] <- lapply(df2[logi_vars], function(x) {
      factor(x, levels = c(FALSE, TRUE), labels = c("No", "Yes"))
    })
  }
  
  # Optional: collapse high-cardinality categoricals (Appendix)
  if (!is.null(collapse_topk) && length(collapse_topk) > 0) {
    df2 <- collapse_df_categoricals(df2, collapse_topk, topk = topk)
  }
  
  # Denominators (exclude missing group labels)
  g2 <- df2[[strata_var]]
  denom_overall <- sum(!is.na(g2))
  denom_g1 <- sum(g2 == groups[1], na.rm = TRUE)
  denom_g2 <- sum(g2 == groups[2], na.rm = TRUE)
  
  out_rows <- list()
  
  for (v in vars) {
    x_raw <- df2[[v]]
    
    # Key: coerce numeric-like factor/character to numeric BEFORE type decision
    x <- coerce_numeric_if_possible(x_raw, min_success_rate = numeric_like_success)
    
    is_cont <- is.numeric(x) || is.integer(x)
    is_cat  <- is.factor(x) || is.character(x) || is.logical(x)
    
    miss_col <- if (add_missing_overall) fmt_missing_n_pct(x, denom_overall) else NULL
    
    if (is_cont) {
      p <- p_continuous_ttest(x, g2)
      row <- data.frame(
        variable = v,
        level = "",
        Overall = fmt_mean_sd(x),
        group1 = fmt_mean_sd(x[g2 == groups[1]]),
        group2 = fmt_mean_sd(x[g2 == groups[2]]),
        p_value = fmt_p(p),
        stringsAsFactors = FALSE
      )
      if (add_missing_overall) row[["Missing n (%) (Overall)"]] <- miss_col
      out_rows[[length(out_rows) + 1]] <- row
      
    } else if (is_cat) {
      # Make sure categorical is factor for stable level order
      x_fac <- as.factor(x)
      
      p <- p_categorical_chi_fisher(x_fac, g2)
      header <- data.frame(
        variable = v,
        level = "",
        Overall = "",
        group1 = "",
        group2 = "",
        p_value = fmt_p(p),
        stringsAsFactors = FALSE
      )
      if (add_missing_overall) header[["Missing n (%) (Overall)"]] <- miss_col
      out_rows[[length(out_rows) + 1]] <- header
      
      lvls <- levels(x_fac)
      for (lv in lvls) {
        row <- data.frame(
          variable = "",
          level = paste0("  ", lv),
          Overall = fmt_n_pct(x_fac == lv, denom_overall),
          group1 = fmt_n_pct(x_fac[g2 == groups[1]] == lv, denom_g1),
          group2 = fmt_n_pct(x_fac[g2 == groups[2]] == lv, denom_g2),
          p_value = "",
          stringsAsFactors = FALSE
        )
        if (add_missing_overall) row[["Missing n (%) (Overall)"]] <- ""
        out_rows[[length(out_rows) + 1]] <- row
      }
      
    } else {
      # Unsupported type -> skip
      next
    }
  }
  
  tbl <- do.call(rbind, out_rows)
  rownames(tbl) <- NULL
  tbl <- as.data.frame(tbl, stringsAsFactors = FALSE)
  
  # Rename group columns
  names(tbl)[names(tbl) == "group1"] <- strata_labels[1]
  names(tbl)[names(tbl) == "group2"] <- strata_labels[2]
  
  tbl
}