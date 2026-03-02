# analysis/table1_utils.R
# -----------------------
# Utility functions to build Table 1 as CSV without gtsummary
# --- formatting helpers ---
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

# --- test helpers ---
p_continuous_ttest <- function(x, g) {
  ok <- !is.na(x) & !is.na(g)
  x <- x[ok]; g <- g[ok]
  if (length(unique(g)) < 2) return(NA_real_)
  tryCatch(t.test(x ~ g)$p.value, error = function(e) NA_real_)
}

p_categorical_chi_fisher <- function(x, g) {
  ok <- !is.na(x) & !is.na(g)
  x <- x[ok]; g <- g[ok]
  if (length(unique(g)) < 2) return(NA_real_)
  x <- as.factor(x)
  tab <- table(x, g)
  if (nrow(tab) < 2) return(NA_real_)
  
  # expected counts from chi-square
  p_chi <- tryCatch(chisq.test(tab)$p.value, error = function(e) NA_real_)
  expct <- tryCatch(chisq.test(tab)$expected, error = function(e) NULL)
  
  # Fisher only for 2x2 with expected < 5
  if (!is.null(expct) && all(dim(tab) == c(2, 2)) && any(expct < 5)) {
    return(tryCatch(fisher.test(tab)$p.value, error = function(e) p_chi))
  }
  p_chi
}

fmt_p <- function(p) {
  if (length(p) == 0 || is.null(p) || is.na(p)) return("")
  p <- as.numeric(p)
  if (is.na(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

# --- collapse high-cardinality categorical variables ---
collapse_topk_other <- function(x, topk = 10) {
  # preserve NA; collapse rare levels
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
    if (is.character(df[[v]]) || is.factor(df[[v]]) || is.logical(df[[v]])) {
      df[[v]] <- collapse_topk_other(df[[v]], topk = topk)
    }
  }
  df
}

# --- core table builder (2-group stratification + overall) ---
make_table1_two_group <- function(df,
                                  vars,
                                  strata_var,
                                  strata_levels,
                                  strata_labels,
                                  add_missing_overall = FALSE,
                                  collapse_topk = NULL,
                                  topk = 10) {
  stopifnot(length(strata_levels) == 2, length(strata_labels) == 2)
  
  # Keep vars present and remove always-exclude handled outside
  vars <- intersect(vars, names(df))
  
  # Prepare output DF with only needed columns
  keep_cols <- intersect(c(strata_var, vars), names(df))
  df2 <- df[, keep_cols, drop = FALSE]
  
  # Coerce strata to factor with desired order/labels
  # --- robust mapping of strata variable to two labels ---
  g_raw <- df2[[strata_var]]
  
  # convert raw strata values to comparable "keys"
  if (is.logical(g_raw)) {
    # FALSE/TRUE -> "0"/"1"
    g_key <- as.character(as.integer(g_raw))
  } else {
    g_key <- as.character(g_raw)
  }
  
  level_keys <- as.character(strata_levels)
  
  g <- ifelse(g_key == level_keys[1], strata_labels[1],
              ifelse(g_key == level_keys[2], strata_labels[2], NA_character_))
  
  # store back (optional) and define group names
  df2[[strata_var]] <- g
  groups <- strata_labels
  
  # Convert character -> factor; logical -> factor (Yes/No)
  chr_vars <- names(df2)[sapply(df2, is.character)]
  if (length(chr_vars) > 0) df2[chr_vars] <- lapply(df2[chr_vars], as.factor)
  
  logi_vars <- names(df2)[sapply(df2, is.logical)]
  if (length(logi_vars) > 0) {
    df2[logi_vars] <- lapply(df2[logi_vars], function(x) factor(x, levels = c(FALSE, TRUE), labels = c("No", "Yes")))
  }
  
  # Optionally collapse high-cardinality categoricals (used in Appendix)
  if (!is.null(collapse_topk) && length(collapse_topk) > 0) {
    df2 <- collapse_df_categoricals(df2, collapse_topk, topk = topk)
  }
  
  g <- as.character(df2[[strata_var]])
  groups <- strata_labels  # exactly two labels in desired order
  
  denom_overall <- sum(!is.na(g))
  denom_g1 <- sum(g == groups[1], na.rm = TRUE)
  denom_g2 <- sum(g == groups[2], na.rm = TRUE)
  
  out_rows <- list()
  
  for (v in vars) {
    x <- df2[[v]]
    is_cont <- is.numeric(x) || is.integer(x)
    is_cat  <- is.factor(x) || is.character(x) || is.logical(x)
    
    miss_col <- if (add_missing_overall) fmt_missing_n_pct(x, denom_overall) else NULL
    
    if (is_cont) {
      p <- p_continuous_ttest(x, g)
      row <- data.frame(
        variable = v,
        level = "",
        Overall = fmt_mean_sd(x),
        group1 = fmt_mean_sd(x[g == groups[1]]),
        group2 = fmt_mean_sd(x[g == groups[2]]),
        p_value = fmt_p(p),
        stringsAsFactors = FALSE
      )
      if (add_missing_overall) row$missing_overall <- miss_col
      out_rows[[length(out_rows) + 1]] <- row
      
    } else if (is_cat) {
      x <- as.factor(x)
      p <- p_categorical_chi_fisher(x, g)
      
      header <- data.frame(
        variable = v,
        level = "",
        Overall = "",
        group1 = "",
        group2 = "",
        p_value = fmt_p(p),
        stringsAsFactors = FALSE
      )
      if (add_missing_overall) header$missing_overall <- miss_col
      out_rows[[length(out_rows) + 1]] <- header
      
      lvls <- levels(x)
      for (lv in lvls) {
        row <- data.frame(
          variable = "",
          level = paste0("  ", lv),
          Overall = fmt_n_pct(x == lv, denom_overall),
          group1 = fmt_n_pct(x[g == groups[1]] == lv, denom_g1),
          group2 = fmt_n_pct(x[g == groups[2]] == lv, denom_g2),
          p_value = "",
          stringsAsFactors = FALSE
        )
        if (add_missing_overall) row$missing_overall <- ""
        out_rows[[length(out_rows) + 1]] <- row
      }
      
    } else {
      # Unsupported type (e.g., Date) -> skip
      next
    }
  }
  
  tbl <- do.call(rbind, out_rows)
  row.names(tbl) <- NULL
  tbl <- as.data.frame(tbl, stringsAsFactors = FALSE)
  
  # Rename group columns to proper labels
  names(tbl)[names(tbl) == "group1"] <- strata_labels[1]
  names(tbl)[names(tbl) == "group2"] <- strata_labels[2]
  
  if (add_missing_overall) {
    names(tbl)[names(tbl) == "missing_overall"] <- "Missing n (%) (Overall)"
  }
  
  tbl
}