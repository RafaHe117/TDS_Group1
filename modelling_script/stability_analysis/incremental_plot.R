library(dplyr)
library(pROC)
library(caret)
library(ggplot2)
library(Matrix)
library(stringr)

setwd("/rds/general/project/hda_25-26/live/TDS/TDS_Group1")

# Reuse helper functions from 
# sex_difference_visual.R: assign_domain(), make_pretty_names()
# stability_sex.R: prepare_xy()

################################################################################
# 0. Maps a dummy term back to its parent variable
#
# This function takes a dummified variable name (e.g., "IMD_quintile3") and a list 
# of original exposure names, returning the matching parent name (e.g., "IMD_quintile"). 
# It sorts the exposures by length descending prior to matching to prevent false 
# positive substring matching errors.
################################################################################

get_parent_var <- function(variable_prefix, exposures) {
  
  exposures <- exposures[order(-nchar(exposures))] 
  
  for (exposure in exposures) {
    if (startsWith(variable_prefix, exposure)) return(exposure)
  }
  
return(variable_prefix)
}

################################################################################
# 1. Calculate Subset AUC
#
# This function extracts the specified variable columns from the design matrices, 
# fits an unpenalised logistic regression model, and computes the Test AUC 
# alongside its DeLong standard error to measure predictive accuracy.
################################################################################
cal_auc <- function(y_tr, X_tr, y_te, X_te, stable_vars) {
  
  # Subset and clean names
  df_train <- cbind(y = y_tr, as.data.frame(X_tr[, stable_vars, drop = FALSE]))
  df_test  <- as.data.frame(X_te[, stable_vars, drop = FALSE])
  
  colnames(df_train) <- make.names(colnames(df_train), unique = TRUE)
  colnames(df_test)  <- make.names(colnames(df_test), unique = TRUE)
  
  # Fit and Predict
  fit <- suppressWarnings(glm(y ~ ., data = df_train, family = binomial()))
  preds_prob <- suppressWarnings(predict(fit, newdata = df_test, type = "response"))
  
  # Calculate AUC and DeLong Confidence Interval
  roc_obj <- roc(y_te, preds_prob, quiet = TRUE)
  auc_val <- as.numeric(auc(roc_obj))
  
  # ci.auc returns c(lower_ci, median_auc, upper_ci). Derive SE from the 95% CI width.
  ci_val <- as.numeric(ci.auc(roc_obj, method = "delong")) 
  se_val <- (ci_val[3] - ci_val[1]) / (2 * qnorm(0.975))
  
  return(data.frame(Mean_AUC = auc_val, SE = se_val))
}

################################################################################
# 2. Compute Incremental AUC
#
# This function iteratively adds exposome variables to a baseline confounder model. 
# It handles both domain-grouped and overall sequential addition, adds them sequentially 
# based on their stability proportion, and applies the 1-SE rule to identify the optimal cutoff.
################################################################################
cal_incremental_performance <- function(train_df, test_df, confounders, selected_csv_path, 
                                        exposure_csv_path, extra_exposures = NULL, group_by_domain) {
  
  exposure_df <- read.csv(exposure_csv_path, stringsAsFactors = FALSE)
  exposures <- exposure_df$exposure
  
  # Append any cohort-specific exposures (for the SBP cohort)
  if (!is.null(extra_exposures)) {
    exposures <- unique(c(exposures, extra_exposures))
  }
  
  exposure_df <- data.frame(exposure = exposures, stringsAsFactors = FALSE)
  
  # Load and assign domains to the selected variables
  sel_df <- read.csv(selected_csv_path, stringsAsFactors = FALSE) |>
    filter(selected == TRUE) |>
    mutate(Domain = assign_domain(term))

  # Apply mapping and aggregate by parent variable
  sel_df$parent_var <- sapply(sel_df$term, get_parent_var, exposures = exposures)
  
  sel_df <- sel_df |> 
    group_by(parent_var, Domain) |> 
    # Use the max selection proportion among the dummy levels to determine the parent variable's rank
    summarise(selection_proportion = max(selection_proportion), .groups = "drop") |> 
    rename(term = parent_var)
  
  # Apply sorting logic based on the selected method (domain-grouped or overall)
  if (group_by_domain) {
    sel_df <- sel_df |> arrange(Domain, desc(selection_proportion))
  } else {
    sel_df <- sel_df |> arrange(desc(selection_proportion))
  }
  
  # Prepare full matrices
  preds <- unique(c(confounders, exposure_df$exposure))
  
  xy <- prepare_xy(train_df, test_df, preds)
  X_tr <- xy$X_tr_mat
  X_te <- xy$X_te_mat
  
  y_tr <- as.integer(train_df$cvd_incident)
  y_te <- as.integer(test_df$cvd_incident)
  
  base_cols <- unlist(lapply(confounders, function(v) grep(paste0("^", v), colnames(X_tr), value = TRUE)))
  all_results <- data.frame()
  
  # Method A: Grouped by Domain
  if (group_by_domain) {
    domains <- unique(sel_df$Domain)
    
    for (domain in domains) {
      dom_features <- sel_df |> filter(Domain == domain)
      cumulative_terms <- c()
      dom_results <- data.frame()
      
      for (i in 1:nrow(dom_features)) {
        current_term <- dom_features$term[i]
        cumulative_terms <- c(cumulative_terms, current_term)
        
        term_cols <- unlist(lapply(cumulative_terms, function(v) grep(paste0("^", v), colnames(X_tr), value = TRUE)))
        active_cols <- c(base_cols, term_cols)
        eval_res <- cal_auc(y_tr, X_tr, y_te, X_te, active_cols)
        
        pretty_name <- make_pretty_names(current_term)
        label_str <- ifelse(i == 1, pretty_name, paste0("+ \n", pretty_name))
        
        dom_results <- rbind(dom_results, data.frame(
          Method = "Grouped", Domain = domain, Step = i, Feature = current_term, 
          Label = label_str, Mean_AUC = eval_res$Mean_AUC, SE = eval_res$SE
        ))
      }
      
      # Apply 1-SE Rule per domain
      max_auc_idx <- which.max(dom_results$Mean_AUC)
      threshold <- dom_results$Mean_AUC[max_auc_idx] - dom_results$SE[max_auc_idx]
      cutoff_step <- min(dom_results$Step[dom_results$Mean_AUC >= threshold])
      
      dom_results$Status <- ifelse(dom_results$Step <= cutoff_step, "Selected", "Unselected")
      dom_results$Cutoff <- cutoff_step
      all_results <- rbind(all_results, dom_results)
    }
    
    # Method B: Overall Sequential Addition
  } else {
    cumulative_terms <- c()
    
    for (i in 1:nrow(sel_df)) {
      current_term <- sel_df$term[i]
      cumulative_terms <- c(cumulative_terms, current_term)
      
      term_cols <- unlist(lapply(cumulative_terms, function(v) grep(paste0("^", v), colnames(X_tr), value = TRUE)))
      active_cols <- c(base_cols, term_cols)
      eval_res <- cal_auc(y_tr, X_tr, y_te, X_te, active_cols)
      
      pretty_name <- make_pretty_names(current_term)
      label_str <- ifelse(i == 1, pretty_name, paste0("+ \n", pretty_name))
      
      all_results <- rbind(all_results, data.frame(
        Method = "Overall", Domain = sel_df$Domain[i], Step = i, Feature = current_term, 
        Label = label_str, Mean_AUC = eval_res$Mean_AUC, SE = eval_res$SE
      ))
    }
    
    # Apply 1-SE Rule overall
    max_auc_idx <- which.max(all_results$Mean_AUC)
    threshold <- all_results$Mean_AUC[max_auc_idx] - all_results$SE[max_auc_idx]
    cutoff_step <- min(all_results$Step[all_results$Mean_AUC >= threshold])
    
    all_results$Status <- ifelse(all_results$Step <= cutoff_step, "Selected", "Unselected")
    all_results$Cutoff <- cutoff_step
  }
  
  return(all_results)
}

################################################################################
# 3. Plotting
#
# This function constructs a faceted scatter-and-error-bar plot to visualise the 
# incremental AUC changes.
################################################################################
create_incremental_plot <- function(results_df, cohort_title, group_by_domain) {
  
  # Lock feature ordering
  results_df$Label <- factor(results_df$Label, levels = unique(results_df$Label))
  
  # Map domain_fill to the Point_Color logic
  domain_colors <- data.frame(
    Domain = c("Socioeconomic", "Lifestyle / Behaviour", "Clinical / Family", "Environmental", "Other"),
    Dom_Color = c("#D4B483", "#48A9A6", "#C1666B", "#4281A4", "#6A5D7B")
  )
  
  df_plot <- merge(results_df, domain_colors, by = "Domain", all.x = TRUE)
  df_plot$Point_Color <- ifelse(df_plot$Status == "Selected", df_plot$Dom_Color, "#B0B0B0")
  
  # Create the base plot
  p <- ggplot(df_plot, aes(x = Label, y = Mean_AUC, color = Point_Color)) +
    geom_errorbar(aes(ymin = Mean_AUC - SE, ymax = Mean_AUC + SE), width = 0.2, linewidth = 0.8) +
    geom_point(size = 3.5) +
    geom_vline(aes(xintercept = Cutoff + 0.5), linetype = "dashed", color = "gray50", linewidth = 0.8) +
    scale_color_identity() +
    theme_bw(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
      axis.title.x = element_blank(),
      legend.position = "none",
      plot.title = element_text(size = 16, face = "bold"),
      panel.grid.minor = element_blank()
    ) 
  
  # Apply faceting if group by domain is True
  if (group_by_domain) {
    p <- p + facet_wrap(~ Domain, scales = "free_x", ncol = 1) + 
      theme(
        strip.background = element_rect(fill = "NA", color = "NA"),
        strip.text = element_text(size = 14, face = "bold", hjust = 0)
      ) +
      labs(y = "Test AUC", title = paste("Incremental Domain AUC:", cohort_title))
  } else {
    p <- p + labs(y = "Test AUC", title = paste("Overall Incremental AUC:", cohort_title))
  }
  
  return(p)
}

################################################################################
# 4. Main Execution and AUC Summary
################################################################################
base_out_dir <- "modelling_script/stability_analysis"

train_full <- readRDS("split_imputed_data/ukb_G1_train_imputed.rds")
test_full  <- readRDS("split_imputed_data/ukb_G1_test_imputed.rds")

calc_plot_height <- function(results_df) {
  n_domains <- length(unique(results_df$Domain))
  return(2 + (n_domains * 2.5))
}

# Define configuration parameters for each cohort to avoid repeated code
cohort_configs <- list(
  
  "Full" = list(
    train = train_full, 
    test = test_full,
    confounders = c("age", "sex", "ethnicity_5cat"),
    sel_csv = file.path(base_out_dir, "subsample_lasso", "selection_proportions_penalized.csv"),
    exp_csv = file.path(base_out_dir, "subsample_lasso", "exposure_list.csv"),
    out_dir = file.path(base_out_dir, "subsample_lasso"),
    suffix = "full"
  ),
  
  "Female" = list(
    train = train_full[train_full$sex == "Female", ], 
    test = test_full[test_full$sex == "Female", ],
    confounders = c("age", "ethnicity_5cat"),
    sel_csv = file.path(base_out_dir, "subsample_lasso_sex", "female", "selection_proportions_penalized_female.csv"),
    exp_csv = file.path(base_out_dir, "subsample_lasso", "exposure_list.csv"),
    out_dir = file.path(base_out_dir, "subsample_lasso_sex"),
    suffix = "female"
  ),
  
  "Male" = list(
    train = train_full[train_full$sex == "Male", ], 
    test = test_full[test_full$sex == "Male", ],
    confounders = c("age", "ethnicity_5cat"),
    sel_csv = file.path(base_out_dir, "subsample_lasso_sex", "male", "selection_proportions_penalized_male.csv"),
    exp_csv = file.path(base_out_dir, "subsample_lasso", "exposure_list.csv"),
    out_dir = file.path(base_out_dir, "subsample_lasso_sex"),
    suffix = "male"
  ),
  
  "SBP" = list(
    train = train_full, 
    test = test_full,
    confounders = c("age", "sex", "ethnicity_5cat"),
    sel_csv = file.path(base_out_dir, "subsample_lasso_SBP", "selection_proportions_penalized.csv"),
    exp_csv = file.path(base_out_dir, "subsample_lasso", "exposure_list.csv"),
    out_dir = file.path(base_out_dir, "subsample_lasso_SBP"),
    suffix = "sbp",
    extra_exposures = c("sbp_mean", "dbp_mean")
  ),
  
  "BP Confounder" = list(
    train = train_full, 
    test = test_full,
    confounders = c("age", "sex", "ethnicity_5cat", "sbp_mean", "dbp_mean"),
    sel_csv = file.path(base_out_dir, "subsample_lasso_bp_confounder", "selection_proportions_penalized.csv"),
    exp_csv = file.path(base_out_dir, "subsample_lasso", "exposure_list.csv"),
    out_dir = file.path(base_out_dir, "subsample_lasso_bp_confounder"),
    suffix = "bp_confounder"
  )
)

# Initialise data structures to collect summary metrics during iteration
summary_rows <- list()
auc_trajectories <- list()

# Iterate through each cohort configuration
for (cohort_name in names(cohort_configs)) {
  cfg <- cohort_configs[[cohort_name]]
  
  cat("Processing:", cohort_name, "Cohort...\n")
  
  # 1. Base AUC Calculation
  fit_base <- suppressWarnings(glm(cvd_incident ~ ., data = cbind(cvd_incident = as.integer(cfg$train$cvd_incident), cfg$train[, cfg$confounders, drop = FALSE]), family = binomial()))
  pred_base <- suppressWarnings(predict(fit_base, newdata = cfg$test, type = "response"))
  base_auc_val <- as.numeric(auc(roc(as.integer(cfg$test$cvd_incident), pred_base, quiet = TRUE)))
  
  # 2. Grouped by Domain Analysis & Plotting
  res_grouped <- cal_incremental_performance(
    train_df = cfg$train, test_df = cfg$test, confounders = cfg$confounders,
    selected_csv_path = cfg$sel_csv, exposure_csv_path = cfg$exp_csv, extra_exposures = cfg$extra_exposures, group_by_domain = TRUE
  )
  plot_grouped <- create_incremental_plot(res_grouped, paste(cohort_name, "Cohort"), group_by_domain = TRUE)
  ggsave(file.path(cfg$out_dir, paste0("domain_incremental_", cfg$suffix, ".png")), plot_grouped, width = 10, height = calc_plot_height(res_grouped), dpi = 300)
  
  write.csv(res_grouped, file.path(cfg$out_dir, paste0("domain_incremental_data_", cfg$suffix, ".csv")), row.names = FALSE)
  
  # 3. Overall Analysis & Plotting
  res_overall <- cal_incremental_performance(
    train_df = cfg$train, test_df = cfg$test, confounders = cfg$confounders,
    selected_csv_path = cfg$sel_csv, exposure_csv_path = cfg$exp_csv, extra_exposures = cfg$extra_exposures, group_by_domain = FALSE
  )
  plot_overall <- create_incremental_plot(res_overall, paste(cohort_name, "Cohort"), group_by_domain = FALSE)
  ggsave(file.path(cfg$out_dir, paste0("overall_incremental_", cfg$suffix, ".png")), plot_overall, width = 12, height = 6, dpi = 300)
  
  write.csv(res_overall, file.path(cfg$out_dir, paste0("overall_incremental_data_", cfg$suffix, ".csv")), row.names = FALSE)
  
  # 4. Extract Metrics for Summary
  summary_rows[[cohort_name]] <- data.frame(
    Cohort = cohort_name,
    Base_AUC = round(base_auc_val, 6),
    Max_Test_AUC = round(res_overall$Mean_AUC[which.max(res_overall$Mean_AUC)], 6),
    Optimal_1SE_AUC = round(res_overall$Mean_AUC[res_overall$Cutoff[1]], 6)
  )
  
  # Store the step-by-step AUC trajectory for dynamic padding later
  auc_trajectories[[cohort_name]] <- res_overall$Mean_AUC
}

# Combine the core metrics into a single dataframe
auc_summary <- do.call(rbind, summary_rows)
rownames(auc_summary) <- NULL

# Determine the maximum number of incremental steps across all cohorts
max_vars <- max(sapply(auc_trajectories, length))

# Dynamically pad and append step-by-step AUC columns
for (i in 1:max_vars) {
  col_name <- paste0("Step_", i, "_AUC")
  
  # Extract the ith value for each cohort, returning NA if the cohort has fewer variables
  step_values <- sapply(names(cohort_configs), function(cohort) {
    if (i <= length(auc_trajectories[[cohort]])) {
      return(round(auc_trajectories[[cohort]][i], 6))
    } else {
      return(NA)
    }
  })
  
  auc_summary[[col_name]] <- step_values
}

print(auc_summary)

write.csv(auc_summary, file.path(base_out_dir, "incremental_test_auc_summary.csv"), row.names = FALSE)