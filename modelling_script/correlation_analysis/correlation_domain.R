library(dplyr)
library(pheatmap)
library(grid)
library(parallel)

############################################################
# 0. Paths
############################################################
in_dir  <- "/rds/general/user/rh725/projects/hda_25-26/live/TDS/TDS_Group1"
out_dir <- file.path(in_dir, "modelling_script")
corr_dir <- file.path(out_dir, "correlation_analysis")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(corr_dir, showWarnings = FALSE, recursive = TRUE)

############################################################
# 1. Load data
############################################################
UK1_1 <- readRDS(file.path(in_dir, "ukb_G1_imputed_final.rds"))

############################################################
# 2. Variable mapping
############################################################
predictor_map <- data.frame(
  variable = c(
    "fruit_intake_fresh","tea_intake","coffee_intake",
    "oily_fish_3cat","salt_3cat","alcohol_combined",
    "redwine_group","sf_score",
    "tv_group","total_met_group",
    "work_hours_week_clean","smoking_3cat",
    "sleep_quality_score",
    "neuroticism_score","mh_satis_mean_score",
    "mood_disorder","stress_count_2yr","stress_group_2yr",
    "no2_2010","pm10_2010","pm2_5_2010","noise_24h",
    "greenspace_pct_1000m","temp_average","urban_rural_2cat",
    "diabetes_bin","famhx_cvd",
    "ICD_Diabetes","ICD_Lipidemia","ICD_Hypertension",
    "ICD_AF","ICD_CKD","ICD_Mental_Health","ICD_Migraine",
    "ICD_Atopy","ICD_Autoimmune","medication_category",
    "sbp_mean","dbp_mean","bmi","creatinine_1","birth_weight_clean",
    "education_2","employment_2cat",
    "living_with_partner","imd_final",
    "imd_source","imd_quintile"
  ),
  field = c(
    rep("Lifestyle_Diet", 13),
    rep("Mental_Health", 5),
    rep("Environment", 7),
    rep("Medical_History", 12),
    rep("Physical_Socio", 11)
  ),
  stringsAsFactors = FALSE
)

biomarker_map <- data.frame(
  variable = c(
    "wbc_count","rbc_count","hemoglobin","mcv","platelet_count","mpv",
    "lymphocyte_count","monocyte_count","neutrophil_count",
    "eosinophil_count","basophil_count","nrbc_count",
    
    "creatinine","cystatin_c","urea","urate","microalbumin","creatinine_in_urine",
    "albumin","alkaline_phosphatase",
    "alanine_aminotransferase","aspartate_aminotransferase",
    "direct_bilirubin","total_bilirubin","total_blood_protein",
    
    "cholesterol","hdl_cholesterol","ldl_cholesterol",
    "lipoprotein_a","apolipoprotein_a","apolipoprotein_b",
    "total_triglyceride","glucose","hba1c","igf_1",
    "crp","blood_vitamin_d","calcium","phosphate"
  ),
  
  field = "Biomarker",
  
  biomarker_subgroup = c(
    rep("Haematology", 12),
    rep("Renal_Liver_Protein", 13),
    rep("Lipid_Metabolic_Inflammation", 14)
  ),
  
  stringsAsFactors = FALSE
)

extra_map <- data.frame(
  variable = c(
    "shift_work","mixed_shift","job_shift_work_clean","ever_nights","work_hours_unified_cat",
    "household_size_clean","living_alone",
    "self_health_bin","gp_anxdep","hrt"
  ),
  field = c(
    rep("Lifestyle_Diet", 5),
    rep("Physical_Socio", 2),
    rep("Medical_History", 3)
  ),
  stringsAsFactors = FALSE
)

if (!"biomarker_subgroup" %in% names(predictor_map)) {
  predictor_map$biomarker_subgroup <- NA_character_
}
if (!"biomarker_subgroup" %in% names(extra_map)) {
  extra_map$biomarker_subgroup <- NA_character_
}

all_map <- bind_rows(predictor_map, biomarker_map, extra_map) %>%
  distinct(variable, .keep_all = TRUE)

available_vars <- intersect(all_map$variable, names(UK1_1))
all_map <- all_map %>% filter(variable %in% available_vars)

data_all <- UK1_1 %>% select(all_of(all_map$variable))

message("Number of selected variables found in dataset: ", ncol(data_all))

############################################################
# 3. Mixed-correlation helper functions
############################################################
correlation_ratio <- function(categories, measurements) {
  measurements <- suppressWarnings(as.numeric(as.character(measurements)))
  valid_idx <- !is.na(measurements) & !is.na(categories)
  categories <- categories[valid_idx]
  measurements <- measurements[valid_idx]
  
  if (length(measurements) == 0) return(NA_real_)
  
  categories <- as.factor(categories)
  cat_levels <- levels(categories)
  if (length(cat_levels) < 2) return(NA_real_)
  
  y_avg <- sapply(cat_levels, function(cat) {
    vals <- measurements[categories == cat]
    if (length(vals) == 0) return(NA_real_)
    mean(vals, na.rm = TRUE)
  })
  
  n_array <- sapply(cat_levels, function(cat) sum(categories == cat))
  y_total_avg <- sum(y_avg * n_array, na.rm = TRUE) / sum(n_array, na.rm = TRUE)
  numerator <- sum(n_array * (y_avg - y_total_avg)^2, na.rm = TRUE)
  denominator <- sum((measurements - y_total_avg)^2, na.rm = TRUE)
  
  if (is.na(denominator) || denominator == 0) return(NA_real_)
  
  sqrt(numerator / denominator)
}

cramers_v <- function(x, y) {
  x <- as.character(x)
  y <- as.character(y)
  x[is.na(x)] <- "MISSING"
  y[is.na(y)] <- "MISSING"
  
  tbl <- table(x, y)
  
  if (sum(tbl) == 0) return(NA_real_)
  if (nrow(tbl) < 2 || ncol(tbl) < 2) return(NA_real_)
  
  chi2 <- suppressWarnings(chisq.test(tbl, correct = FALSE)$statistic)
  n <- sum(tbl)
  r <- nrow(tbl)
  k <- ncol(tbl)
  
  if (is.na(chi2) || n == 0 || min(r - 1, k - 1) == 0) return(NA_real_)
  
  as.numeric(sqrt(chi2 / (n * min(r - 1, k - 1))))
}

is_numeric_like <- function(x) {
  if (is.numeric(x) || is.integer(x)) return(TRUE)
  
  x2 <- suppressWarnings(as.numeric(as.character(x)))
  non_missing_original <- sum(!is.na(x))
  non_missing_converted <- sum(!is.na(x2))
  
  if (non_missing_original == 0) return(FALSE)
  
  (non_missing_converted / non_missing_original) >= 0.9
}

############################################################
# 4. Continuous-variable detector for heatmaps
############################################################
is_continuous_for_heatmap <- function(x, varname = NULL, min_unique = 8) {
  if (!is_numeric_like(x)) return(FALSE)
  
  x_num <- suppressWarnings(as.numeric(as.character(x)))
  x_num <- x_num[!is.na(x_num)]
  
  if (length(x_num) <= 1) return(FALSE)
  if (sd(x_num) == 0) return(FALSE)
  
  n_unique <- length(unique(x_num))
  if (n_unique < min_unique) return(FALSE)
  
  if (!is.null(varname)) {
    if (grepl("^ICD_", varname)) return(FALSE)
    if (grepl("_bin$", varname)) return(FALSE)
    if (grepl("_2cat$", varname)) return(FALSE)
    if (grepl("_3cat$", varname)) return(FALSE)
    if (grepl("_group$", varname)) return(FALSE)
    if (grepl("_group_", varname)) return(FALSE)
    if (grepl("_quintile$", varname)) return(FALSE)
    if (grepl("_category$", varname)) return(FALSE)
    if (grepl("_source$", varname)) return(FALSE)
    
    explicit_exclusions <- c(
      "education_2",
      "employment_2cat",
      "living_with_partner",
      "shift_work",
      "mixed_shift",
      "job_shift_work_clean",
      "ever_nights",
      "work_hours_unified_cat",
      "living_alone",
      "self_health_bin",
      "gp_anxdep",
      "hrt",
      "famhx_cvd",
      "mood_disorder",
      "medication_category",
      "urban_rural_2cat",
      "imd_quintile"
    )
    
    if (varname %in% explicit_exclusions) return(FALSE)
  }
  
  TRUE
}

############################################################
# 5. Parallel mixed correlation with progress bar
############################################################
mixed_corr_parallel <- function(df, categorical_cols = NULL, ncores = 8, chunk_size = 2) {
  
  if (is.null(categorical_cols)) {
    numeric_flags <- sapply(df, is_numeric_like)
    categorical_cols <- names(df)[!numeric_flags]
  }
  
  cols <- colnames(df)
  n <- length(cols)
  
  calc_one_pair <- function(col1, col2, df, categorical_cols) {
    if (col1 == col2) {
      return(1)
    }
    
    if (col1 %in% categorical_cols && col2 %in% categorical_cols) {
      return(cramers_v(df[[col1]], df[[col2]]))
    }
    
    if (col1 %in% categorical_cols && !(col2 %in% categorical_cols)) {
      return(correlation_ratio(df[[col1]], df[[col2]]))
    }
    
    if (!(col1 %in% categorical_cols) && col2 %in% categorical_cols) {
      return(correlation_ratio(df[[col2]], df[[col1]]))
    }
    
    x <- suppressWarnings(as.numeric(as.character(df[[col1]])))
    y <- suppressWarnings(as.numeric(as.character(df[[col2]])))
    
    if (all(is.na(x)) || all(is.na(y))) {
      return(NA_real_)
    } else {
      return(cor(x, y, use = "pairwise.complete.obs"))
    }
  }
  
  row_worker <- function(i_vec, df, cols, categorical_cols) {
    out <- lapply(i_vec, function(i) {
      sapply(seq_along(cols), function(j) {
        calc_one_pair(cols[i], cols[j], df, categorical_cols)
      })
    })
    names(out) <- cols[i_vec]
    out
  }
  
  ncores <- max(1, min(ncores, detectCores()))
  message("Starting mixed correlation using ", ncores, " cores...")
  
  row_index <- seq_len(n)
  chunks <- split(row_index, ceiling(seq_along(row_index) / chunk_size))
  
  cl <- makeCluster(ncores)
  on.exit(stopCluster(cl), add = TRUE)
  
  clusterExport(
    cl,
    varlist = c(
      "df", "cols", "categorical_cols",
      "correlation_ratio", "cramers_v", "is_numeric_like",
      "row_worker"
    ),
    envir = environment()
  )
  
  pb <- txtProgressBar(min = 0, max = length(chunks), style = 3)
  
  result_list <- vector("list", length(chunks))
  
  for (k in seq_along(chunks)) {
    result_list[[k]] <- parLapply(
      cl,
      list(chunks[[k]]),
      row_worker,
      df = df,
      cols = cols,
      categorical_cols = categorical_cols
    )[[1]]
    
    setTxtProgressBar(pb, k)
  }
  
  close(pb)
  
  row_results <- unlist(result_list, recursive = FALSE)
  row_results <- row_results[cols]
  
  corr <- do.call(rbind, row_results)
  corr <- as.matrix(corr)
  rownames(corr) <- cols
  colnames(corr) <- cols
  
  as.data.frame(corr)
}

############################################################
# 6. Run mixed correlation
############################################################
ncores <- as.integer(Sys.getenv("MY_NCORES", "16"))
chunk_size <- as.integer(Sys.getenv("MY_CHUNK_SIZE", "2"))

message("MY_NCORES from environment: ", ncores)
message("MY_CHUNK_SIZE from environment: ", chunk_size)

numeric_flags_all <- sapply(data_all, is_numeric_like)
categorical_cols <- names(data_all)[!numeric_flags_all]

message("Categorical variables detected: ", length(categorical_cols))
message("Numeric-like variables detected: ", sum(numeric_flags_all))

correlation_matrix_all <- mixed_corr_parallel(
  df = data_all,
  categorical_cols = categorical_cols,
  ncores = ncores,
  chunk_size = chunk_size
)

write.csv(
  correlation_matrix_all,
  file.path(corr_dir, "correlation_matrix.csv"),
  row.names = TRUE
)

message("Full mixed correlation matrix saved.")

############################################################
# 7. Prepare continuous variables for heatmaps only
############################################################
continuous_flags <- sapply(
  names(data_all),
  function(v) is_continuous_for_heatmap(data_all[[v]], varname = v, min_unique = 8)
)

continuous_vars <- names(data_all)[continuous_flags]

data_cont <- data_all[, continuous_vars, drop = FALSE] %>%
  mutate(across(everything(), ~ suppressWarnings(as.numeric(as.character(.)))))

valid_continuous <- sapply(
  data_cont,
  function(x) {
    non_missing_n <- sum(!is.na(x))
    if (non_missing_n <= 1) return(FALSE)
    if (all(is.na(x))) return(FALSE)
    if (isTRUE(sd(x, na.rm = TRUE) == 0)) return(FALSE)
    TRUE
  }
)

data_cont <- data_cont[, valid_continuous, drop = FALSE]
all_map_cont <- all_map %>% filter(variable %in% colnames(data_cont))

message("Continuous variables kept for heatmap: ", ncol(data_cont))
message("Variables kept for heatmap: ", paste(colnames(data_cont), collapse = ", "))

############################################################
# 8. Overall continuous heatmap
############################################################
if (ncol(data_cont) >= 2) {
  
  field_order <- c(
    "Lifestyle_Diet",
    "Mental_Health",
    "Environment",
    "Medical_History",
    "Physical_Socio",
    "Biomarker"
  )
  
  biomarker_subgroup_order <- c(
    "Haematology",
    "Renal_Liver_Protein",
    "Lipid_Metabolic_Inflammation"
  )
  
  all_map_cont <- all_map_cont %>%
    mutate(
      field = factor(field, levels = field_order),
      biomarker_subgroup = factor(biomarker_subgroup, levels = biomarker_subgroup_order)
    ) %>%
    arrange(field, biomarker_subgroup, variable)
  
  data_cont <- data_cont[, all_map_cont$variable, drop = FALSE]
  
  cor_mat_cont <- cor(data_cont, use = "pairwise.complete.obs")
  
  annotation <- data.frame(Field = all_map_cont$field)
  rownames(annotation) <- all_map_cont$variable
  
  annotation_colors <- list(
    Field = c(
      Lifestyle_Diet = "#F4A261",
      Mental_Health = "#C77DFF",
      Environment = "#1FBAD6",
      Medical_History = "#D9B44A",
      Physical_Socio = "#43AA8B",
      Biomarker = "#4D96FF"
    )
  )
  
  heat_overall <- pheatmap(
    cor_mat_cont,
    clustering_method = "ward.D2",
    annotation_row = annotation,
    annotation_col = annotation,
    annotation_colors = annotation_colors,
    annotation_legend = TRUE,
    color = colorRampPalette(c("#3B4CC0", "white", "#B40426"))(100),
    breaks = seq(-1, 1, length.out = 101),
    fontsize = 7,
    border_color = NA,
    show_colnames = TRUE,
    show_rownames = TRUE,
    angle_col = 90,
    treeheight_row = 60,
    treeheight_col = 60,
    main = "Overall Correlation Heatmap",
    silent = TRUE
  )
  
  pdf(file.path(corr_dir, "correlation_heatmap_overall.pdf"), width = 12, height = 12)
  grid.newpage()
  grid.draw(heat_overall$gtable)
  dev.off()
  
  png(file.path(corr_dir, "correlation_heatmap_overall.png"), width = 3600, height = 3600, res = 300)
  grid.newpage()
  grid.draw(heat_overall$gtable)
  dev.off()
  
  message("Overall continuous heatmap saved.")
} else {
  message("Skipped overall heatmap: fewer than 2 valid continuous variables.")
}

############################################################
# 9. Biomarker continuous heatmap
############################################################
biomarker_vars_cont <- all_map_cont %>%
  filter(field == "Biomarker") %>%
  pull(variable)

biomarker_vars_cont <- intersect(biomarker_vars_cont, colnames(data_cont))

if (length(biomarker_vars_cont) >= 2) {
  
  biomarker_map_cont <- all_map_cont %>%
    filter(variable %in% biomarker_vars_cont) %>%
    mutate(
      biomarker_subgroup = factor(
        biomarker_subgroup,
        levels = c("Haematology", "Renal_Liver_Protein", "Lipid_Metabolic_Inflammation")
      )
    ) %>%
    arrange(biomarker_subgroup, variable)
  
  biomarker_vars_cont <- biomarker_map_cont$variable
  biomarker_data <- data_cont[, biomarker_vars_cont, drop = FALSE]
  biomarker_cor <- cor(biomarker_data, use = "pairwise.complete.obs")
  
  biomarker_annotation <- data.frame(
    Biomarker_Domain = biomarker_map_cont$biomarker_subgroup
  )
  rownames(biomarker_annotation) <- biomarker_vars_cont
  
  biomarker_annotation_colors <- list(
    Biomarker_Domain = c(
      Haematology = "#4D96FF",
      Renal_Liver_Protein = "#43AA8B",
      Lipid_Metabolic_Inflammation = "#F4A261"
    )
  )
  
  heat_biomarker <- pheatmap(
    biomarker_cor,
    clustering_method = "ward.D2",
    annotation_row = biomarker_annotation,
    annotation_col = biomarker_annotation,
    annotation_colors = biomarker_annotation_colors,
    annotation_legend = TRUE,
    color = colorRampPalette(c("#3B4CC0", "white", "#B40426"))(100),
    breaks = seq(-1, 1, length.out = 101),
    fontsize = 8,
    border_color = NA,
    show_colnames = TRUE,
    show_rownames = TRUE,
    angle_col = 90,
    treeheight_row = 60,
    treeheight_col = 60,
    main = "Biomarker Correlation Heatmap",
    silent = TRUE
  )
  
  pdf(file.path(corr_dir, "correlation_heatmap_biomarker.pdf"), width = 10, height = 10)
  grid.newpage()
  grid.draw(heat_biomarker$gtable)
  dev.off()
  
  png(file.path(corr_dir, "correlation_heatmap_biomarker.png"), width = 3000, height = 3000, res = 300)
  grid.newpage()
  grid.draw(heat_biomarker$gtable)
  dev.off()
  
  message("Biomarker continuous heatmap saved.")
} else {
  message("Skipped biomarker heatmap: fewer than 2 biomarker continuous variables.")
}

message("Analysis completed.")