# ============================================================
# render_table1_png.R
# Final publication-style renderer for 6 Table 1 CSV outputs
# ============================================================

# ---------------- packages ----------------
# Only base/grid used to avoid package installation issues
library(grid)

# ---------------- paths ----------------
indir  <- file.path("analysis", "output")
outdir <- file.path("analysis", "output_rendered_png")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---------------- input files ----------------
files <- c(
  "table1_main_outcome_before.csv",
  "table1_main_outcome_after.csv",
  "table1_appendix_outcome_before_missing.csv",
  "table1_appendix_outcome_after.csv",
  "table1_bio_bysex_before.csv",
  "table1_bio_bysex_after.csv"
)

# ============================================================
# 1. LABEL DICTIONARY
# ============================================================

pretty_var <- c(
  age = "Age, years",
  sex = "Sex",
  ethnicity_5cat = "Ethnicity",
  imd_quintile = "Townsend deprivation quintile",
  urban_rural_2cat = "Urban/rural residence",
  education_2 = "Education",
  employment_2cat = "Employment status",
  living_with_partner = "Living with partner",
  work_hours_week_clean = "Working hours per week",
  bmi = "Body mass index",
  sbp_mean = "Systolic blood pressure",
  dbp_mean = "Diastolic blood pressure",
  diabetes_bin = "Diabetes",
  famhx_cvd = "Family history of CVD",
  smoking_3cat = "Smoking status",
  alcohol_combined = "Alcohol consumption",
  birth_weight_clean = "Birth weight",
  total_met_group = "Physical activity",
  tv_group = "Television viewing",
  sleep_quality_score = "Sleep quality score",
  fruit_intake_fresh = "Fresh fruit intake",
  oily_fish_3cat = "Oily fish intake",
  salt_3cat = "Salt added to food",
  coffee_intake = "Coffee intake",
  tea_intake = "Tea intake",
  redwine_group = "Red wine intake",
  neuroticism_score = "Neuroticism score",
  sf_score = "Social support score",
  mh_satis_mean_score = "Mental health satisfaction score",
  stress_group_2yr = "Stress in past 2 years",
  mood_disorder = "Mood disorder",
  medication_category = "Medication category",
  pm2_5_2010 = "PM2.5",
  pm10_2010 = "PM10",
  no2_2010 = "NO2",
  noise_24h = "24-hour noise",
  greenspace_pct_1000m = "Greenspace within 1000m",
  temp_average = "Average temperature",
  crp = "C-reactive protein",
  hba1c = "HbA1c",
  glucose = "Glucose",
  hdl_cholesterol = "HDL cholesterol",
  ldl_cholesterol = "LDL cholesterol",
  cholesterol = "Total cholesterol",
  total_triglyceride = "Triglycerides",
  creatinine = "Creatinine",
  cystatin_c = "Cystatin C",
  urea = "Urea",
  urate = "Urate",
  albumin = "Albumin",
  total_bilirubin = "Total bilirubin",
  direct_bilirubin = "Direct bilirubin",
  alkaline_phosphatase = "Alkaline phosphatase",
  alanine_aminotransferase = "Alanine aminotransferase",
  aspartate_aminotransferase = "Aspartate aminotransferase",
  calcium = "Calcium",
  phosphate = "Phosphate",
  total_blood_protein = "Total protein",
  blood_vitamin_d = "Vitamin D",
  apolipoprotein_a = "Apolipoprotein A",
  apolipoprotein_b = "Apolipoprotein B",
  lipoprotein_a = "Lipoprotein(a)",
  wbc_count = "White blood cell count",
  rbc_count = "Red blood cell count",
  hemoglobin = "Haemoglobin",
  platelet_count = "Platelet count",
  mpv = "Mean platelet volume",
  mcv = "Mean corpuscular volume",
  lymphocyte_count = "Lymphocyte count",
  monocyte_count = "Monocyte count",
  neutrophil_count = "Neutrophil count",
  eosinophil_count = "Eosinophil count",
  basophil_count = "Basophil count",
  ICD_Diabetes = "ICD diabetes",
  ICD_Lipidemia = "ICD dyslipidaemia",
  ICD_Hypertension = "ICD hypertension",
  ICD_AF = "ICD atrial fibrillation",
  ICD_CKD = "ICD chronic kidney disease",
  ICD_Mental_Health = "ICD mental health condition",
  ICD_Migraine = "ICD migraine",
  ICD_Atopy = "ICD atopy",
  ICD_Autoimmune = "ICD autoimmune disease"
)

# ============================================================
# 2. GENERAL HELPERS
# ============================================================

wrap_text <- function(x, width = 28) {
  x <- ifelse(is.na(x), "", as.character(x))
  vapply(
    x,
    function(s) {
      if (!nzchar(s)) return("")
      paste(strwrap(s, width = width), collapse = "\n")
    },
    character(1)
  )
}

fmt_header <- function(x) {
  x <- as.character(x)
  x[x == "variable"] <- "Variable"
  x[x == "level"] <- ""
  x[x == "Overall"] <- "Overall"
  x[x == "p_value"] <- "P-value"
  x[x == "missing_overall"] <- "Missing, n (%)"
  x[x == "No incident CVD"] <- "No incident CVD"
  x[x == "Incident CVD"] <- "Incident CVD"
  x[x == "Female"] <- "Female"
  x[x == "Male"] <- "Male"
  x
}

fmt_pval <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA")] <- ""
  
  out <- x
  is_num <- suppressWarnings(!is.na(as.numeric(x)))
  
  num <- suppressWarnings(as.numeric(x[is_num]))
  out[is_num] <- ifelse(num < 0.001, "<0.001", sprintf("%.3f", num))
  
  out
}

clean_level <- function(x) {
  x <- ifelse(is.na(x), "", as.character(x))
  x <- gsub("^\\s+", "", x)
  x <- gsub("_", " ", x)
  x
}

clean_variable <- function(x) {
  x <- ifelse(is.na(x), "", as.character(x))
  x2 <- x
  idx <- nzchar(x2)
  x2[idx] <- ifelse(x2[idx] %in% names(pretty_var), pretty_var[x2[idx]], x2[idx])
  x2
}

# detect rows with huge level labels and shrink them a bit more
is_long_level_row <- function(x) {
  nchar(gsub("\n", " ", x, fixed = TRUE)) > 45
}

row_line_count <- function(df) {
  cols <- lapply(df, function(x) {
    x <- ifelse(is.na(x), "", as.character(x))
    sapply(strsplit(x, "\n", fixed = TRUE), length)
  })
  Reduce(pmax, cols)
}

# ============================================================
# 3. TABLE PREPARATION
# ============================================================

prepare_table_df <- function(df, file_name) {
  names(df) <- fmt_header(names(df))
  
  if ("Variable" %in% names(df)) {
    df$Variable <- clean_variable(df$Variable)
    df$Variable <- wrap_text(df$Variable, width = 24)
  }
  
  if ("" %in% names(df)) {
    df[[which(names(df) == "")[1]]] <- clean_level(df[[which(names(df) == "")[1]]])
    df[[which(names(df) == "")[1]]] <- wrap_text(df[[which(names(df) == "")[1]]], width = 34)
  }
  
  if ("P-value" %in% names(df)) {
    df[["P-value"]] <- fmt_pval(df[["P-value"]])
  }
  
  # Wrap other headers if needed
  names(df) <- wrap_text(names(df), width = 18)
  
  # For bio tables, give a bit more wrapping in level column
  if (grepl("bio", file_name, ignore.case = TRUE) && "" %in% names(df)) {
    df[[which(names(df) == "")[1]]] <- wrap_text(clean_level(df[[which(names(df) == "")[1]]]), width = 28)
  }
  
  # For appendix, be more aggressive in wrapping long levels
  if (grepl("appendix", file_name, ignore.case = TRUE) && "" %in% names(df)) {
    df[[which(names(df) == "")[1]]] <- wrap_text(clean_level(df[[which(names(df) == "")[1]]]), width = 26)
  }
  
  df
}

# ============================================================
# 4. PAGINATION
# ============================================================

paginate_df <- function(df, max_lines_per_page = 38) {
  pages <- list()
  current_idx <- integer(0)
  current_lines <- 2  # header budget
  
  rl <- row_line_count(df)
  
  for (i in seq_len(nrow(df))) {
    need <- max(1, rl[i])
    
    if (length(current_idx) > 0 && current_lines + need > max_lines_per_page) {
      pages[[length(pages) + 1]] <- df[current_idx, , drop = FALSE]
      current_idx <- i
      current_lines <- 2 + need
    } else {
      current_idx <- c(current_idx, i)
      current_lines <- current_lines + need
    }
  }
  
  if (length(current_idx) > 0) {
    pages[[length(pages) + 1]] <- df[current_idx, , drop = FALSE]
  }
  
  pages
}

# ============================================================
# 5. COLUMN LAYOUT RULES
# ============================================================

get_col_widths <- function(df, file_name) {
  nms <- names(df)
  
  # with missing column
  if ("Missing, n (%)" %in% nms) {
    # Variable | level | overall | g1 | g2 | missing | p
    return(c(0.16, 0.25, 0.13, 0.13, 0.13, 0.14, 0.06))
  }
  
  # 6-column standard
  if (grepl("appendix", file_name, ignore.case = TRUE)) {
    # more space to level column for appendix
    return(c(0.16, 0.31, 0.15, 0.15, 0.15, 0.08))
  }
  
  if (grepl("bio", file_name, ignore.case = TRUE)) {
    return(c(0.18, 0.30, 0.15, 0.15, 0.15, 0.07))
  }
  
  # main table default
  c(0.18, 0.24, 0.17, 0.17, 0.17, 0.07)
}

# ============================================================
# 6. DRAW SINGLE PAGE
# ============================================================

render_one_page <- function(df_page, outfile, file_name,
                            page_width = 2600,
                            page_height = 3200,
                            font_family = "Times") {
  
  png(outfile, width = page_width, height = page_height, res = 220, bg = "white")
  grid.newpage()
  
  nr <- nrow(df_page)
  nc <- ncol(df_page)
  
  col_widths <- get_col_widths(df_page, file_name)
  if (length(col_widths) != nc) {
    col_widths <- rep(1 / nc, nc)
  }
  col_widths <- col_widths / sum(col_widths)
  
  x_left <- cumsum(c(0, head(col_widths, -1)))
  
  # dynamic row heights
  body_lines <- row_line_count(df_page)
  row_weights <- c(1.4, pmax(1, body_lines * 0.95))
  row_heights <- row_weights / sum(row_weights)
  
  y_bottoms <- 1 - cumsum(row_heights)
  
  headers <- names(df_page)
  
  # draw grid and text
  for (r in 0:nr) {
    for (c in seq_len(nc)) {
      x0 <- x_left[c]
      w  <- col_widths[c]
      y0 <- y_bottoms[r + 1]
      h  <- row_heights[r + 1]
      
      fill_col <- if (r == 0) "#F3F3F3" else "white"
      
      grid.rect(
        x = x0, y = y0,
        width = w, height = h,
        just = c("left", "bottom"),
        gp = gpar(col = "grey75", fill = fill_col, lwd = 0.8)
      )
      
      txt <- if (r == 0) headers[c] else as.character(df_page[r, c])
      
      # first two columns left-aligned, others right-aligned
      left_align <- c <= 2
      x_txt <- if (left_align) x0 + 0.006 else x0 + w - 0.006
      
      # slightly smaller font for appendix/bio long rows
      fontsize_body <- 10
      if (grepl("appendix", file_name, ignore.case = TRUE)) fontsize_body <- 9.2
      if (grepl("bio", file_name, ignore.case = TRUE)) fontsize_body <- 9.0
      
      # if level cell is very long, make it smaller
      if (r > 0 && c == 2 && is_long_level_row(txt)) {
        fontsize_body <- fontsize_body - 0.8
      }
      
      grid.text(
        label = txt,
        x = x_txt,
        y = y0 + h / 2,
        just = c(if (left_align) "left" else "right", "center"),
        gp = gpar(
          fontsize   = if (r == 0) 11 else fontsize_body,
          fontface   = if (r == 0) "bold" else "plain",
          fontfamily = font_family,
          lineheight = 1.05,
          col = "black"
        )
      )
    }
  }
  
  dev.off()
}

# ============================================================
# 7. RENDER CSV -> PNG(S)
# ============================================================

render_csv_paginated <- function(file, indir, outdir, font_family = "Times") {
  inpath <- file.path(indir, file)
  if (!file.exists(inpath)) {
    message("Skip (not found): ", inpath)
    return(invisible(NULL))
  }
  
  df <- read.csv(inpath, stringsAsFactors = FALSE, check.names = FALSE)
  df <- prepare_table_df(df, file)
  
  is_appendix <- grepl("appendix", file, ignore.case = TRUE)
  is_bio      <- grepl("bio", file, ignore.case = TRUE)
  is_main     <- grepl("main", file, ignore.case = TRUE)
  
  max_lines <- if (is_appendix) 24 else if (is_bio) 28 else 40
  pages <- paginate_df(df, max_lines_per_page = max_lines)
  
  base <- sub("\\.csv$", "", file, ignore.case = TRUE)
  
  # remove old same-prefix pngs
  old_pngs <- list.files(outdir, pattern = paste0("^", gsub("\\.", "\\\\.", base), "(_page[0-9]+)?\\.png$"), full.names = TRUE)
  if (length(old_pngs) > 0) unlink(old_pngs)
  
  page_width  <- if (is_appendix) 3000 else if (is_bio) 2800 else 2400
  page_height <- if (is_appendix) 3800 else if (is_bio) 3200 else 2600
  
  if (length(pages) == 1) {
    outfile <- file.path(outdir, paste0(base, ".png"))
    message("Rendering: ", file, " -> ", outfile)
    render_one_page(
      df_page = pages[[1]],
      outfile = outfile,
      file_name = file,
      page_width = page_width,
      page_height = page_height,
      font_family = font_family
    )
  } else {
    for (i in seq_along(pages)) {
      outfile <- file.path(outdir, paste0(base, "_page", i, ".png"))
      message("Rendering: ", file, " -> ", outfile)
      render_one_page(
        df_page = pages[[i]],
        outfile = outfile,
        file_name = file,
        page_width = page_width,
        page_height = page_height,
        font_family = font_family
      )
    }
  }
}

# ============================================================
# 8. RUN ALL
# ============================================================

for (f in files) {
  render_csv_paginated(
    file = f,
    indir = indir,
    outdir = outdir,
    font_family = "Times"
  )
}

message("Done. Publication-style PNG tables saved to: ", outdir)