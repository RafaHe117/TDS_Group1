setwd("/rds/general/project/hda_25-26/live/TDS/TDS_Group1")

input_file  <- "/rds/general/project/hda_25-26/live/TDS/TDS_Group1/extraction_and_recoding/outputs/recoded.rds"
output_file <- "/rds/general/project/hda_25-26/live/TDS/TDS_Group1/ukb_G1_raw.rds"

############################################################
# CHUNK 1: Load recoded data, keep baseline (.0.0) fields and
# retain all columns for multi-response variables
############################################################

data <- readRDS(input_file)

# use rownames as the true UKB eid
data$eid <- as.integer(rownames(data))

multi_vars <- c(
  "education",
  "household_relationship",
  "employment_status",
  "iibs_2yr",
  "chol_bp_diabetes_hormone_medication",
  "chol_bp_diabetes_medication",
  "father_illness",
  "mother_illness",
  "sibling_illness",
  "primary_diagnoses_icd10",
  "primary_diagnoses_icd9",
  "secondary_diagnoses_icd10",
  "secondary_diagnoses_icd9"
)

cn <- colnames(data)

# keep all array columns for multi-response variables
multi_pattern <- paste0("^(", paste(multi_vars, collapse = "|"), ")\\.0\\.\\d+$")
multi_cols <- grep(multi_pattern, cn, value = TRUE)

# keep only baseline instance (.0.0) for other variables
baseline_cols <- grep("\\.0\\.0$", setdiff(cn, multi_cols), value = TRUE)

# subset dataset
cols_to_keep <- unique(c("eid", multi_cols, baseline_cols))
data_out <- data[, cols_to_keep, drop = FALSE]

# remove ".0.0" suffix for non-multi variables
new_names <- colnames(data_out)
is_multi <- sapply(new_names, function(x) any(startsWith(x, paste0(multi_vars, "."))))
new_names[!is_multi] <- sub("\\.0\\.0$", "", new_names[!is_multi])

if (anyDuplicated(new_names)) stop("Duplicate column names after renaming")

colnames(data_out) <- new_names


############################################################
# CHUNK 2: Merge CVD outcomes and derive event indicators
############################################################

cvd <- readRDS("cvd_events.rds")

data_out$eid <- as.integer(data_out$eid)
cvd$eid <- as.integer(cvd$eid)

cvd$date <- as.Date(cvd$date)
data_out$date_recr <- as.Date(data_out$date_recr)

# use earliest CVD date per eid
cvd_first <- aggregate(date ~ eid, data = cvd, FUN = min)
names(cvd_first)[2] <- "cvd_first_date"

data_out$cvd_first_date <- cvd_first$cvd_first_date[match(data_out$eid, cvd_first$eid)]

# derive event indicators as logical TRUE/FALSE
data_out$cvd_event     <- !is.na(data_out$cvd_first_date)
data_out$cvd_prevalent <- data_out$cvd_event & data_out$cvd_first_date <  data_out$date_recr
data_out$cvd_incident  <- data_out$cvd_event & data_out$cvd_first_date >= data_out$date_recr

stopifnot(all(data_out$cvd_prevalent <= data_out$cvd_event))
stopifnot(all(data_out$cvd_incident  <= data_out$cvd_event))
stopifnot(all(data_out$cvd_prevalent + data_out$cvd_incident == data_out$cvd_event))

# move outcome variables to the front
data_out <- data_out[, c(
  "eid",
  "cvd_event",
  "cvd_prevalent",
  "cvd_incident",
  "cvd_first_date",
  setdiff(names(data_out), c(
    "eid",
    "cvd_event",
    "cvd_prevalent",
    "cvd_incident",
    "cvd_first_date"
  ))
), drop = FALSE]

saveRDS(data_out, output_file)

cat("Saved:", output_file, "\n")
cat("Final:", nrow(data_out), "rows x", ncol(data_out), "cols\n")