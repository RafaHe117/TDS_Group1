dir.create("analysis/mediation/inputs", recursive = TRUE, showWarnings = FALSE)

biomarker_candidates <- data.frame(
  biomarker = c(
    "crp",
    "wbc_count",
    "rbc_count",
    "hemoglobin",
    "mcv",
    "platelet_count",
    "mpv",
    "lymphocyte_count",
    "monocyte_count",
    "neutrophil_count",
    "eosinophil_count",
    "basophil_count",
    "creatinine",
    "cystatin_c",
    "urea",
    "albumin",
    "alkaline_phosphatase",
    "alanine_aminotransferase",
    "aspartate_aminotransferase",
    "total_blood_protein",
    "calcium",
    "phosphate",
    "urate",
    "apolipoprotein_a",
    "apolipoprotein_b",
    "cholesterol",
    "hdl_cholesterol",
    "ldl_cholesterol",
    "lipoprotein_a",
    "total_triglyceride",
    "glucose",
    "hba1c",
    "igf_1",
    "blood_vitamin_d"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  biomarker_candidates,
  "analysis/mediation/inputs/biomarker_candidates.csv",
  row.names = FALSE,
  quote = TRUE
)