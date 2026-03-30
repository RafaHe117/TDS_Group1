#!/bin/bash
#$ -N table1_pipeline
#$ -cwd
#$ -o analysis/table1/output/table1_pipeline.out
#$ -e analysis/table1/output/table1_pipeline.err
#$ -l h_rt=02:00:00
#$ -l mem=16G
#$ -pe smp 1

set -euo pipefail

PROJECT_ROOT="/rds/general/project/hda_25-26/live/TDS/TDS_Group1"
TABLE1_DIR="${PROJECT_ROOT}/analysis/table1"
OUT_DIR="${TABLE1_DIR}/output"
ICD_DIR="${TABLE1_DIR}/icd_cat_output"
BYSEX_DIR="${TABLE1_DIR}/output_main_bysex_before_after"
BYSEX_PUB_DIR="${TABLE1_DIR}/output_main_bysex_before_after_publication"
RENDER_DIR="${OUT_DIR}/output_rendered_png"
LOG_DIR="${TABLE1_DIR}/logs"

mkdir -p "${OUT_DIR}" "${LOG_DIR}"

TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="${LOG_DIR}/table1_pipeline_${TIMESTAMP}.log"

exec > >(tee -a "${LOG_FILE}") 2>&1

echo "=================================================="
echo "Table1 pipeline started"
echo "Time: $(date)"
echo "Host: $(hostname)"
echo "Project root: ${PROJECT_ROOT}"
echo "Log file: ${LOG_FILE}"
echo "=================================================="

cd "${PROJECT_ROOT}"

module purge
module load tools/prod
module load R/4.3.2-gfbf-2023a
module load R-bundle-CRAN/2024.11-foss-2024a

export R_LIBS_USER=~/R/x86_64-pc-linux-gnu-library/4.3

run_step () {
  echo
  echo "--------------------------------------------------"
  echo "Running: $1"
  echo "Start: $(date)"
  echo "--------------------------------------------------"
  eval "$2"
  echo "Finished: $1"
  echo "End: $(date)"
}

check_file () {
  if [ ! -f "$1" ]; then
    echo "ERROR: Required file not found:"
    echo "  $1"
    exit 1
  fi
}

check_dir () {
  if [ ! -d "$1" ]; then
    echo "ERROR: Required directory not found:"
    echo "  $1"
    exit 1
  fi
}

echo
echo "Checking required inputs and scripts..."
check_file "${PROJECT_ROOT}/ukb_G1_cleaned.rds"
check_file "${PROJECT_ROOT}/ukb_G1_imputed_final.rds"

check_file "${TABLE1_DIR}/table1_config.R"
check_file "${TABLE1_DIR}/table1_utils.R"
check_file "${TABLE1_DIR}/table1_run.R"
check_file "${TABLE1_DIR}/table1_icd_run.R"
check_file "${TABLE1_DIR}/table1_main_bysex_before_after_run.R"
check_file "${TABLE1_DIR}/table1_main_bysex_before_after_publication.R"
check_file "${TABLE1_DIR}/render_table1_png_base.R"

echo "Pre-checks passed."

run_step "Stage 1 - main table1 core run" \
  "Rscript analysis/table1/table1_run.R"

echo
echo "Checking Stage 1 outputs..."
check_file "${OUT_DIR}/table1_main_outcome_before.csv"
check_file "${OUT_DIR}/table1_main_outcome_after.csv"
check_file "${OUT_DIR}/table1_appendix_outcome_before_missing.csv"
check_file "${OUT_DIR}/table1_appendix_outcome_after.csv"
check_file "${OUT_DIR}/table1_bio_bysex_before.csv"
check_file "${OUT_DIR}/table1_bio_bysex_after.csv"

run_step "Stage 2 - ICD table1 run" \
  "Rscript analysis/table1/table1_icd_run.R"

echo
echo "Checking Stage 2 outputs..."
check_dir "${ICD_DIR}"
check_file "${ICD_DIR}/table1_icd_outcome_before_missing.csv"
check_file "${ICD_DIR}/table1_icd_outcome_after.csv"
check_file "${ICD_DIR}/table1_icd_outcome_before_missing.png"
check_file "${ICD_DIR}/table1_icd_outcome_after.png"

run_step "Stage 3 - by-sex main table1 run" \
  "Rscript analysis/table1/table1_main_bysex_before_after_run.R"

echo
echo "Checking Stage 3 outputs..."
check_dir "${BYSEX_DIR}"
check_file "${BYSEX_DIR}/table1_main_female_outcome_before.csv"
check_file "${BYSEX_DIR}/table1_main_female_outcome_after.csv"
check_file "${BYSEX_DIR}/table1_main_male_outcome_before.csv"
check_file "${BYSEX_DIR}/table1_main_male_outcome_after.csv"
check_file "${BYSEX_DIR}/table1_main_female_outcome_before.png"
check_file "${BYSEX_DIR}/table1_main_female_outcome_after.png"
check_file "${BYSEX_DIR}/table1_main_male_outcome_before.png"
check_file "${BYSEX_DIR}/table1_main_male_outcome_after.png"

run_step "Stage 4 - by-sex publication table1 run" \
  "Rscript analysis/table1/table1_main_bysex_before_after_publication.R"

echo
echo "Checking Stage 4 outputs..."
check_dir "${BYSEX_PUB_DIR}"
check_file "${BYSEX_PUB_DIR}/table1_main_female_outcome_before_publication.csv"
check_file "${BYSEX_PUB_DIR}/table1_main_female_outcome_after_publication.csv"
check_file "${BYSEX_PUB_DIR}/table1_main_male_outcome_before_publication.csv"
check_file "${BYSEX_PUB_DIR}/table1_main_male_outcome_after_publication.csv"
check_file "${BYSEX_PUB_DIR}/table1_main_female_outcome_before_publication.png"
check_file "${BYSEX_PUB_DIR}/table1_main_female_outcome_after_publication.png"
check_file "${BYSEX_PUB_DIR}/table1_main_male_outcome_before_publication.png"
check_file "${BYSEX_PUB_DIR}/table1_main_male_outcome_after_publication.png"

run_step "Stage 5 - render core output CSVs to PNG" \
  "Rscript analysis/table1/render_table1_png_base.R"

echo
echo "Checking Stage 5 outputs..."
check_dir "${RENDER_DIR}"

check_file "${RENDER_DIR}/table1_main_outcome_before_page01.png"
check_file "${RENDER_DIR}/table1_main_outcome_after_page01.png"
check_file "${RENDER_DIR}/table1_appendix_outcome_before_missing_page01.png"
check_file "${RENDER_DIR}/table1_appendix_outcome_after_page01.png"
check_file "${RENDER_DIR}/table1_bio_bysex_before_page01.png"
check_file "${RENDER_DIR}/table1_bio_bysex_after_page01.png"

check_file "${RENDER_DIR}/table1_main_outcome_before_LONG.png"
check_file "${RENDER_DIR}/table1_main_outcome_after_LONG.png"
check_file "${RENDER_DIR}/table1_appendix_outcome_before_missing_LONG.png"
check_file "${RENDER_DIR}/table1_appendix_outcome_after_LONG.png"
check_file "${RENDER_DIR}/table1_bio_bysex_before_LONG.png"
check_file "${RENDER_DIR}/table1_bio_bysex_after_LONG.png"

echo
echo "=================================================="
echo "Table1 pipeline completed successfully"
echo "Time: $(date)"
echo "Core output directory: ${OUT_DIR}"
echo "ICD output directory: ${ICD_DIR}"
echo "By-sex output directory: ${BYSEX_DIR}"
echo "Publication by-sex output directory: ${BYSEX_PUB_DIR}"
echo "Rendered PNG directory: ${RENDER_DIR}"
echo "Log file: ${LOG_FILE}"
echo "=================================================="
