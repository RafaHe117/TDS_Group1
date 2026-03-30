#!/bin/bash
#$ -N famd_pipeline
#$ -cwd
#$ -o analysis/FAMD/famd_pipeline.out
#$ -e analysis/FAMD/famd_pipeline.err
#$ -l h_rt=08:00:00
#$ -l mem=24G
#$ -pe smp 2

set -euo pipefail

PROJECT_ROOT="/rds/general/project/hda_25-26/live/TDS/TDS_Group1"
FAMD_DIR="${PROJECT_ROOT}/analysis/FAMD"
OUT_DIR="${FAMD_DIR}/FAMD_output"
LOG_DIR="${FAMD_DIR}/logs"

mkdir -p "${LOG_DIR}"

TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="${LOG_DIR}/famd_pipeline_${TIMESTAMP}.log"

exec > >(tee -a "${LOG_FILE}") 2>&1

echo "=================================================="
echo "FAMD pipeline started"
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
echo "Checking required inputs..."
check_file "${PROJECT_ROOT}/ukb_G1_imputed_final.rds"
check_file "${PROJECT_ROOT}/split_imputed_data/ukb_G1_train_imputed.rds"
check_file "${FAMD_DIR}/check_famd_variance_only.R"
check_file "${FAMD_DIR}/pseudo_famd_run.R"

echo "Pre-checks passed."

run_step "Stage 1 - variance check" \
  "Rscript analysis/FAMD/check_famd_variance_only.R"

run_step "Stage 2 - pseudo FAMD main run" \
  "Rscript analysis/FAMD/pseudo_famd_run.R"

echo
echo "Checking expected output structure..."
check_dir "${OUT_DIR}"
check_dir "${OUT_DIR}/final"
check_dir "${OUT_DIR}/train"

# -----------------------------
# FINAL dataset outputs
# -----------------------------
check_file "${OUT_DIR}/final/pseudo_famd1a_baseline_main_no_alcohol_scree.png"
check_file "${OUT_DIR}/final/pseudo_famd1a_baseline_main_no_alcohol_indiv_by_outcome.png"
check_file "${OUT_DIR}/final/pseudo_famd1a_baseline_main_no_alcohol_indiv_by_sex.png"
check_file "${OUT_DIR}/final/pseudo_famd1a_baseline_main_no_alcohol_contrib_pc1.png"
check_file "${OUT_DIR}/final/pseudo_famd1a_baseline_main_no_alcohol_contrib_pc2.png"
check_file "${OUT_DIR}/final/pseudo_famd1a_baseline_main_no_alcohol_variance_explained.csv"
check_file "${OUT_DIR}/final/pseudo_famd1a_baseline_main_no_alcohol_scores.csv"
check_file "${OUT_DIR}/final/pseudo_famd1a_baseline_main_no_alcohol_contrib_pc1.csv"
check_file "${OUT_DIR}/final/pseudo_famd1a_baseline_main_no_alcohol_contrib_pc2.csv"

check_file "${OUT_DIR}/final/pseudo_famd1b_baseline_plus_alcohol_sensitivity_scree.png"
check_file "${OUT_DIR}/final/pseudo_famd1b_baseline_plus_alcohol_sensitivity_indiv_by_outcome.png"
check_file "${OUT_DIR}/final/pseudo_famd1b_baseline_plus_alcohol_sensitivity_indiv_by_sex.png"
check_file "${OUT_DIR}/final/pseudo_famd1b_baseline_plus_alcohol_sensitivity_contrib_pc1.png"
check_file "${OUT_DIR}/final/pseudo_famd1b_baseline_plus_alcohol_sensitivity_contrib_pc2.png"
check_file "${OUT_DIR}/final/pseudo_famd1b_baseline_plus_alcohol_sensitivity_variance_explained.csv"
check_file "${OUT_DIR}/final/pseudo_famd1b_baseline_plus_alcohol_sensitivity_scores.csv"
check_file "${OUT_DIR}/final/pseudo_famd1b_baseline_plus_alcohol_sensitivity_contrib_pc1.csv"
check_file "${OUT_DIR}/final/pseudo_famd1b_baseline_plus_alcohol_sensitivity_contrib_pc2.csv"

check_file "${OUT_DIR}/final/pseudo_famd2_biological_supplementary_scree.png"
check_file "${OUT_DIR}/final/pseudo_famd2_biological_supplementary_indiv_by_outcome.png"
check_file "${OUT_DIR}/final/pseudo_famd2_biological_supplementary_indiv_by_sex.png"
check_file "${OUT_DIR}/final/pseudo_famd2_biological_supplementary_contrib_pc1.png"
check_file "${OUT_DIR}/final/pseudo_famd2_biological_supplementary_contrib_pc2.png"
check_file "${OUT_DIR}/final/pseudo_famd2_biological_supplementary_variance_explained.csv"
check_file "${OUT_DIR}/final/pseudo_famd2_biological_supplementary_scores.csv"
check_file "${OUT_DIR}/final/pseudo_famd2_biological_supplementary_contrib_pc1.csv"
check_file "${OUT_DIR}/final/pseudo_famd2_biological_supplementary_contrib_pc2.csv"

# -----------------------------
# TRAIN dataset outputs
# -----------------------------
check_file "${OUT_DIR}/train/pseudo_famd1a_baseline_main_no_alcohol_scree.png"
check_file "${OUT_DIR}/train/pseudo_famd1a_baseline_main_no_alcohol_indiv_by_outcome.png"
check_file "${OUT_DIR}/train/pseudo_famd1a_baseline_main_no_alcohol_indiv_by_sex.png"
check_file "${OUT_DIR}/train/pseudo_famd1a_baseline_main_no_alcohol_contrib_pc1.png"
check_file "${OUT_DIR}/train/pseudo_famd1a_baseline_main_no_alcohol_contrib_pc2.png"
check_file "${OUT_DIR}/train/pseudo_famd1a_baseline_main_no_alcohol_variance_explained.csv"
check_file "${OUT_DIR}/train/pseudo_famd1a_baseline_main_no_alcohol_scores.csv"
check_file "${OUT_DIR}/train/pseudo_famd1a_baseline_main_no_alcohol_contrib_pc1.csv"
check_file "${OUT_DIR}/train/pseudo_famd1a_baseline_main_no_alcohol_contrib_pc2.csv"

check_file "${OUT_DIR}/train/pseudo_famd1b_baseline_plus_alcohol_sensitivity_scree.png"
check_file "${OUT_DIR}/train/pseudo_famd1b_baseline_plus_alcohol_sensitivity_indiv_by_outcome.png"
check_file "${OUT_DIR}/train/pseudo_famd1b_baseline_plus_alcohol_sensitivity_indiv_by_sex.png"
check_file "${OUT_DIR}/train/pseudo_famd1b_baseline_plus_alcohol_sensitivity_contrib_pc1.png"
check_file "${OUT_DIR}/train/pseudo_famd1b_baseline_plus_alcohol_sensitivity_contrib_pc2.png"
check_file "${OUT_DIR}/train/pseudo_famd1b_baseline_plus_alcohol_sensitivity_variance_explained.csv"
check_file "${OUT_DIR}/train/pseudo_famd1b_baseline_plus_alcohol_sensitivity_scores.csv"
check_file "${OUT_DIR}/train/pseudo_famd1b_baseline_plus_alcohol_sensitivity_contrib_pc1.csv"
check_file "${OUT_DIR}/train/pseudo_famd1b_baseline_plus_alcohol_sensitivity_contrib_pc2.csv"

check_file "${OUT_DIR}/train/pseudo_famd2_biological_supplementary_scree.png"
check_file "${OUT_DIR}/train/pseudo_famd2_biological_supplementary_indiv_by_outcome.png"
check_file "${OUT_DIR}/train/pseudo_famd2_biological_supplementary_indiv_by_sex.png"
check_file "${OUT_DIR}/train/pseudo_famd2_biological_supplementary_contrib_pc1.png"
check_file "${OUT_DIR}/train/pseudo_famd2_biological_supplementary_contrib_pc2.png"
check_file "${OUT_DIR}/train/pseudo_famd2_biological_supplementary_variance_explained.csv"
check_file "${OUT_DIR}/train/pseudo_famd2_biological_supplementary_scores.csv"
check_file "${OUT_DIR}/train/pseudo_famd2_biological_supplementary_contrib_pc1.csv"
check_file "${OUT_DIR}/train/pseudo_famd2_biological_supplementary_contrib_pc2.csv"

echo
echo "=================================================="
echo "FAMD pipeline completed successfully"
echo "Time: $(date)"
echo "Main output directory: ${OUT_DIR}"
echo "Log file: ${LOG_FILE}"
echo "=================================================="
