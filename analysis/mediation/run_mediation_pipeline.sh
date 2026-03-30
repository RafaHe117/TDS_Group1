#!/bin/bash
#$ -N med_pipeline
#$ -cwd
#$ -o analysis/mediation/mediation_pipeline.out
#$ -e analysis/mediation/mediation_pipeline.err
#$ -l h_rt=24:00:00
#$ -l mem=32G
#$ -pe smp 4

set -euo pipefail

PROJECT_ROOT="/rds/general/project/hda_25-26/live/TDS/TDS_Group1"
MED_DIR="${PROJECT_ROOT}/analysis/mediation"
INPUT_DIR="${MED_DIR}/inputs"
OUTPUT_DIR="${MED_DIR}/outputs"
LOG_DIR="${MED_DIR}/logs"

mkdir -p "${LOG_DIR}"

TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="${LOG_DIR}/mediation_pipeline_${TIMESTAMP}.log"

exec > >(tee -a "${LOG_FILE}") 2>&1

echo "=================================================="
echo "Mediation pipeline started"
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

echo
echo "Checking required directories..."
[ -d "${MED_DIR}" ] || { echo "ERROR: Missing mediation directory: ${MED_DIR}"; exit 1; }
[ -d "${INPUT_DIR}" ] || { echo "ERROR: Missing inputs directory: ${INPUT_DIR}"; exit 1; }

echo "Checking required static input files..."
check_file "${INPUT_DIR}/formal_mediation_config.csv"

echo "Pre-checks passed."

run_step "Stage 1a - create biomarker candidates" \
  "Rscript analysis/mediation/create_biomarker_candidates.R"

run_step "Stage 1b - build exposure list" \
  "Rscript analysis/mediation/build_exposure_list.R"

echo
echo "Checking generated Stage 1 files..."
check_file "${INPUT_DIR}/biomarker_candidates.csv"
check_file "${INPUT_DIR}/raw_exposure_universe.csv"
check_file "${INPUT_DIR}/selected_terms_main.csv"
check_file "${INPUT_DIR}/selected_terms_female.csv"
check_file "${INPUT_DIR}/selected_terms_male.csv"

run_step "Stage 2a - pooled biomarker screening" \
  "Rscript analysis/mediation/run_biomarker_models_main.R"

run_step "Stage 2b - female biomarker screening" \
  "Rscript analysis/mediation/run_biomarker_models_sex.R female 0"

run_step "Stage 2c - male biomarker screening" \
  "Rscript analysis/mediation/run_biomarker_models_sex.R male 1"

echo
echo "Checking Stage 2 outputs..."
check_file "${OUTPUT_DIR}/main/stable_links_main.csv"
check_file "${OUTPUT_DIR}/female/stable_links_female.csv"
check_file "${OUTPUT_DIR}/male/stable_links_male.csv"

run_step "Stage 3 - final outcome models" \
  "Rscript analysis/mediation/run_final_outcome_models.R"

echo
echo "Checking Stage 3 outputs..."
check_file "${OUTPUT_DIR}/final_outcome_models/main/final_model_coefficients_main.csv"
check_file "${OUTPUT_DIR}/final_outcome_models/female/final_model_coefficients_female.csv"
check_file "${OUTPUT_DIR}/final_outcome_models/male/final_model_coefficients_male.csv"

run_step "Stage 4 - build final shortlist" \
  "Rscript analysis/mediation/build_final_shortlist.R"

echo
echo "Checking Stage 4 outputs..."
check_file "${OUTPUT_DIR}/final_shortlists/strict_shortlist_main.csv"
check_file "${OUTPUT_DIR}/final_shortlists/strict_shortlist_female.csv"
check_file "${OUTPUT_DIR}/final_shortlists/strict_shortlist_male.csv"
check_file "${OUTPUT_DIR}/final_shortlists/summary_shortlist_main.csv"
check_file "${OUTPUT_DIR}/final_shortlists/summary_shortlist_female.csv"
check_file "${OUTPUT_DIR}/final_shortlists/summary_shortlist_male.csv"

echo
echo "Checking formal mediation config and pair files..."
check_file "${INPUT_DIR}/formal_mediation_config.csv"
check_file "${INPUT_DIR}/pairs_main_chunk1.csv"
check_file "${INPUT_DIR}/pairs_main_chunk2.csv"
check_file "${INPUT_DIR}/pairs_female_chunk1.csv"
check_file "${INPUT_DIR}/pairs_female_chunk2.csv"

run_step "Stage 5a - formal mediation refit main chunk 1" \
  "Rscript analysis/mediation/run_formal_mediation_refit_chunk.R analysis/mediation/inputs/pairs_main_chunk1.csv main_chunk1 all"

run_step "Stage 5b - formal mediation refit main chunk 2" \
  "Rscript analysis/mediation/run_formal_mediation_refit_chunk.R analysis/mediation/inputs/pairs_main_chunk2.csv main_chunk2 all"

run_step "Stage 5c - formal mediation refit female chunk 1" \
  "Rscript analysis/mediation/run_formal_mediation_refit_chunk.R analysis/mediation/inputs/pairs_female_chunk1.csv female_chunk1 female"

run_step "Stage 5d - formal mediation refit female chunk 2" \
  "Rscript analysis/mediation/run_formal_mediation_refit_chunk.R analysis/mediation/inputs/pairs_female_chunk2.csv female_chunk2 female"

echo
echo "Checking Stage 5 outputs..."
check_file "${OUTPUT_DIR}/formal_mediation_refit_split/main_chunk1/formal_mediation_results_main_chunk1.csv"
check_file "${OUTPUT_DIR}/formal_mediation_refit_split/main_chunk2/formal_mediation_results_main_chunk2.csv"
check_file "${OUTPUT_DIR}/formal_mediation_refit_split/female_chunk1/formal_mediation_results_female_chunk1.csv"
check_file "${OUTPUT_DIR}/formal_mediation_refit_split/female_chunk2/formal_mediation_results_female_chunk2.csv"

echo
echo "Male formal mediation refit is not run in this pipeline because no strict shortlisted male pathway entered the formal refit stage."

run_step "Stage 6 - make mediation outputs" \
  "Rscript analysis/mediation/make_mediation_outputs.R"

echo
echo "Checking Stage 6 outputs..."
check_file "${OUTPUT_DIR}/final_mediation_outputs/final_mediation_results_main.csv"
check_file "${OUTPUT_DIR}/final_mediation_outputs/final_mediation_results_female.csv"
check_file "${OUTPUT_DIR}/final_mediation_outputs/display_table_main.csv"
check_file "${OUTPUT_DIR}/final_mediation_outputs/display_table_female.csv"
check_file "${OUTPUT_DIR}/final_mediation_outputs/pathway_counts_summary.csv"

echo
echo "=================================================="
echo "Mediation pipeline completed successfully"
echo "Time: $(date)"
echo "Main output directory: ${OUTPUT_DIR}/final_mediation_outputs"
echo "Log file: ${LOG_FILE}"
echo "=================================================="
