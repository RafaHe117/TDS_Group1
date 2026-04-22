#!/bin/bash
#PBS -N med_male_refit
#PBS -o /rds/general/project/hda_25-26/live/TDS/anw16/TDS_Group1/modelling_script/08_mediation_outdated/male_formal_mediation.out
#PBS -e /rds/general/project/hda_25-26/live/TDS/anw16/TDS_Group1/modelling_script/08_mediation_outdated/male_formal_mediation.err
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=2:mem=16gb

set -euo pipefail

PROJECT_ROOT="/rds/general/project/hda_25-26/live/TDS/anw16/TDS_Group1"
MED_DIR="${PROJECT_ROOT}/modelling_script/08_mediation_outdated"
INPUT_DIR="${MED_DIR}/inputs"
OUT_DIR="${MED_DIR}/outputs"
LOG_DIR="${MED_DIR}/logs"

mkdir -p "${LOG_DIR}"

TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="${LOG_DIR}/male_formal_mediation_${TIMESTAMP}.log"

exec > >(tee -a "${LOG_FILE}") 2>&1

cd "${PROJECT_ROOT}"

module purge
module load tools/prod
module load R/4.3.2-gfbf-2023a
module load R-bundle-CRAN/2024.11-foss-2024a

export R_LIBS_USER=~/R/x86_64-pc-linux-gnu-library/4.3

check_file () {
  if [ ! -f "$1" ]; then
    echo "ERROR: Required file not found:"
    echo "  $1"
    exit 1
  fi
}

echo "Building male pair inputs..."
Rscript modelling_script/08_mediation_outdated/build_formal_mediation_inputs_male.R

check_file "${INPUT_DIR}/formal_mediation_config.csv"
check_file "${INPUT_DIR}/pairs_male_chunk1.csv"
check_file "${INPUT_DIR}/pairs_male_chunk2.csv"

echo "Running male chunk 1..."
Rscript modelling_script/08_mediation_outdated/run_formal_mediation_refit_chunk.R \
  "${INPUT_DIR}/pairs_male_chunk1.csv" \
  male_chunk1 \
  male

echo "Running male chunk 2..."
Rscript modelling_script/08_mediation_outdated/run_formal_mediation_refit_chunk.R \
  "${INPUT_DIR}/pairs_male_chunk2.csv" \
  male_chunk2 \
  male

check_file "${OUT_DIR}/formal_mediation_refit_split/male_chunk1/formal_mediation_results_male_chunk1.csv"
check_file "${OUT_DIR}/formal_mediation_refit_split/male_chunk2/formal_mediation_results_male_chunk2.csv"

echo "Male formal mediation refit completed successfully."
echo "Log file: ${LOG_FILE}"
