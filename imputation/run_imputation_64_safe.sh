#!/bin/bash
#PBS -N ukb_impute
#PBS -q v1_medium72
#PBS -l select=1:ncpus=64:ompthreads=64:mem=100gb
#PBS -l walltime=24:00:00
#PBS -j n
#PBS -o /dev/null
#PBS -e /dev/null

set -euo pipefail

############################################
# 1️⃣ Force correct working directory
############################################
WORKDIR="/rds/general/project/hda_25-26/live/TDS/TDS_Group1/imputation"
cd "$WORKDIR" || { echo "Cannot cd to $WORKDIR"; exit 1; }

mkdir -p logs

############################################
# 2️⃣ Unique log files per job
############################################
JOBID="${PBS_JOBID:-manual}"
LOG_OUT="logs/ukb_impute_${JOBID}.out"
LOG_ERR="logs/ukb_impute_${JOBID}.err"

exec > "$LOG_OUT" 2> "$LOG_ERR"

echo "========================================"
echo "JobID:   ${JOBID}"
echo "Host:    $(hostname)"
echo "Start:   $(date)"
echo "Workdir: $(pwd)"
echo "========================================"

############################################
# 3️⃣ Load environment
############################################
source ~/miniforge3/bin/activate
conda activate Renv

export LC_ALL=C
export LANG=C

############################################
# 4️⃣ Thread control (prevents thread conflicts)
############################################
export MY_NCORES=64
export OMP_NUM_THREADS=64
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "MY_NCORES=$MY_NCORES"
echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
echo "----------------------------------------"

############################################
# 5️⃣ Run R (line-buffered for real-time logs)
############################################
stdbuf -oL -eL ~/miniforge3/envs/Renv/bin/Rscript --vanilla ukb_G1_imputation.R

############################################
# 6️⃣ End message
############################################
echo "----------------------------------------"
echo "Finished at: $(date)"
echo "========================================"
