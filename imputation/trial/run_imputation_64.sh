#!/bin/bash
#PBS -N ukb_impute
#PBS -q v1_medium72
#PBS -l select=1:ncpus=64:ompthreads=64:mem=100gb
#PBS -l walltime=24:00:00
#PBS -j oe

set -euo pipefail

# go to submit directory (or fallback)
cd "${PBS_O_WORKDIR:-/rds/general/project/hda_25-26/live/TDS/TDS_Group1/imputation}"

mkdir -p logs

# ---- REAL-TIME LOGS (this expands PBS_JOBID correctly) ----
JOBID="${PBS_JOBID:-manual}"
exec > "logs/ukb_impute_${JOBID}.out" 2> "logs/ukb_impute_${JOBID}.err"

# ---- environment ----
source ~/miniforge3/bin/activate
conda activate Renv

export LC_ALL=C
export LANG=C

# ---- IMPORTANT: your cluster does NOT set PBS_NCPUS ----
# We pass cores explicitly to R via MY_NCORES:
export MY_NCORES=64

# Some libs respect OMP threads too
export OMP_NUM_THREADS="$MY_NCORES"

echo "JobID:   ${PBS_JOBID:-NA}"
echo "Host:    $(hostname)"
echo "Date:    $(date)"
echo "Workdir: $(pwd)"
echo "MY_NCORES=$MY_NCORES"
echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
echo "----------------------------------------"

# ---- run R ----
~/miniforge3/envs/Renv/bin/Rscript ukb_G1_imputation_V2.R

echo "----------------------------------------"
echo "Finished at: $(date)"
