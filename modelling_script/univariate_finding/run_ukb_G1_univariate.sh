#!/bin/bash
#PBS -N ukb_G1_uni
#PBS -q v1_medium72
#PBS -l select=1:ncpus=4:mem=32gb
#PBS -l walltime=12:00:00
#PBS -j oe

set -euo pipefail

WORKDIR="/rds/general/project/hda_25-26/live/TDS/TDS_Group1/modelling_script/univariate_finding"
SCRIPT="ukb_G1_univariate.R"

cd "$WORKDIR"

mkdir -p logs

JOBID="${PBS_JOBID:-manual}"
LOG="logs/ukb_G1_univariate_${JOBID}.log"

exec > "$LOG" 2>&1

echo "========================================"
echo "JobID:   $JOBID"
echo "Host:    $(hostname)"
echo "Start:   $(date)"
echo "Workdir: $(pwd)"
echo "Script:  $SCRIPT"
echo "========================================"

source ~/miniforge3/bin/activate
conda activate Renv

export LC_ALL=C
export LANG=C

export MY_NCORES=4
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "MY_NCORES=$MY_NCORES"
echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
echo "OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS"
echo "MKL_NUM_THREADS=$MKL_NUM_THREADS"
echo "----------------------------------------"

stdbuf -oL -eL ~/miniforge3/envs/Renv/bin/Rscript --vanilla "$SCRIPT"

echo "----------------------------------------"
echo "End: $(date)"
echo "========================================"
