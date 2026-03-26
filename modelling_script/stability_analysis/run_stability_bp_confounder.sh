#!/bin/bash
#PBS -N G1_stab_bp_conf
#PBS -q v1_medium72
#PBS -l select=1:ncpus=32:ompthreads=1:mem=120gb
#PBS -l walltime=24:00:00
#PBS -j n
#PBS -o /dev/null
#PBS -e /dev/null

set -euo pipefail

WORKDIR="/rds/general/project/hda_25-26/live/TDS/TDS_Group1/modelling_script/stability_analysis"
cd "$WORKDIR"

mkdir -p logs

JOBID="${PBS_JOBID:-manual}"
LOG_OUT="logs/stability_bp_confounder_${JOBID}.out"
LOG_ERR="logs/stability_bp_confounder_${JOBID}.err"

exec > "$LOG_OUT" 2> "$LOG_ERR"

echo "====================================="
echo "JobID: $JOBID"
echo "Host: $(hostname)"
echo "Start: $(date)"
echo "Working dir: $(pwd)"
echo "Script: stability_bp_confounder.R"
echo "====================================="

source ~/miniforge3/etc/profile.d/conda.sh
conda activate Renv

export MY_NCORES=32
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "Using cores: $MY_NCORES"

Rscript --vanilla stability_bp_confounder.R

echo "-------------------------------------"
echo "Finished: $(date)"
echo "====================================="
