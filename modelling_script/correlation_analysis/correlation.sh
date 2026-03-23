#!/bin/bash
#PBS -N corr_domain_24c
#PBS -q v1_medium72
#PBS -l select=1:ncpus=24:ompthreads=24:mem=64gb
#PBS -l walltime=24:00:00
#PBS -j n
#PBS -o /dev/null
#PBS -e /dev/null

set -euo pipefail

############################################
# 1. Working directory
############################################
WORKDIR="/rds/general/user/rh725/projects/hda_25-26/live/TDS/TDS_Group1/modelling_script/correlation_analysis"
SCRIPT="correlation_domain.R"

cd "$WORKDIR" || { echo "Cannot cd to $WORKDIR"; exit 1; }

mkdir -p logs

############################################
# 2. Unique logs
############################################
JOBID="${PBS_JOBID:-manual}"
LOG_OUT="logs/correlation_domain_${JOBID}.out"
LOG_ERR="logs/correlation_domain_${JOBID}.err"

exec > "$LOG_OUT" 2> "$LOG_ERR"

echo "========================================"
echo "JobID:      ${JOBID}"
echo "Host:       $(hostname)"
echo "Start:      $(date)"
echo "Workdir:    $(pwd)"
echo "Script:     ${SCRIPT}"
echo "User:       $(whoami)"
echo "========================================"

############################################
# 3. Environment
############################################
source ~/miniforge3/bin/activate
conda activate Renv

export LC_ALL=C
export LANG=C

echo "R location: $(which R)"
echo "Rscript:    $(which Rscript)"
echo "----------------------------------------"

############################################
# 4. Thread / parallel control
############################################
export MY_NCORES=24
export MY_CHUNK_SIZE=3

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1

echo "PBS_NCPUS=${PBS_NCPUS:-NA}"
echo "MY_NCORES=$MY_NCORES"
echo "MY_CHUNK_SIZE=$MY_CHUNK_SIZE"
echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
echo "OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS"
echo "MKL_NUM_THREADS=$MKL_NUM_THREADS"
echo "----------------------------------------"

############################################
# 5. Run R script
############################################
stdbuf -oL -eL ~/miniforge3/envs/Renv/bin/Rscript --vanilla "$SCRIPT"

############################################
# 6. End
############################################
echo "----------------------------------------"
echo "Finished at: $(date)"
echo "========================================"
