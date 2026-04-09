#!/bin/bash
#PBS -N ukb_impute
#PBS -q v1_small24
#PBS -l select=1:ncpus=16:ompthreads=16:mem=60gb
#PBS -l walltime=24:00:00
#PBS -j oe

cd /rds/general/project/hda_25-26/live/TDS/TDS_Group1/imputation

source ~/miniforge3/bin/activate
conda activate Renv

export LC_ALL=C
export LANG=C

~/miniforge3/envs/Renv/bin/Rscript ukb_G1_imputation.R
