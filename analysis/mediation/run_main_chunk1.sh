#!/bin/bash
#PBS -N med_m1
#PBS -o run_main_chunk1.out
#PBS -e run_main_chunk1.err

cd /rds/general/project/hda_25-26/live/TDS/TDS_Group1 || exit 1

module purge
module load tools/prod
module load R/4.3.2-gfbf-2023a
module load R-bundle-CRAN/2024.11-foss-2024a

export R_LIBS_USER=~/R/x86_64-pc-linux-gnu-library/4.3

Rscript analysis/mediation/run_formal_mediation_refit_chunk.R \
  analysis/mediation/inputs/pairs_main_chunk1.csv \
  main_chunk1 \
  all