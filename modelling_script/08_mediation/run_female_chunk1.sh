#!/bin/bash
#$ -N med_f1
#$ -cwd
#$ -o analysis/mediation/run_female_chunk1.out
#$ -e analysis/mediation/run_female_chunk1.err
#$ -l h_rt=04:00:00
#$ -l mem=16G
#$ -pe smp 2

module purge
module load tools/prod
module load R/4.3.2-gfbf-2023a
module load R-bundle-CRAN/2024.11-foss-2024a

export R_LIBS_USER=~/R/x86_64-pc-linux-gnu-library/4.3

cd /rds/general/project/hda_25-26/live/TDS/TDS_Group1

Rscript analysis/mediation/run_formal_mediation_refit_chunk.R \
  analysis/mediation/inputs/pairs_female_chunk1.csv \
  female_chunk1 \
  female