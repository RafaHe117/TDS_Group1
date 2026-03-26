#!/bin/bash
#$ -N med_shortlist
#$ -cwd
#$ -o analysis/mediation/build_final_shortlist.out
#$ -e analysis/mediation/build_final_shortlist.err
#$ -l h_rt=02:00:00
#$ -l mem=12G
#$ -pe smp 2

module purge
module load tools/prod
module load R/4.3.2-gfbf-2023a
module load R-bundle-CRAN/2024.11-foss-2024a

export R_LIBS_USER=~/R/x86_64-pc-linux-gnu-library/4.3

cd /rds/general/project/hda_25-26/live/TDS/TDS_Group1

Rscript analysis/mediation/build_final_shortlist.R