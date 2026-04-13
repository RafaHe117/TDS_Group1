#!/bin/bash
#$ -N med_final
#$ -cwd
#$ -o modelling_script/08_mediation/final_outcome_models.out
#$ -e modelling_script/08_mediation/final_outcome_models.err
#$ -l h_rt=06:00:00
#$ -l mem=24G
#$ -pe smp 4

module purge
module load tools/prod
module load R/4.3.2-gfbf-2023a
module load R-bundle-CRAN/2024.11-foss-2024a

export R_LIBS_USER=~/R/x86_64-pc-linux-gnu-library/4.3

cd /rds/general/project/hda_25-26/live/TDS/anw16/TDS_Group1

Rscript modelling_script/08_mediation/run_final_outcome_models.R