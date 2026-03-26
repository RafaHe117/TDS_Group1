#!/bin/bash
#$ -N med_male
#$ -cwd
#$ -o analysis/mediation/male.out
#$ -e analysis/mediation/male.err
#$ -l h_rt=06:00:00
#$ -l mem=32G
#$ -pe smp 4

module purge
module load tools/prod
module load R/4.3.2-gfbf-2023a

cd /rds/general/project/hda_25-26/live/TDS/TDS_Group1

Rscript analysis/mediation/run_biomarker_models_sex.R male 1