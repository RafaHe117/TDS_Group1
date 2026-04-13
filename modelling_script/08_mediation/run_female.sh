#!/bin/bash
#$ -N med_female
#$ -cwd
#$ -o modelling_script/08_mediation/female.out
#$ -e modelling_script/08_mediation/female.err
#$ -l h_rt=06:00:00
#$ -l mem=32G
#$ -pe smp 4

module purge
module load tools/prod
module load R/4.3.2-gfbf-2023a

cd /rds/general/project/hda_25-26/live/TDS/anw16/TDS_Group1

Rscript modelling_script/08_mediation/run_biomarker_models_sex.R female 0