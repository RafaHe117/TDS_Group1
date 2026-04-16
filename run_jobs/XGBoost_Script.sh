#!/bin/bash
#PBS -l walltime=06:00:00
#PBS -l select=1:ncpus=16:mem=32gb
#PBS -N XGBoost

cd /rds/general/project/hda_25-26/live/TDS/anw16/TDS_Group1/modelling_script/07_XGBoost # ADAPT THIS 

module load Python/3.12.3-GCCcore-13.3.0

export PYTHONPATH="/rds/general/user/anw16/home/my_python_libs:$PYTHONPATH"

python XGBoost_Full_Script.py
