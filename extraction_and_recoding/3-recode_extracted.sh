#PBS -l walltime=2:00:00
#PBS -l select=1:ncpus=1:mem=50gb
#PBS -N recoding

cd /rds/general/project/hda_25-26/live/TDS/fg520/TDS_Group1/extraction_and_recoding/scripts

eval "$(~/anaconda3/bin/conda shell.bash hook)"
source activate phd_r

Rscript 3-recode_variables.R

