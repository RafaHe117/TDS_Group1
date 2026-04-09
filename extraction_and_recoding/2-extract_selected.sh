#PBS -l walltime=5:00:00
#PBS -l select=1:ncpus=1:mem=100gb
#PBS -N extraction

eval "$(~/anaconda3/bin/conda shell.bash hook)"
source activate phd_r

cd /rds/general/project/hda_25-26/live/TDS/fg520/TDS_Group1/extraction_and_recoding/scripts

ukb_path=/rds/general/project/chadeau_ukbb_folder/live/data/project_data/UKB_677583/ukb677583.csv

Rscript 2-extract_selected.R $ukb_path

conda deactivate
