#PBS -l walltime=20:00:00
#PBS -l select=1:ncpus=64:ompthreads=64:mem=100gb
#PBS -N imputation

cd /rds/general/project/hda_25-26/live/TDS/fg520/TDS_Group1/imputation/

eval "$(~/anaconda3/bin/conda shell.bash hook)"
source activate phd_r

Rscript ukb_G1_imputation.R
