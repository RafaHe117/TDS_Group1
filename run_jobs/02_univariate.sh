#PBS -l walltime=04:00:00
#PBS -l select=1:ncpus=32:ompthreads=32:mem=64gb
#PBS -N univariate
#PBS -o univariate.log
#PBS -j oe

cd /rds/general/project/hda_25-26/live/TDS/fg520/TDS_Group1/modelling_script/02_univariate_finding/

eval "$(~/anaconda3/bin/conda shell.bash hook)"
source activate phd_r

export OMP_NUM_THREADS=32
export MY_NCORES=32

Rscript ukb_G1_univariate.R
