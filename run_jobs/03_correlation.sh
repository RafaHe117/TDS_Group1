#PBS -l walltime=06:00:00
#PBS -l select=1:ncpus=16:ompthreads=16:mem=64gb
#PBS -N correlation_domain
#PBS -o correlation_domain.log
#PBS -j oe

cd /rds/general/project/hda_25-26/live/TDS/fg520/TDS_Group1/modelling_script/03_correlation_analysis/

eval "$(~/anaconda3/bin/conda shell.bash hook)"
source activate phd_r

export OMP_NUM_THREADS=16
export MY_NCORES=16
export MY_CHUNK_SIZE=2

echo "Job started at $(date)"
Rscript correlation_domain.R
echo "Job finished at $(date)"
