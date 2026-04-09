#PBS -l walltime=10:00:00
#PBS -l select=1:ncpus=32:ompthreads=32:mem=80gb
#PBS -N stability_subsample
#PBS -o stability_subsample.log
#PBS -j oe
#PBS -V

cd /rds/general/project/hda_25-26/live/TDS/fg520/TDS_Group1/modelling_script/09_stability_analysis/

eval "$(~/anaconda3/bin/conda shell.bash hook)"
source activate phd_r

export OMP_NUM_THREADS=32
export MY_NCORES=32

echo "==============================="
echo "Job started at $(date)"
echo "Working directory:"
pwd
echo "==============================="

Rscript stability.R

echo "==============================="
echo "Job finished at $(date)"
echo "==============================="
