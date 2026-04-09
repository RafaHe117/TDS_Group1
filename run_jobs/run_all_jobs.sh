#!/bin/bash

echo "=================================="
echo "Submitting TDS pipeline jobs"
echo "Started at $(date)"
echo "=================================="

cd /rds/general/project/hda_25-26/live/TDS/fg520/TDS_Group1/run_jobs || exit

echo "Submitting imputation split..."
qsub 01_imputation_split.sh

echo "Submitting univariate..."
qsub 02_univariate.sh

echo "Submitting correlation..."
qsub 03_correlation.sh

echo "Submitting stability main..."
qsub 04_stability.sh

echo "Submitting stability BP..."
qsub 04_stability_bp.sh

echo "Submitting stability BP confounder..."
qsub 04_stability_bp_confounder.sh

echo "Submitting stability sex..."
qsub 04_stability_sex.sh

echo "Submitting Table 1..."
qsub 05_table1.sh

echo "Submitting FAMD..."
qsub 06_famd.sh

echo "Submitting mediation..."
qsub 07_mediation.sh

echo "=================================="
echo "All jobs submitted at $(date)"
echo "=================================="
