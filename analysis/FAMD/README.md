# FAMD pipeline

## Overview
This folder contains the pseudo-FAMD / mixed-data PCA workflow for the TDS Group 1 project.

The workflow covers two datasets:
1. `ukb_G1_imputed_final.rds`
2. `split_imputed_data/ukb_G1_train_imputed.rds`

For each dataset, three pseudo-FAMD versions are run:
- `pseudo_famd1a`: baseline main, no alcohol
- `pseudo_famd1b`: baseline + alcohol sensitivity
- `pseudo_famd2`: biological supplementary

## Main execution entry point
The preferred single-entry script is:

```bash
qsub analysis/FAMD/run_famd_pipeline.sh
Core scripts
check_famd_variance_only.R
Checks the variance explained by the first principal components for the baseline pseudo-FAMD setup.
pseudo_famd_run.R
Runs the full pseudo-FAMD workflow for both datasets and all three analysis versions.
run_famd_pipeline.sh
Master bash script that runs the variance check, runs the main pseudo-FAMD workflow, and checks that expected outputs exist.
Outputs

Main output directory:

analysis/FAMD/FAMD_output/

Outputs are generated separately for:

analysis/FAMD/FAMD_output/final/
analysis/FAMD/FAMD_output/train/

For each dataset and each pseudo-FAMD version, the workflow exports:

scree plot PNG
individuals-by-outcome PNG
individuals-by-sex PNG
contribution PC1 PNG
contribution PC2 PNG

In addition, figure-ready CSV outputs are exported so plots can be regenerated without rerunning the full analysis:

variance explained CSV
scores CSV
contribution PC1 CSV
contribution PC2 CSV
Notes
The output folders for each dataset are overwritten on each run.
The CSV outputs are intended to preserve the numbers needed to regenerate the figures later.
The single bash script is the preferred reproducible entry point for this section.
