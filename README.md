# TDS_Group1

## Instructions
### STEP 1: Produce missingness heatmap
Run `ukb_G1_build_raw_dataset.R` (Update path if necessary)

Outputs are saved in the `figure/` directory.

### STEP 2: Run analysis (part 1)
Run `run_jobs\run_all_jobs.sh` (Update path if necessary)

This bash script runs the following analyses in sequential order:

1) Imputation Split
2) Univariate Analysis
3) Correlation Analysis
4) Stability Analysis
5) Table 1 Generation
6) FAMD
7) Mediation Analysis

Note: Mediation Analysis may take a few hours to complete.

### STEP 3: Run analysis (part 2)
3.1 Run `modelling_script\05_pca\pca_script-copy.R` for PCA

3.2 Run `modelling_script\09_stability_analysis\incremental_plot.R` to produce incremental plots

3.3 Run the following for XGBoost: (Update path if necessary)

- `modelling_script\07_XGBoost\XGBoost All.ipynb`
- `modelling_script\07_XGBoost\XGBoost Female.ipynb`
- `modelling_script\07_XGBoost\XGBoost Male.ipynb`


All outputs (plots, tables, CSV files, etc.) from STEP 2 and 3 are saved in their respective folders in `modelling_script`.