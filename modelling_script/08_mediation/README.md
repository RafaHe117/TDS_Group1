# Mediation analysis pipeline

## Overview
This folder contains the mediation analysis workflow for the TDS Group 1 project.

The final pipeline is organised into six stages:

1. Define candidate biomarker and exposure inputs
2. Run exposure–biomarker screening using stability-based penalised regression
3. Fit additional multivariable final outcome models
4. Build final strict and medium shortlists
5. Refit formal mediation models for shortlisted exposure–mediator pairs
6. Generate final tables and figures

## Main execution entry point
The preferred single-entry script is:

```bash
qsub analysis/mediation/run_mediation_pipeline.sh
Core scripts
create_biomarker_candidates.R
Defines the candidate biomarker list used on the mediator side of the mediation analysis.
build_exposure_list.R
Builds the candidate exposure list / exposure universe used in the mediation workflow.
run_biomarker_models_main.R
Runs the exposure–biomarker screening model in the pooled dataset.
run_biomarker_models_sex.R
Runs the exposure–biomarker screening model separately in female and male subsets.
run_final_outcome_models.R
Fits the additional multivariable final outcome models including selected exposures, selected biomarkers, and confounders.
build_final_shortlist.R
Integrates stable links and final outcome model support to generate the final strict and medium shortlists.
run_formal_mediation_refit_chunk.R
Refits shortlisted exposure–mediator pairs formally and estimates a path, b path, direct effect, indirect effect, total effect, and bootstrap confidence intervals.
make_mediation_outputs.R
Aggregates formal mediation refit outputs into final tables and figures.
run_mediation_pipeline.sh
Master bash script that runs the full mediation pipeline in the correct order and checks that expected outputs exist.
Inputs

Main input directory:

analysis/mediation/inputs/

Key inputs include:

biomarker_candidates.csv
raw_exposure_universe.csv
selected_terms_main.csv
selected_terms_female.csv
selected_terms_male.csv
formal_mediation_config.csv
pairs_main_chunk1.csv
pairs_main_chunk2.csv
pairs_female_chunk1.csv
pairs_female_chunk2.csv
Outputs

Main output directory:

analysis/mediation/outputs/

Main output components include:

screening outputs for main, female, and male
final outcome model outputs
final shortlist outputs
formal mediation refit outputs
final mediation tables and figure outputs
Notes
The current formal refit stage uses pre-defined chunked pair files.
Male formal mediation refit is not run in the current pipeline because no strict shortlisted male pathway proceeded to that stage.
Major stages export CSV outputs so figures and tables can be regenerated without rerunning the full analysis.
The single bash script is the preferred reproducible entry point for this section.
