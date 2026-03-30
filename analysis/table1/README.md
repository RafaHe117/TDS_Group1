Table 1 pipeline
Overview

This folder contains the Table 1 generation workflow for the TDS Group 1 project.

The current Table 1 pipeline includes several components:

Core Table 1 outputs (main, appendix, biological)
ICD-specific Table 1 outputs
Sex-specific main Table 1 outputs
Publication-ready sex-specific main Table 1 outputs
PNG rendering of the core CSV outputs
Main execution entry point

The preferred single-entry script is:

qsub analysis/table1/run_table1_pipeline.sh
Core scripts
table1_config.R
Defines data paths, output paths, labels, variable sets, and exclusion rules.
table1_utils.R
Contains helper functions for formatting, statistical testing, missingness summaries, high-cardinality collapsing, and core table construction.
table1_run.R
Generates the six core Table 1 CSV outputs.
table1_icd_run.R
Generates ICD-specific before/after outcome-stratified tables and their PNG renders.
table1_main_bysex_before_after_run.R
Generates sex-specific main Table 1 outputs (female/male, before/after).
table1_main_bysex_before_after_publication.R
Generates publication-ready sex-specific main Table 1 outputs with relabelled rows/columns.
render_table1_png_base.R
Renders the core table1_*.csv files in analysis/table1/output/ into paginated PNGs and long PNGs.
run_table1_pipeline.sh
Master bash script that runs the full Table 1 workflow and checks that expected outputs exist.
Core outputs

Main output directory:

analysis/table1/output/

Core CSV outputs:

table1_main_outcome_before.csv
table1_main_outcome_after.csv
table1_appendix_outcome_before_missing.csv
table1_appendix_outcome_after.csv
table1_bio_bysex_before.csv
table1_bio_bysex_after.csv
Additional outputs
ICD outputs

Directory:

analysis/table1/icd_cat_output/

Includes:

ICD before/after CSVs
ICD before/after PNGs
Sex-specific main outputs

Directory:

analysis/table1/output_main_bysex_before_after/

Includes:

female before/after CSVs
male before/after CSVs
corresponding PNGs
Publication-ready sex-specific outputs

Directory:

analysis/table1/output_main_bysex_before_after_publication/

Includes:

female before/after publication CSVs
male before/after publication CSVs
corresponding PNGs
Rendered PNG outputs for core tables

Directory:

analysis/table1/output/output_rendered_png/

Includes:

paginated PNG renders of the core table1_*.csv files
single long PNG renders of the same core files
Notes
Before-imputation tables include a missingness column when requested by the pipeline.
Appendix tables use TopK + Other collapsing for high-cardinality categorical variables.
Sex-specific main tables remove sex as a row variable because sex is used to split the tables.
The single bash script is the preferred reproducible entry point for this section.
