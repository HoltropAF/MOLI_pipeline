# Oligo pipeline
This repository contains resources and code to run the pipeline as described in our article.


## Required data format for the pipeline
The pipeline only accepts a (Pandas) dataframe with the following requirements/format:

| sample_ID | compound | sample_amt | peak_area | creatinine_mmol_liter | age_in_years |
|--------------|------------|------------|------------|------------|------------|
| Sample_1  | 904_α-mann  | 25 | 15000 | 1.5 | 2.2 |
| ..  | ..  | .. | .. | .. | .. |
| Sample_n  | 1071_GM1  | 35 | 25000 | 2.1 | 10.1 |

where **sample_amt** is the semi-quantitative concentration of the compound, and **peak_area** is the area under the chromatographic peak.

A test notebook can be found in this repository on how to use the pipeline (see Pipeline_test.ipynb)

