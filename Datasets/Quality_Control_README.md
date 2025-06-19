# Quality Control

This document outlines the Quality Control (QC) procedures applied to both the **Reference Set** and the **Dataset**, focusing on intra-batch variability, inter-batch stability, and outlier detection. The QC process ensures the reliability and consistency of data prior to downstream analysis.

These steps are included for transparency and to provide users with clear, reproducible documentation. After reading the associated publication, users can refer to this document to process their own data following the same QC strategy. Thresholds are either:
- **Suggested**: added by external users or based on experience,
- **Default**: defined within the original methodology.

---

## Reference Set

The reference set undergoes a structured multi-step QC pipeline to identify and exclude unreliable data.

### 1. Internal Standard Selection
Only **Internal Standard (ISTD)** compounds are selected from the preprocessed reference data for QC evaluation.

### 2. Intra-Batch Variability
- **Z-Score Detection**: For each compound within a batch, z-scores are calculated to detect individual outliers.  
  _Suggested threshold: 3.5_
- **Batch Size Filtering**: Batches with fewer than a minimum number of samples are flagged and excluded.  
  _Default: 10 samples_  
  Small batch sizes can disproportionately affect QC metrics and should be removed early.
- **Coefficient of Variation (CV) Evaluation**: The CV is calculated per compound per batch.  
  Batches with high CVs are flagged as unstable.  
  _Default threshold: 60%_

### 3. Inter-Batch Stability
- **Batch Means**: Mean peak areas are calculated per compound per batch over time.
- **Outlier Detection**: Batches are flagged or removed based on deviation from global mean values:
  - \>3 SD: automatic removal
  - 2–3 SD: flagged for review  
  - \>3 consecutive warnings: requires manual investigation and possibly targeted batch removal
- **Iterative Filtering**: This detection and filtering process is repeated up to a defined number of times to refine batch consistency.  
  _Suggested iterations: 5_
- **Levey-Jennings Plots**: Generated for each compound to visualize deviations using both z-score and SD-based methods.  
  These plots follow principles from Westgard Rules (not all rules implemented yet).

### 4. Manual Evaluation of Positive Controls
- All expected disease markers should be **detectable** in these controls.  
- Absence may indicate technical errors or sample degradation.

---

## Dataset

The dataset QC process applies the same logic and thresholds as established from the reference set, focused primarily on internal standard consistency.

### 1. Outlier Detection
- **Modified Z-Score Method**: Outliers are detected across all samples and compounds using the modified z-score.
- **Visualization**: Plots per compound highlight detected outliers and samples where ISTD peak area is zero.  
  _(Not described in the original paper, but highly recommended!)_

### 2. Missing Internal Standards
Each sample is checked for the presence of **all expected internal standards**.  
Samples with missing ISTDs are flagged for exclusion or further review.

