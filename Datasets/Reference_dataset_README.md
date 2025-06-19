# Reference Set

Careful selection and documentation of the reference set are critically important. While the processing pipeline for the reference set is similar to that used for other datasets, the use of z-scores as a comparative metric makes transparency and consistency absolutely essential.

The reference set serves as the **baseline** for evaluating patient samples — z-scores measure how much a patient's value deviates from this reference population. Therefore, the reference set must represent **only individuals who are unaffected by inborn errors of metabolism (IMDs)** or other confounding factors that could artificially shift the distribution.

To ensure this, only data from individuals with no known IMD-related complications should be included. Any deviations or special conditions should be explicitly documented.

> For a detailed description of the technical QC steps applied to the reference set — including intra-batch and inter-batch filtering, outlier removal, and visualization — see the [Quality Control documentation](./quality_control.md).


### Examples of Potential Deviations:
1. **Medication use or transplantations** – may influence specific metabolic markers.
2. **Receiving breastfeeding** – can significantly elevate or alter the levels of various compounds, especially in infants, leading to noisy or misleading profiles. 
3. **Dietary influences** – certain foods or dietary habits may affect metabolite levels and should be accounted for.

### Age-Related Considerations
It is well known that younger children naturally have higher metabolite concentrations than adults. This reference framework accounts for age-based variation through regression modeling, as developed in the associated publication. However, age distribution within the reference set must still be well characterized and balanced to avoid skewed baselines.

