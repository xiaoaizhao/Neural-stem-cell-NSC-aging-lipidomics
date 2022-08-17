# Primary NSC culture #2 NSC lipidomics on LC-MS on LC-MS (2019)

This folder contains all processing script for this data set.

Please run scripts in order (0-12), which includes all pre-processing steps and downstream analysis.

### Pre-processing steps are in the following order:

* Manually annotate spike-in standards and organize into a dataframe
* Correct lipids with multiple annotation by only keeping one with the highest mscore
* For each lipid class, only keep lipid with the most abundant ion adduct (previously determined)
* Remove duplicated lipids (i.e. lipid with the exact same headgroup and side chain) by only keeping one with the highest intensity across all samples
* Calculate endogenous lipid concentration based on spike-in standard for each lipid class
* Normalization with spike-in standard
* Use median concentration normalization to normalize input material.
* Perform imputation to replace missing value

### Free fatty acid quantification
* Normalize raw FFA intensity with the normalization factor based on median concentration of complx lipid species
* Calculate FFA concentration use deuterated standard Oleic acid (d17)

### Figure panels generated:
Significant lipids with age in qNSCs (Fig. S1e)

Double bond composition heatmap (Fig. S1i)

PCA on individual lipids (Fig. 4c)

Assess KO efficiency by specific substrate and product level of target enzyme (Fig. S6b-c)