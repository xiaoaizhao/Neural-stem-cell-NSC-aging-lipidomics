# Giant plasma membrane vesicles (GPMV) lipidomics on LC-MS/MS

This folder contains all processing script for this data set.

Please run scripts in order (0-7), which includes all pre-processing steps.

### Pre-processing steps are in the following order:

* Manually annotate spike-in standards and organize into a dataframe
* Correct lipids with multiple annotation by only keeping one with the highest mscore
* For each lipid class, only keep lipid with the most abundant ion adduct (previously determined)
* Remove duplicated lipids (i.e. lipid with the exact same headgroup and side chain) by only keeping one with the highest intensity across all samples
* Calculate endogenous lipid concentration based on spike-in standard for each lipid class
* Normalization with spike-in standard
* Use median concentration normalization to normalize input material
* Perform imputation to replace missing value
