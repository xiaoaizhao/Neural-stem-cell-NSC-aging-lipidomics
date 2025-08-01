# Primary NSC culture #1 NSC lipidomics on LC-MS/MS

This folder contains all processing script for this data set.

Please run scripts in order (0-9), which includes all pre-processing steps.

### Pre-processing steps are in the following order:

* Correct lipids with multiple annotation by only keeping one with the highest mscore

* For each lipid class, only keep lipid with the most abundant ion adduct (previously determined)

* Remove duplicated lipids (i.e. lipid with the exact same headgroup and side chain) by only keeping one with the highest intensity acrosss all samples

* Normalization with spike-in standard

* Average 2 technical replicates for a subgroup of samples

* Use median intensity normalization to normalize input material

* Perform imputation to replace missing value



