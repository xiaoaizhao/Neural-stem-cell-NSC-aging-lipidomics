# Primary NSC culture #1 NSC lipidomics on LC-MS (2017)

This folder contains all processing script for this data set.

Please run scripts in order (0-10), which includes all pre-processing steps and downstream analysis.

### Pre-processing steps are in the following order:

* Correct lipids with multiple annotation by only keeping one with the highest mscore

* For each lipid class, only keep lipid with the most abundant ion adduct (previously determined)

* Remove duplicated lipids (i.e. lipid with the exact same headgroup and side chain) by only keeping one with the highest intensity acrosss all samples

* Normalization with spike-in standard

* Average 2 technical replicates for a subgroup of samples

* Use median intensity normalization to normalize input material.

* Perform imputation to replace missing value

### Figure panels generated:

PCA plot (Fig. 1b)

Significant lipids with age in qNSCs (Fig. 1c)

Bar chart summarizing significant lipids classification (Fig. 1d)

Pie chart on lipid class composition (Fig. S1a)

Double bond composition heatmap (Fig. 1e)



