# _Mboat2_ overexpression and young plasma membrane lipid supplementation lipidomics

This folder contains all processing script for this data set.

### Pre-processing steps are in the following order:

* Organize deuterated standard peaks from all samples
* Append concentration to each standard
* Calculate endogenous lipid concentration based on spike-in standard for each lipid class
* Normalization with spike-in standard
* Use median concentration normalization to normalize input material
* Perform imputation to replace missing value
