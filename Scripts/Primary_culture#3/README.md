# Primary NSC culture #3 NSC lipidomics on LC-MS/MS

This folder contains all processing script for this data set.

Please run scripts in order (0-12), which includes all pre-processing steps.

### Pre-processing steps are in the following order:

* Extract deuterated spike-in standards from all samples
*  Match Progenesis extracted peaks with LipidSearch annotation on peaks
* Remove multiple annotation from positive mode data and clean based on ions used for quantification for each class.
* Remove multiple annotation from negative mode data and clean based on ions used for quantification for each class.
* Remove duplicated lipids (lipids that have the identical headgroup and side chain identification) by only keeping one lipid that has the highest average intensity.
* Match deuterated spike-in standards with their respective concentration, in preparation for lipid quantification.
* Calculate endogenous lipid concentration based on spike-in standard for each lipid class
* Normalization with spike-in standard
* Use median concentration normalization to normalize input material
* Perform imputation to replace missing value
