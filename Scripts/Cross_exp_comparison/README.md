# Cross-lipidomic dataset comparison

This folder contains scripts for cross-dataset comparison.

Scripts make most sense to run in order.

### 6 Datasets included, all samples from each dataset are independent from other datasets

* Primary culture: NSC culture run on LC-MS platform. aNSCs + qNSCs from young and old cultures.
* Primary culture #2: NSC culture run on LC-MS platform with 6 different knockout conditions. qNSCs from young and old cultures.
* Lipidyzer: NSC culture run on Lipidyzer platform. aNSCs + qNSCs from young and old cultures.
* GPMV: Giant plasma membrane vesicles harvested from qNSCs of young and old cultures.
* _In vivo_: FACS-sorted qNSCs from yougn and old mice.
* DESI: Imaging mass spec from young and old SVZ niche.

### Figure panels generated:

PCA (combine Primary culture #1 + Primary culture #2) (Fig. S1c)

Effect size (on age difference) correlation between primary culture #1 vs. in vivo (Fig. 1g)

Effect size (on age difference) correlation between primary culture #1 vs. GPMV (Fig. 3c)

Effect size (on age difference) correlation between primary culture #1 vs. Lipidyzer (Fig. S1h)

Effect size (on cell type difference) correlation between primary culture #1 vs. Lipidyzer (Fig. S1g)

Effect size (on age difference) correlation between primary culture #1 vs. primary culture #2 (Fig. S1d)

Consistent double bond composition effect size (on age difference) across 3 datasets (2 primary culture + 1 in vivo datasets) (Fig. 1h)

Consistent double bond composition effect size (on age difference) across 4 datasets (2 primary culture + 1 in vivo datasets + GPMV) (Fig. S5d)

Class level change with age (Log2 fold change), compare whole cell (Primary culture #1) vs. GPMV data (Fig. S5c)

Correlation on lipid changes with age between meta-analysis across all LC-MS/MS experiment and DESI-MSI (Fig. 2e)

LION - lipid ontology enrichment analysis on Primary culture #1 dataset between young and old qNSCs (Fig. 3a)

LION - lipid ontology enrichment analysis between Primary culture #1 dataset and GPMV lipidomics (Fig. S5b)

All lipidomic aging features - including side chain unsaturation features and individual lipid features from meta-analysis (Fig. 3e)

Lipidomic aging score on all aging features in KO samples from Primary culture #2 (Fig. 4d)

Lipidomic aging score on side chain unsaturation features alone in KO samples from Primary culture #2 (Fig. S6d)

Lipidomic aging score on individual lipid features alone in KO samples from Primary culture #2 (Fig. S6e)



