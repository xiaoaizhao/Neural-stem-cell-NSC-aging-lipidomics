# DESI-MSI

This folder contains scripts for analyzing data from Desorption Electrospray Ionization Mass Spectrometry Imaging (DESI-MSI) .


* `mass_spec_deconvolution_run.R` is the script to obtain cell type-specific lipidomic profiling from SVZ of young and old mice. This script requires function scripts `data_loader.R` and `deconvolution_functions.R`.

* The other R notebooks are downstream analysis using following deconvolution
1. `In_silico_admixure_deconvolution.Rmd` - tests the accuracy of our _in silico_ deconvolution method
2. `PCA_qNSCs_DESI-MSI.Rmd` - PCA on quiescent NSC lipidomic profiles only between young and old mice
3. `sPLS-DA.Rmd` - feature selection using sPLS-DA to obtain cell type or age-specific lipids/metabolites.
4. `Aging_metabolites_heatmap_DESI-MSI.Rmd` - Heatmap on age-specific lipids/metabolites across all cell types.
5. `ANOVA_CellType_Age_interaction.Rmd` - Lipids/metbolites that show unique age vs. cell type interaction

### Figure panels generated:

Cell type-specific deconvolution from _in silico_ mixture lipidomic data (Fig. S3a)

Deconvolution from _in silico_ mixture lipidomic data with mismatched cell type proportions (Fig.S3b)

Sparse partial least squares-discriminant analysis (sPLS-DA) on DESI-MSI metabolomic profiling of SVZ cells. sPLS-DA was performed to obtain cell type-specific (Fig. S3c) and age-specific (Fig. S3d) metabolic signatures.

Heatmap on age-specific lipids/metabolites across all cell types (Fig. S3e)

Lipids/metabolites that show unique age vs. cell type interaction (Fig. S3f)

