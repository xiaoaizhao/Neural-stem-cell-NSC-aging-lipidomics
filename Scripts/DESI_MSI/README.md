# DESI-MSI

This folder contains scripts for analyzing data from Desorption Electrospray Ionization Mass Spectrometry Imaging (DESI-MSI) .


* `mass_spec_deconvolution_run.R` is the script to obtain cell type-specific lipidomic profiling from SVZ of young and old mice. This script requires function scripts `data_loader.R` and `deconvolution_functions.R`.
1. `In_silico_admixure_deconvolution.Rmd` - tests the accuracy of our _in silico_ deconvolution method

2. `OPLS-DA_in_DESI.R` - OPLS-DA analysis on identifying cell type-specific metabolites between aNSCs and qNSCs after deconvolution. As well as OPLS-DA on age-specific metabolites in qNSCs between young and old mice.

3. `DESI_effect_size_on_annotated_lipids.R` - Comparison of age-related effect size changes in lipids from DESI-MSI, compared to _in vitro_ summary and _in vivo_ qNSC lipidomics.

   
