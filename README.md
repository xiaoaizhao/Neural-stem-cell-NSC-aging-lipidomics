# Neural stem cell (NSC) aging lipidomics

This is the repository containing code for pre-processing, normalization, quantification and data presentation of lipidomics data included in the following manuscript:

>Reference: <br>**Lipidomic profiling reveals age-dependent changes in plasma membrane lipids that regulate neural stem cell aging** <br>
>Authors: Xiaoai Zhao#, Ryan M. Feitzinger@, Jeeyoon Na@, Xin Yan@, Andrew Erickson@, Kévin Contrepois@, Olivia Y. Zhou, Francesco Vallania, Mathew Ellenberger, Chloe M. Kashiwagi, Stephanie D. Gagnon, Cynthia J. Siebrand, Matias Cabruja, Gavin M. Traber, Andrew McKay, Daniel Hornburg, Purvesh Khatri, Michael P. Snyder^, Richard N. Zare^ and Anne Brunet#
>
><br>
>Preprint: https://www.biorxiv.org/content/10.1101/2022.08.18.503095v1.article-metrics

## Repository content:

[Input_Data](https://github.com/xiaoaizhao/Neural-stem-cell-NSC-aging-lipidomics/tree/main/Input_Data) - Raw input data

[Output_Data](https://github.com/xiaoaizhao/Neural-stem-cell-NSC-aging-lipidomics/tree/main/Output_Data) - Intermediary and final data output

[Scripts](https://github.com/xiaoaizhao/Neural-stem-cell-NSC-aging-lipidomics/tree/main/Scripts) - Analysis scripts for individual studies

1. [_In vitro_ #1](https://github.com/xiaoaizhao/Neural-stem-cell-NSC-aging-lipidomics/tree/main/Scripts/In_vitro_%231) - Untargeted lipidomics on aNSCs and qNSCs
2. [_In vitro_ #2](https://github.com/xiaoaizhao/Neural-stem-cell-NSC-aging-lipidomics/tree/main/Scripts/In_vitro_%232) - Untargeted lipidomics on qNSCs with genetic knockout
3. [_In vivo_ lipidomics](https://github.com/xiaoaizhao/Neural-stem-cell-NSC-aging-lipidomics/tree/main/Scripts/In_vivo_lipidomics) - Untargeted lipidomics on _in vivo_ isolated quiescent NSCs
4. [GPMV lipidomics](https://github.com/xiaoaizhao/Neural-stem-cell-NSC-aging-lipidomics/tree/main/Scripts/GPMV) - Untargeted lipidomics on giant plasma membrane vesicles (GPMVs) of quiescent NSCs
5. [DESI-MSI](https://github.com/xiaoaizhao/Neural-stem-cell-NSC-aging-lipidomics/tree/main/Scripts/DESI_MSI) - Untargetd lipidomics on the subventricular zone (SVZ) neurogenic niche _in situ_.
6. [_Mboat2_ overexpression and young plasma membrane lipid supplementation lipidomics](https://github.com/xiaoaizhao/Neural-stem-cell-NSC-aging-lipidomics/tree/main/Scripts/Mboat2_OE_and_PM_lipid_supplementation_lipidomics) - Untargeted lipidomics on qNSCs with _Mboat2_ overexpression or with young plasma membrane lipid supplementation
7. [Function scripts](https://github.com/xiaoaizhao/Neural-stem-cell-NSC-aging-lipidomics/tree/main/Scripts/Function_scripts) - Function scripts that are used in multiple different analyses
8. [Overlap and cross-study analyses](https://github.com/xiaoaizhao/Neural-stem-cell-NSC-aging-lipidomics/tree/main/Scripts/Overlap_or_across_studies) - Analyses comparing 2 or more datasets from this study


## Technical note

This repository is created using [`renv`](https://rstudio.github.io/renv/index.html) to promote and facilitate reproducibility. The R environment associated with this project is built under R version 4.4.2 and macOS 13.7.6.

To set up the environment, please install R version 4.4.2 and the package `renv` version 0.12.5. Use the `.Rproj` file to initiate this project. Then in R/Rstudio enter `renv::restore`. After the R environment is restored, all the required packages associated with this project will be installed with the correct version.



### Technical requirements:

macOS: [xcode](https://mac.install.guide/commandlinetools/4.html), [gfortran](https://stackoverflow.com/questions/35999874/mac-os-x-r-error-ld-warning-directory-not-found-for-option)

Windows: [Rtools](https://cran.r-project.org/bin/windows/Rtools/)