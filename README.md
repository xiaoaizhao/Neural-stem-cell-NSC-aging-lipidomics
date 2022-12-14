# Neural stem cell (NSC) aging lipidomics

This is the repository containing code for pre-processing, normalization, quantification and data presentation of lipidomics data included in the following manuscript:

>Reference: <br>**Lipidomic profiling reveals age-dependent changes in complex plasma membrane lipids that regulate neural stem cell aging** <br>
>Authors: Xiaoai Zhao, Xin Yan#, Kévin Contrepois#, Francesco Vallania, Mathew Ellenberger, Chloe M. Kashiwagi, Stephanie D. Gagnon, Cynthia J. Siebrand, Matias Cabruja, Gavin M. Traber, Andrew McKay, Daniel Hornburg, Purvesh Khatri, Michael P. Snyder^, Richard N. Zare^ and Anne Brunet <br>
>Preprint: https://www.biorxiv.org/content/10.1101/2022.08.18.503095v1.article-metrics

## Repository content:

[Input_Data](https://github.com/xiaoaizhao/Neural-stem-cell-NSC-aging-lipidomics/tree/main/Input_Data) - Raw input data

[Output_Data](https://github.com/xiaoaizhao/Neural-stem-cell-NSC-aging-lipidomics/tree/main/Output_Data) - Intermediary and final data output

[Scripts](https://github.com/xiaoaizhao/Neural-stem-cell-NSC-aging-lipidomics/tree/main/Scripts) - Analysis scripts for individual studies

1. [Primary culture#1](https://github.com/xiaoaizhao/Neural-stem-cell-NSC-aging-lipidomics/tree/main/Scripts/Primary_culture%231) - Untargeted lipidomics on activated and quiescent NSCs

2. [Primary culture#2](https://github.com/xiaoaizhao/Neural-stem-cell-NSC-aging-lipidomics/tree/main/Scripts/Primary_culture%232) - Untargeted lipidomics on quiescent NSCs with genetic knockout

3. [In vivo lipidomics](https://github.com/xiaoaizhao/Neural-stem-cell-NSC-aging-lipidomics/tree/main/Scripts/In_vivo_lipidomics) - Untargeted lipidomics on _in vivo_ isolated quiescent NSCs

4. [GPMV lipidomics](https://github.com/xiaoaizhao/Neural-stem-cell-NSC-aging-lipidomics/tree/main/Scripts/GPMV) - Untargeted lipidomics on giant plasma membrane vesicles (GPMVs) of quiescent NSCs

5. [Lipidyzer lipidomics](https://github.com/xiaoaizhao/Neural-stem-cell-NSC-aging-lipidomics/tree/main/Scripts/Lipidyzer) - Targeted lipidomics on activated and quiescent NSCs

6. [DESI-MSI](https://github.com/xiaoaizhao/Neural-stem-cell-NSC-aging-lipidomics/tree/main/Scripts/DESI_MSI) - Untargetd lipidomics on the subventricular zone (SVZ) neurogenic niche _in situ_.

7. [Function scripts](https://github.com/xiaoaizhao/Neural-stem-cell-NSC-aging-lipidomics/tree/main/Scripts/Function_scripts) - Function scripts that are used in multiple different scripts

[Figure Panels](https://github.com/xiaoaizhao/Neural-stem-cell-NSC-aging-lipidomics/tree/main/Figure_Panels) - Figures included in the manuscript



## Technical note

This repository is created using [`renv`](https://rstudio.github.io/renv/index.html) to promote and facilitate reproducibility. The R environment associated with this project is built under R version 4.0.2 and macOS 12.3.1.

To set up the environment, please install R version 4.0.2 and the package `renv`. Use the `.Rproj` file to initiate this project. Then in R/Rstudio enter `renv::restore`. After the R environment is restored, all the required packages associated with this project will be installed with the correct version.



### Technical requirements:

macOS: [xcode](https://mac.install.guide/commandlinetools/4.html), [gfortran](https://stackoverflow.com/questions/35999874/mac-os-x-r-error-ld-warning-directory-not-found-for-option)

Windows: [Rtools](https://cran.r-project.org/bin/windows/Rtools/)