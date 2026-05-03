# _KANSL1_: A Potential Suppressor of Metastasis in Muscle-Invasive Bladder Cancer
---
This code includes the analysis for my final-year BSc project. <br> 

It is an investigation into _KANSL1_ mutations in The Cancer Genome Atlas (TCGA) Muscle-Invasive Bladder Cancer (MIBC) Cohort. <br>

---
## Useful Links

[Project Brief](https://asmasonomics.github.io/courses/BSc_dissertation_2025 "Dr Andrew Mason's Project Brief") <br> 
[Cleaned Data - Google Drive](https://drive.google.com/drive/folders/1inlKf9uXYdE9skNrJz5yBoElZmteu8CG "Cleaned Data - Google Drive") <br> 
[cBioPortal - TCGA MIBC Dashboard](https://www.cbioportal.org/study/summary?id=blca_tcga_pub_2017 "TCGA MIBC cohort dashboard - cBioPortal") <br> 
[Robertson et al. (2017), Introduction to TCGA MIBC cohort](https://www.sciencedirect.com/science/article/pii/S0092867417310565?via%3Dihub "MIBC Cohort Introduction") <br>

---
## Contact Information

This project was carried out by Caleb. *Email address:* jfn510@york.ac.uk <br> 

---
## Installation Instructions

Prerequisite software: RStudio: [Download RStudio Desktop Here](https://posit.co/download/rstudio-desktop).
The following packages should be installed:

```
install.packages(c("tidyverse", "data.table", "googledrive", "ggrepel"))

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("DESeq2", "fgsea", "maftools", "BiocParallel"))

```

---
## Contents of notable files
| File | Description|
| -----|------------|
| analysis.R | main analysis code|
| results/KANSL1_mutvsWT_DEA_results.csv | DEA results for _KANSL1_ mutant versus _KANSL1_ wildtype tumours |
| results/FGSEA_KANSL1_muts_vs_WT.csv | GSEA results for_KANSL1_ mutant versus _KANSL1_ wildtype tumours |
| results/KANSL1_mutvsWT_DEA_results_excl_BaSq.csv | DEA results for _KANSL1_ mutant versus _KANSL1_ wildtype tumours, excluding Ba/Sq tumours |
| results/FGSEA_KANSL1_muts_vs_WT_excl_BaSq.csv | GSEA results for_KANSL1_ mutant versus _KANSL1_ wildtype tumours, excluding Ba/Sq tumours |


---
## That's the end of this README!
Thanks for reading to the end - let me end with one of the author's favourite jokes:
> "How do you turn a duck into a soul singer?
> 
> Put it in the oven until its bill withers."
