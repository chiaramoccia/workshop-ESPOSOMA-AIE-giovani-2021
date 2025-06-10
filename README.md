




# 🧪 Workshop: Introduction to the Exposome Approach in Epidemiology Using R

📅 **AIE (Italian Epidemiology Association) Young Researchers Conference – Autumn 2021**  
👩‍🔬 Chiara Moccia & Antonio D'Errico  
🏛️ University of Turin – Department of Medical Sciences, Cancer Epidemiology Unit
This repository contains R scripts and materials from a hands-on workshop introducing the exposome approach in epidemiological research, with practical examples using the `rexposome` package.

---

## Requirements

- R (version 4.1.0 used during the workshop, year 2021)  
- R packages:
  - `rexposome` (version 1.15.0)  
  - `ggplot2`  
  - `BiocManager`  

To install `rexposome`:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("rexposome")
```


---

## 🎯 Workshop Objectives

- Introduce the concept of the **exposome** and its relevance in epidemiology
- Provide a practical guide to using the **`rexposome`** R package
- Explore data preprocessing, visualization, dimensionality reduction (PCA), clustering, and Exposome-Wide Association Studies (ExWAS)



## 📁 Dataset Structure

The analysis uses three input datasets:

- **Description data**: metadata for exposures
- **Exposure data**: quantitative/categorical exposure variables
- **Phenotype data**: individual-level outcome/phenotype information

## 📊 Analysis Workflow

1. Load and inspect the datasets
2. Create an `ExposomeSet` object
3. Handle missing data and test exposure normality
4. Explore exposure distributions by family and phenotype
5. Perform PCA on exposures
6. Visualize PCA results and cluster individuals
7. Assess correlation between exposures (inter- and intra-family)
8. Run univariate (ExWAS) and multivariate (ENET) association analyses

## 📌 Notes

- Some plots require an interactive R session (e.g. `win.graph()` on Windows)
- Example data paths should be updated to match your local setup
- The script includes both visual and statistical outputs

## 🔄 Based on Original Work by ISGlobal

This workshop was adapted from materials developed by the **ISGlobal Bioinformatics Research Group in Epidemiology (BRGE)**.  
The original tutorial and code are available at:  
👉 [https://isglobal-brge.github.io/rexposome/](https://isglobal-brge.github.io/rexposome/)

We acknowledge their work and thank them for making the resources openly available.

---

Feel free to clone or adapt the code for educational or research purposes.
