# miRSCAPE - Inferring miRNA expression in single-cell clusters

## About

Micro-RNAs (miRNA) are critical in development, homeostasis, and diseases, including cancer. However, our understanding of miRNA function at cellular resolution is thwarted by the inability of the standard single cell RNA-seq protocols to capture miRNAs. Here we introduce a machine learning tool -- miRSCAPE -- to infer miRNA expression in a sample from its RNA-seq profile. We establish miRSCAPE's accuracy separately in 10 tissues comprising ~10,000 tumor and normal bulk samples and demonstrate that miRSCAPE accurately infers cell type-specific miRNA activities (predicted vs observed fold-difference correlation ~ 0.81) in two independent datasets where miRNA profiles of specific cell types are available (HEK-GBM, Kidney-Breast-Skin). When trained on human hematopoietic cancers, miRSCAPE can identify active miRNAs in 8 hematopoietic cell lines in mouse with a reasonable accuracy (auROC = 0.67). Finally, we apply miRSCAPE to infer miRNA activities in scRNA clusters in Pancreatic and Lung cancers, as well as in 56 cell types in the Human Cell Landscape (HCL). Across the board, miRSCAPE recapitulates and provides a refined view of known miRNA biology. miRSCAPE promises to substantially expand our understanding of gene regulatory networks at cellular resolution.

## How to run

R 4.1.0 is required

Depends : xgboost, Seurat


# Step 0: Initialization: Data and library load

```
setwd("~/miRSCAPE-main")
rm(list = ls())
require(Seurat)
source("~/miRSCAPE-main/code/miRSCAPE.R")
```

Rows are genes, columns are patient values. In order to properly run the code, gene names as row names must be provided. 
```
example = readRDS('example data/example.rds')
pdac_mirna <- read.delim("~/miRSCAPE-main/example data/pdac_mirna.txt", row.names=1)
pdac_mrna <- read.delim("~/miRSCAPE-main/example data/pdac_mrna.txt", row.names=1)
```

To specify the cell types in scRNA, define the cell clusters of interest. 
```
clustId = c('Ductal cell type 1')
```

## Step 1. scRNA preparation
To prepare the data which is stored in Seurat object, `modifySeuratObject` function should be used.

ScRNA data preparation for only on the specific cell clusters
```
denem = modifySeuratObject(pbmc = example, clusterId = clustId)

```

ScRNA data preparation on all cell clusters
```
denem = modifySeuratObject(pbmc = example)
```

## Step 2. Bulk data preparation
Please use bulkTransform function. `justNorm` parameter is used only to log-transform the data. Don't scale the miRNA data, please set `justNorm` parameter to `TRUE` in miRNA data.

```
bulkk_mirna = bulkTransform(pdac_mirna, justNorm = TRUE)
bulkk_mrna = bulkTransform(pdac_mrna)
```

## Step 3. Predict miRNA
To predict the miRNAs, `miRSCAPE` function should be used. 
```
pred = miRSCAPE(bulkmRNA = bulkk_mrna, bulkmiRNA = bulkk_mirna, scmRNA = denem, nrnds = 20)

```
