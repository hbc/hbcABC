# Run clustering and markers quickly
library(dplyr)
# load("data/cycle.rda")
library(SingleCellExperiment)
library(Seurat)

## General Clustering
se=readRDS("data/se_filtered.rds")
# sub  = se[,grepl("group1", colData(se)[["sample"]])] # this for subsampling
seurat <- CreateSeuratObject(
    raw.data = assay(se), meta.data = as.data.frame(colData(se)))
prefix = "all"
source("scripts/clustering.R")


seurat=readRDS("data/all_seurat_tsne.rds")
prefix = "all"
source("scripts/makers.R")
