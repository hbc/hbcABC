# Rscript  DE_across_condition_in_clusters.R cluster_n
set.seed(1454944673L)

args = commandArgs(trailingOnly=TRUE)
c = args[1]
cat("Doing cluster ", c, "\n")
library(Seurat)
library(zinbwave)
library(tidyverse)
data_dir = "data"
column = "condition"
BiocParallel::register(BiocParallel::SerialParam())


seurat = readRDS(file.path(data_dir, "seurat_filter_tsne.rds"))
sc = Seurat::as.SingleCellExperiment(seurat)
se = SummarizedExperiment(assays = SimpleList(counts=as.matrix(assays(sc)[["counts"]])),
                          colData = colData(sc))


se_small = se[,colData(se)[["ident"]] == c]
genes = sapply(unique(colData(se)[[column]]), function(condition){
    cells = colnames(se_small)[colData(se_small)[[column]] == condition]
    genes = rowSums(assay(se_small)[,cells]>0) / length(cells) > 0.25
    names(genes[genes])
}) %>% unlist() %>%  unique
se_small = se_small[genes,]

if (file.exists(file.path(data_dir, paste0("zb", c, ".rds")))){
    cat("\nLoad zb file\n")
    zb = readRDS(file.path(data_dir, paste0("zb", c, ".rds")))
}else{
    zb = zinbwave(se_small, K=2, epsilon=100000)
    saveRDS(zb, file = file.path(data_dir, paste0("zb", c, ".rds")))
}

library(edgeR)
dge <- DGEList(assay(zb))
dge <- calcNormFactors(dge)
design <- model.matrix(~0 + condition, data = colData(se_small))
dge$weights <- assay(zb, "weights")
dge <- estimateDisp(dge, design)
fit <- glmFit(dge, design)
lrt <- glmWeightedF(fit, coef = 2) 
# use this as needed. You can use contrast as well the same way it is used in edgeR or limma.

save(fit, lrt, file = file.path(data_dir, paste0("fit", c, ".rda")))
#rsync you data if you work in the cluster and local computer for the same data.