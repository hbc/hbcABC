# single cell find markers with seurat

library(Seurat)
library(tidyverse)
data_dir = "data"
set.seed(1454944673L)
seurat_object = "seurat_tsne.rds"

readRDS(seurat, file = file.path(data_dir, seurat_object))
markers = FindAllMarkers(object = seurat, only.pos = TRUE, min.pct = 0.25,
                         thresh.use = 0.25)
saveRDS(markers, file.path(data_dir, "markers.rds"))
#rsync you data if you work in the cluster and local computer for the same data.
