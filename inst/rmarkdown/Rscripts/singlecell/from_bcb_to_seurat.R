## Run in local computer
specie = "mus musculus"
data_dir = "data"
bcb = "bcb_filtered.rds"
library(bcbioSingleCell)
cell_cycle_markers <- bcbioSingleCell::cellCycleMarkers[[camel(specie)]]
s_genes <- cell_cycle_markers %>%
    filter(phase == "S") %>%
    pull("geneID")
g2m_genes <- cell_cycle_markers %>%
    filter(phase == "G2/M") %>%
    pull("geneID")
save(g2m_genes, s_genes, file = file.path(data_dir,"cycle.rda"))

bcb_filtered = readRDS(file.path(data_dir, bcb))
saveRDS(rowData(bcb_filtered), file = file.path(data_dir,"rowData.rds"))
seurat <- CreateSeuratObject(
  raw.data = counts(bcb_filtered), meta.data = metrics(bcb_filtered))
saveRDS(seurat, file = file.path(data_dir,"seurat.rds"))
