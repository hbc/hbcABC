library(Seurat)
library(tidyverse)
data_dir = "data"
load(file.path(data_dir, "cycle.rda"))

readRDS(file.path(data_dir, "seurat.rds"))
seurat <- seurat %>%  NormalizeData(
  normalization.method = "LogNormalize",
  scale.factor = 10000)
seurat =  seurat %>%
  FindVariableGenes(
    mean.function = ExpMean,
    dispersion.function = LogVMR,
    do.plot = FALSE)
seurat = seurat %>%
  ScaleData(model.use = "linear")
VariableGenePlot(seurat, do.text = F)

seurat = CellCycleScoring(
  seurat,
  g2m.genes = g2m_genes,
  s.genes = s_genes)
seurat = RunPCA(
  seurat,
  pc.genes = c(s_genes, g2m_genes),
  do.print = FALSE)
saveRDS(seurat, file = file.path(data_dir, "seurat_pre_regress.rds"))

vars_to_regress = c("nUMI", "S.Score", "G2M.Score")
seurat <- ScaleData(seurat, vars.to.regress = vars_to_regress)
seurat = RunPCA(
  seurat,
  pc.genes = c(s_genes, g2m_genes),
  do.print = FALSE)

seurat = seurat %>%
  RunPCA(do.print = FALSE) %>%
  ProjectPCA(do.print = FALSE)

pc_use <- PCElbowPlot(seurat)
pcs = 15 # change to any other number
seurat <- FindClusters(
  seurat,
  dims.use = 1:pcs,
  force.recalc = TRUE,
  print.output = TRUE,
  resolution = 1,
  save.SNN = TRUE)
seurat <- RunTSNE(
  seurat,
  dims.use = 1:pcs,
  do.fast = TRUE)
saveRDS(seurat, file = file.path(data_dir, "seurat_tsne.rds"))
# dev.off() # if you want to save the figures generated.
# Use this objects in local computer to make the report with figures.
