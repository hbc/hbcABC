# Seurat singlecell RNA-Seq clustering analysis

This is a clustering analysis workflow to be run mostly on O2 using the output from the QC which is the filtered Seurat object. This workflow incorporates Lorena's scripts available within this same `Rscripts` folder. 



To run the clustering analysis on O2, be sure to have X11 forwarding working if you want to visualize any of the images. To do this, you may need to have XQuartz running on your local machine and log onto O2 with the terminal:

```bash
ssh -XY username@o2.hms.harvard.edu
```

Then start an interactive session with extra memory and x11:

```bash
srun --pty -p interactive -t 0-12:00 --x11 --mem 64G /bin/bash
```

After starting the interactive session, load the necessary R modules:

```bash
module load gcc/6.2.0 R/3.4.1
```

# single cell clustering with seurat

library(Seurat)
library(tidyverse)
data_dir = "data"
load(file.path(data_dir, "cycle.rda"))
set.seed(1454944673L)

readRDS(file.path(data_dir, "seurat_raw.rds"))
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

# this is an automatation of the PC to be used but
# it is better to show as many as possible to the client and
# decide together if more are needed.
PCElbowPlot(seurat)
pct = seurat@dr$pca@sdev / sum(seurat@dr$pca@sdev) * 100
cum = cumsum(pct)
co1 = which(cum > 90 & pct < 5)[1]
co2 = sort(which((pct[1:length(pct)-1] - pct[2:length(pct)]) > 0.1),
           decreasing = T)[1] + 1 # last point where change of % of variation is more than 0.1%.
pcs = min(co1, co2) # change to any other number

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
#rsync you data if you work in the cluster and local computer for the same data.
