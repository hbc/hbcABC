library(Seurat)
library(tidyverse)
data_dir = "data"
clusters = c(1)
seurat_obj  = "seurat_filtered_tsne.rds"
name = "subset1"
set.seed(1454944673L)

### Subclusters
seurat_init = readRDS(file.path(params$data_dir, seurat_obj))

tsne_label = FetchData(seurat_init, vars.all = c("ident", "tSNE_1", "tSNE_2"))  %>% 
    as.data.frame() %>% 
    group_by(ident) %>%
    summarise(x=mean(tSNE_1), y=mean(tSNE_2))

# check cluster are the ones we want
FetchData(seurat_init, c("tSNE_1", "tSNE_2","ident")) %>%
    rownames_to_column("cell") %>% 
    mutate(subset = cell  %in% cells) %>% 
    ggplot(aes(tSNE_1, tSNE_2)) + 
    geom_point(aes(color = subset)) +
    geom_text(data=tsne_label, aes(x, y, label = ident))

cells = FetchData(seurat_init, c("ident")) %>%
    rownames_to_column("cells") %>% 
    filter(ident  %in% clusters) %>% .[["cells"]]
seurat = SubsetData(seurat_init, cells.use = cells)

seurat =  seurat %>%
    FindVariableGenes(
        mean.function = ExpMean,
        dispersion.function = LogVMR,
        do.plot = FALSE)

seurat = seurat %>%
    RunPCA(do.print = FALSE, pcs.compute = 20) %>%
    ProjectPCA(do.print = FALSE)

# Allow for user-defined PCs to use, otherwise calculate using a PC elbow plot
pct = seurat@dr$pca@sdev / sum(seurat@dr$pca@sdev) * 100
cum = cumsum(pct)
co1 = which(cum > 90 & pct < 5)[1]
co2 = sort(which((pct[1:length(pct)-1] - pct[2:length(pct)]) >0.1),
           decreasing = T)[1] + 1
co = min(co1, co2)
data.frame(pct = pct, cumsum = cum, rank = 1:length(pct)) %>% 
    ggplot(aes(cumsum, pct, label = rank, color = rank>=co)) +
    geom_text() +
    geom_vline(xintercept = 90, color = "grey") +
    geom_hline(yintercept = min(pct[pct>5]), color = "grey")
pcs <- co

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

seurat <- RunUMAP(seurat, dims.use = 1:pcs)

saveRDS(seurat, file = file.path(data_dir,
                                     paste0("seurat_", name,"_tsne.rds")))
#rsync you data if you work in the cluster and local computer for the same data.