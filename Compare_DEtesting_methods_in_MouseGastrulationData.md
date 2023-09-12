
# Comparison of different DE testing methods

First comparison was done in library(MouseGastrulationData)

in R:

### Setup
```{r}
library(tidyverse)
library(miloDE)
suppressMessages(library(MouseGastrulationData))
library(scuttle)
suppressMessages(library(miloR))
suppressMessages(library(uwot))
library(scran)
library(reshape2)
#library(scWGCNA)
suppressMessages(library(Seurat))
library(viridis)
library(ggpubr)


library(BiocParallel)

ncores = 4
mcparam = MulticoreParam(workers = ncores)
register(mcparam)

# If matrix packge goes crazy again:
# remotes::install_version("Matrix", version ="1.5.4.1")

```

## 1. miloDE
### 1.1. Load data
```{r}
sce = suppressMessages(MouseGastrulationData::Tal1ChimeraData())
# downsample to few selected cell types
cts = c("Spinal cord" , "Haematoendothelial progenitors", "Endothelium" , "Blood progenitors 1" , "Blood progenitors 2")
sce = sce[, sce$celltype.mapped %in% cts]
sce$celltype.mapped[sce$celltype.mapped == "Haematoendothelial progenitors"] = "Haem. prog-s."
sce = sce[!rownames(sce) == "tomato-td" , ]
sce = logNormCounts(sce)
sce$tomato = sapply(sce$tomato , function(x) ifelse(isTRUE(x) , "Tal1_KO" , "WT"))

# for this exercise, we focus on 3000 highly variable genes (for computational efficiency)
dec.sce = modelGeneVar(sce)
hvg.genes = getTopHVGs(dec.sce, n = 3000)
sce = sce[hvg.genes , ]
rowdata = as.data.frame(rowData(sce))
rownames(sce) = rowdata$SYMBOL

set.seed(32)
umaps = as.data.frame(uwot::umap(reducedDim(sce , "pca.corrected")))
reducedDim(sce , "UMAP") = umaps

umaps = cbind(as.data.frame(colData(sce)) , reducedDim(sce , "UMAP"))
names(EmbryoCelltypeColours)[names(EmbryoCelltypeColours) == "Haematoendothelial progenitors"] = "Haem. prog-s."
cols_ct = EmbryoCelltypeColours[names(EmbryoCelltypeColours) %in% unique(umaps$celltype.mapped)]

p = ggplot(umaps , aes(x = V1 , y = V2 , col = celltype.mapped)) +
  geom_point() + 
  scale_color_manual(values = cols_ct) +
  facet_wrap(~tomato) +
  theme_bw() + 
  labs(x = "UMAP-1", y = "UMAP-2")
p

```

### 1.2. Assign neighbourhoods
```{r}
# 2.1 Estimate k -> neighbourhood size
stat_k = estimate_neighbourhood_sizes(sce, k_grid = seq(10,40,5) , order = 2, prop = 0.1 , filtering = TRUE, reducedDim_name = "pca.corrected" , plot_stat = TRUE)
stat_k
# We will use k=20, order=2 –> that returns an average neighbourhood size ~400 cells.
kable(stat_k , caption = "Neighbourhood size distribution ~ k")


# 2.2 Assign neighborhoods
set.seed(32)
sce_milo = assign_neighbourhoods(sce , k = 20 , order = 2, filtering = TRUE , reducedDim_name = "pca.corrected" , verbose = F)


nhoods_sce = nhoods(sce_milo)
# assign cell types for nhoods 
nhood_stat_ct = data.frame(Nhood = 1:ncol(nhoods_sce) , Nhood_center = colnames(nhoods_sce))
nhood_stat_ct = miloR::annotateNhoods(sce_milo , nhood_stat_ct , coldata_col = "celltype.mapped")

p = plot_milo_by_single_metric(sce_milo, nhood_stat_ct, colour_by = "celltype.mapped" , layout = "UMAP" , size_range = c(1.5,3) , edge_width = c(0.2,0.5)) + 
  scale_fill_manual(values = cols_ct , name = "Cell type")
p


```

### 1.3. DE testing
```{r}
# 1. Calculate AUC per neighbourhood
stat_auc = suppressWarnings(calc_AUC_per_neighbourhood(sce_milo , sample_id = "sample" , condition_id = "tomato", min_n_cells_per_sample = 1, BPPARAM = mcparam))
p = plot_milo_by_single_metric(sce_milo, stat_auc, colour_by = "auc" , layout = "UMAP" , size_range = c(1.5,3) , edge_width = c(0.2,0.5)) + 
  scale_fill_viridis(name = "AUC")
p

# 2. DE testing
de_stat = de_test_neighbourhoods(sce_milo ,
                             sample_id = "sample",
                             design = ~tomato,
                             covariates = c("tomato"),
                             subset_nhoods = stat_auc$Nhood[!is.na(stat_auc$auc)],
                             output_type = "SCE",
                             plot_summary_stat = TRUE,
                             layout = "UMAP", BPPARAM = mcparam , 
                             verbose = T)


```

### 1.4. Take a look at the results
```{r}
# 1 Get neighbourhood ranking by the extent of DE
stat_de_magnitude = rank_neighbourhoods_by_DE_magnitude(de_stat)

p1 = plot_milo_by_single_metric(sce_milo, stat_de_magnitude, colour_by = "n_DE_genes" , layout = "UMAP" , size_range = c(1.5,3) , edge_width = c(0.2,0.5)) + 
  scale_fill_viridis(name = "# DE genes")

p2 = plot_milo_by_single_metric(sce_milo, stat_de_magnitude, colour_by = "n_specific_DE_genes" , layout = "UMAP" , size_range = c(1.5,3) , edge_width = c(0.2,0.5)) + 
  scale_fill_viridis(name = "# specific\nDE genes" , option = "inferno")

p = ggarrange(p1,p2)
p

```

### 1.5 Get a table
```{r}
de_stat_df = convert_de_stat(de_stat)
de_stat_df = merge(de_stat_df , nhood_stat_ct , by = c("Nhood" , "Nhood_center"))
de_stat_df
write.table(de_stat_df, "/home/twatzeels/VIB_Hackathon/Hackathon_2023/results/MouseGastrulationData/de_stat_df_MouseGastrulationData.miloDE.tsv", row.names = FALSE, quote = FALSE, sep = "\t")

plots = lapply(sort(unique(modules_wgcna$cluster)) , function(cluster){
  p = plot_beeswarm_gene_set(de_stat_df, 
                             genes = modules_wgcna$gene[modules_wgcna$cluster == cluster], 
                             nhoodGroup = "celltype.mapped") + 
    ggtitle(paste0("Module ", cluster))
  return(p)
})

p
```

## 2. EdgeR on pseudobulk
```{r}
# convert to Seurat object and perform standard processing and DE
obj <- as.Seurat(sce, counts = "counts", data = "logcounts")
DimPlot(obj, reduction = "UMAP", group.by = "celltype.mapped")
obj <- RenameAssays(obj, "originalexp"="RNA")
y <- Seurat2PB(obj, sample="tomato", cluster="celltype.mapped")
dim(y)
head(y$samples, n=10L)


keep.genes <- filterByExpr(y, group=y$samples$cluster, min.count=10, min.total.count=20)
table(keep.genes)
#y <- y[keep.genes, , keep=FALSE]

# Explore y
cluster <- as.factor(y$samples$cluster)
plotMDS(y, pch=16, col=c(2:8)[cluster], main="MDS")
legend("bottomright", legend=paste0("cluster",levels(cluster)), pch=16, col=2:8, cex=0.8)

# Design matrix
donor <- factor(y$samples$sample)
design <- model.matrix(~ cluster + donor)
colnames(design) <- gsub("donor", "", colnames(design))
colnames(design)[1] <- "Int"
head(design)

y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion
plotBCV(y)

fit <- glmQLFit(y, design, robust=TRUE)
plotQLDisp(fit)

# Marker gene identification
ncls <- nlevels(cluster)
contr <- rbind( matrix(1/(1-ncls), ncls, ncls), matrix(0, ncol(design)-ncls, ncls) )
diag(contr) <- 1
contr[1,] <- 0
rownames(contr) <- colnames(design)
colnames(contr) <- paste0("cluster", levels(cluster))
contr

qlf <- list()
for(i in 1:ncls){
 qlf[[i]] <- glmQLFTest(fit, contrast=contr[,i])
 qlf[[i]]$comparison <- paste0("cluster", levels(cluster)[i], "_vs_others")
}

topTags(qlf[[1]], n=10L)
dt <- lapply(lapply(qlf, decideTestsDGE), summary)
dt.all <- do.call("cbind", dt)
dt.all

# Gather all marker genes in a table
edgeR_markers_df <- data.frame()
for (i in 1:ncls) {
  cluster_table <- qlf[[i]]$table
  cluster_table <- cluster_table[order(cluster_table$PValue, decreasing = FALSE),]
  cluster_table$Cluster <- sub("cluster", "", sub("_vs_others", "", qlf[[i]][["comparison"]]))
  edgeR_markers_df <- rbind(edgeR_markers_df, cluster_table)
}
edgeR_markers_df <- edgeR_markers_df %>% rownames_to_column("Gene")
write.table(edgeR_markers_df, "/home/twatzeels/VIB_Hackathon/Hackathon_2023/results/MouseGastrulationData/de_df_MouseGastrulationData.edgeR.tsv", row.names = FALSE, quote = FALSE, sep = "\t")

# Get top n genes for plotting
top <- 50
topMarkers <- list()
for(i in 1:ncls) {
  ord <- order(qlf[[i]]$table$PValue, decreasing=FALSE)
  up <- qlf[[i]]$table$logFC > 0
  topMarkers[[i]] <- rownames(y)[ord[up][1:top]]
}
topMarkers <- unique(unlist(topMarkers))
topMarkers


# Plot heatmap
lcpm <- cpm(y, log=TRUE)
annot <- data.frame(cluster=paste0("cluster ", cluster))
rownames(annot) <- colnames(y)
ann_colors <- list(cluster=2:6)
names(ann_colors$cluster) <- paste0("cluster ", levels(cluster))
pheatmap::pheatmap(lcpm[topMarkers, ], breaks=seq(-2,2,length.out=101), color=colorRampPalette(c("blue","white","red"))(100), scale="row", 
                   cluster_cols=TRUE, border_color="NA", fontsize_row=5,
                   treeheight_row=70, treeheight_col=70, cutree_cols=8,
                   clustering_method="ward.D2", show_colnames=FALSE, 
                   annotation_col=annot, annotation_colors=ann_colors)

```

## 3. DESeq2 and MAST (standard implementations in latest Seurat)
```{r}
head(obj@meta.data)
Idents(obj) <- obj@meta.data$celltype.mapped
MAST_markers_df <- FindAllMarkers(obj, min.pct = 0.1, test.use = "MAST")
DESeq2_markers_df <- FindAllMarkers(obj, test.use = "DESeq2")
write.table(MAST_markers_df, "/home/twatzeels/VIB_Hackathon/Hackathon_2023/results/MouseGastrulationData/de_df_MouseGastrulationData.MAST.tsv", row.names = FALSE, quote = FALSE, sep = "\t")
write.table(DESeq2_markers_df, "/home/twatzeels/VIB_Hackathon/Hackathon_2023/results/MouseGastrulationData/de_df_MouseGastrulationData.DESeq2.tsv", row.names = FALSE, quote = FALSE, sep = "\t")

```

## 4. Merge results and compare













