#Mouse Gastrulation public data analysis, https://bioconductor.org/packages/release/data/experiment/html/MouseGastrulationData.html

#Lets try to plot the DE results sent by Tijs and Tim:

#In R:

library(data.table)
library(UpSetR)
library(tidyr)
library(dplyr)  
library(ggplot2)

deseq2 <- fread("de_df_MouseGastrulationData.DESeq2_CaseVsCon.tsv")
#edgeR <- fread("de_df_MouseGastrulationData.edgeR.tsv") #Still not between conditions, ignore it for now
MAST <- fread("de_df_MouseGastrulationData.MAST_CaseVsCon.tsv")
miloDE <- fread("de_stat_df_MouseGastrulationData.miloDE.tsv")

processed_deseq2 <- deseq2[,c("Gene","celltype","p_val_adj","avg_log2FC")]
#processed_edgeR <- edgeR[,c("Gene","Cluster","PValue","logFC")]
processed_MAST <- MAST[,c("Gene","celltype","p_val_adj","avg_log2FC")]
processed_miloDE <- miloDE[,c("gene","celltype.mapped","pval_corrected_across_genes","logFC")] #here consider btween pval_corrected between genes or nhoods.

#Do FDR and bonf corrections for edgeR:

#processed_edgeR$p_corrected_FDR <- p.adjust(processed_edgeR$PValue, method = "fdr", n = length(processed_edgeR$PValue))
#processed_edgeR$p_corrected_Bonf <- p.adjust(processed_edgeR$PValue, method = "bonferroni", n = length(processed_edgeR$PValue))  

#6670 FDR hits vs 284 Bonf hits -> use Bonf maybe.

#processed2_edgeR <- processed_edgeR[,c(c("Gene","Cluster","p_corrected_Bonf","logFC"))]

processed_deseq2$category <- "DESeq2"
#processed2_edgeR$category <- "edgeR"
processed_MAST$category <- "MAST"
processed_miloDE$category <- "miloDE"

colnames(processed_deseq2) <- c("gene","cluster","p_val_adj","logFC","category")
#colnames(processed2_edgeR) <- c("gene","cluster","p_val_adj","logFC","category")
colnames(processed_MAST) <- c("gene","cluster","p_val_adj","logFC","category")
colnames(processed_miloDE) <- c("gene","cluster","p_val_adj","logFC","category")

#Potential problem: in miloDE there are multiple occurences of gene-cluster pairs, likely because of neighborhood parameter. Choose the most significant to plot.

processed2_miloDE <- as.data.table(processed_miloDE %>% group_by(gene,cluster) %>% top_n(-1, p_val_adj))

processed2_miloDE <- as.data.table(processed_miloDE %>% group_by(gene,cluster) %>% top_n(1, logFC))

merged_DE <- rbind(processed_deseq2,processed_MAST,processed2_miloDE)

signif_merged_DE <- subset(merged_DE, p_val_adj < 0.05)

summary(as.factor(signif_merged_DE$cluster))

#For all clusters

allSignif_gene_sets <- subset(merged_DE, gene %in% unique(signif_merged_DE$gene))

allSignif_gene_sets <- allSignif_gene_sets[,c("gene","p_val_adj","category","cluster")]

wide_allSignif_gene_sets <- allSignif_gene_sets %>% spread(category, p_val_adj)

wide_allSignif_gene_sets[is.na(wide_allSignif_gene_sets)] <- 1

wide_allSignif_gene_sets$DeSeq2 <- ifelse(wide_allSignif_gene_sets$DESeq2 < 0.05, 1, 0)
#wide_allSignif_gene_sets$edgeR <- ifelse(wide_allSignif_gene_sets$edgeR < 0.05, 1, 0)
wide_allSignif_gene_sets$MAST <- ifelse(wide_allSignif_gene_sets$MAST < 0.05, 1, 0)
wide_allSignif_gene_sets$miloDE <- ifelse(wide_allSignif_gene_sets$miloDE < 0.05, 1, 0)

wide_allSignif_gene_sets$clusterOrdered <- factor(wide_allSignif_gene_sets$cluster, level = c("Endothelium", "Spinal cord", "Haem. prog-s.", "Blood progenitors 1", "Blood progenitors 2"))

#UpsetR plot below, but cannot be combined with facet_grid

plot_all_DE <- upset(wide_allSignif_gene_sets, nsets = 4, number.angles = 30, point.size = 3.5, line.size = 2, mainbar.y.label = "DE Method Hit Intersections", sets.x.label = "Hits Per DE Method", text.scale = c(1.3, 1.3, 1, 1, 2, 0.75),order.by = "freq")

plot_all_DE_wPanels <- plot_all_DE

png("all_comparisonDE_final_withoutEdgeR.png", height=6, width=9, res=600, units="in")
plot_all_DE_wPanels
dev.off()

