# Peter van Galen, 220205
# Use Seurat and Harmony to classify normal BM cells

library(tidyverse)
library(Seurat)
library(harmony)

setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/Github/1_Seurat_Harmony")

rm(list=ls())

# Functions
source("../Single-cell_BPDCN_Functions.R")

# Load Seq-Well BM02-2 data and make colnames more compatible
load("../../ExprStar/190514.190301.BM02-1.star/190514.190301.BM02-1.filter.RData", verbose = T)
colnames(CM.df) <- paste0(cutf(colnames(CM.df), d = "_", f = 2), "-", gsub("-", ".", cutf(colnames(CM.df), d = "_", f = 1)))
# Load 10X BM data (BM191119 & BM191227)
load("../../CellRanger/BM_agg/BM_agg.RData", verbose = T)
colnames(CM.dgm) <- gsub("-1", "-BM191119.1",
                         gsub("-2", "-BM191119.2",
                              gsub("-3", "-BM191119.3",
                                   gsub("-4", "-BM191227.1",
                                        gsub("-5", "-BM191227.2", colnames(CM.dgm))))))

# Merge. You lose some genes here: setdiff(rownames(CM.dgm), rownames(CM.df)) and setdiff(rownames(CM.df), rownames(CM.dgm))
allgenes <- intersect(rownames(CM.dgm), rownames(CM.df))
BMall.mat <- as.matrix(cbind(CM.dgm[allgenes,], CM.df[allgenes,]))


# Create Seurat object ----------------------------------------------------------------------------

bm <- CreateSeuratObject(counts = BMall.mat, project = "BM")

# Add replicate and tech information to bm object metadata
bm$replicate <- cutf(colnames(BMall.mat), d = "-", f = 2)
bm$tech <- ifelse(bm$replicate == "BM02.1", yes = "SW", no = "TenX")

# Normalize (log, transcript per 10K)
bm <- NormalizeData(bm)

# Identify variable features
bm <- FindVariableFeatures(bm)
LabelPoints(plot = VariableFeaturePlot(bm), points = head(VariableFeatures(bm), 20), repel = T, xnudge = 0, ynudge = 0) + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) + ggtitle("Variable genes")

# Scale data: (1) shift the expression of each gene, so that the mean expression across cells is 0 and (2) scales the expression of each gene, so that the variance across cells is 1; this step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
bm <- ScaleData(bm, features = rownames(bm))


# Dimension reduction with Harmony ----------------------------------------------------------------

# Run Harmony, see also http://htmlpreview.github.io/?https://github.com/immunogenomics/harmony/blob/master/docs/SeuratV3.html
bm <- RunPCA(bm, features = VariableFeatures(object = bm))
bm <- RunHarmony(object = bm, group.by.vars = "tech", reduction = "pca", plot_convergence = T)
DimPlot(bm, reduction = "harmony", group.by = "replicate") + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) + ggtitle("Harmony PCA")

# Run UMAP on Harmony distance
bm <- RunUMAP(bm, reduction = "harmony", dims = 1:50)

# Cluster. Recommended resolution 0.4-1.2 for 3,000 cells
bm <- FindNeighbors(bm, reduction = "harmony", dims = 1:50)
bm <- FindClusters(bm, resolution = 0.8)

pdf("Seurat_Harmony.pdf", width = 6, height = 6)
DimPlot(bm, reduction = "umap", group.by = "replicate") + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) + ggtitle("Harmony PCA + Seurat UMAP")
DimPlot(bm, reduction = "umap") + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) + ggtitle("Harmony PCA + Seurat clusters")
FeaturePlot(bm, features = "CD34", reduction = "umap", cols = colItay(c(1:100))) + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
FeaturePlot(bm, features = "CD14", reduction = "umap", cols = colItay(c(1:100))) + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
FeaturePlot(bm, features = "FCER1A", reduction = "umap", cols = colItay(c(1:100))) + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
FeaturePlot(bm, features = "MME", reduction = "umap", cols = colItay(c(1:100))) + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
FeaturePlot(bm, features = "MS4A1", reduction = "umap", cols = colItay(c(1:100))) + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
FeaturePlot(bm, features = "CD3D", reduction = "umap", cols = colItay(c(1:100))) + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
FeaturePlot(bm, features = "GNLY", reduction = "umap", cols = colItay(c(1:100))) + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
FeaturePlot(bm, features = "HBB", reduction = "umap", cols = colItay(c(1:100))) + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
dev.off()

# Rotate umap for easier visualization
bm.umap <- data.frame(bm@reductions$umap@cell.embeddings)
bm@reductions$umap@cell.embeddings[,1] <- bm.umap[,2]
bm@reductions$umap@cell.embeddings[,2] <- -bm.umap[,1]

# Clean up & save
bm@commands <- list()
bm@meta.data$RNA_snn_res.0.8 <- NULL
bm <- DietSeurat(bm, assays = "RNA", dimreducs = c("pca", "harmony", "umap"))
saveRDS(bm, file = "BM_Seurat_clusters.rds")



