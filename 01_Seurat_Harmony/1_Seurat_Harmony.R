# Peter van Galen, 220205
# Use Seurat and Harmony to classify normal BM cells

library(tidyverse)
library(Seurat)
library(harmony)
library(ggforce)

setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/01_Seurat_Harmony")

rm(list=ls())

# Functions
source("../Single-cell_BPDCN_Functions.R")

# The following links are broken as of 241010; the count data is in "DropboxMGB/Peter van Galen/Projects/Archive/2024_Single-cell_BPDCN/AnalysisPeter"
# Load Seq-Well BM02-1 data. This already has undergone some standard QC filtering steps.
load("../../ExprStar/190514.190301.BM02-1.star/190514.190301.BM02-1.filter.RData", verbose = T)
colnames(CM.df) <- paste0(cutf(colnames(CM.df), d = "_", f = 2), "-", gsub("BM02-", "BM6.", cutf(colnames(CM.df), d = "_", f = 1)))
# Load 10x BM data. This already has undergone some standard QC filtering steps.
TenX_BM_Files <- list.files("../../CellRanger", recursive = T, pattern = "BM.*.RData", full.names = T) 
CM_merge.dgm <- NULL
for (x in TenX_BM_Files) {
  print(x)
  load(x)
  n <- 
  colnames(CM.dgm) <- gsub(1, gsub("-", ".", gsub(".RData", "", basename(x))), colnames(CM.dgm))
  CM_merge.dgm <- cbind(CM_merge.dgm, CM.dgm)
}
dim( CM_merge.dgm )

# Merge. You lose some genes here: setdiff(rownames(CM.dgm), rownames(CM.df)) and setdiff(rownames(CM.df), rownames(CM.dgm))
allgenes <- intersect(rownames(CM_merge.dgm), rownames(CM.df))
BMall.mat <- as.matrix(cbind(CM_merge.dgm[allgenes,], CM.df[allgenes,]))


# Create Seurat object ----------------------------------------------------------------------------

# Optional for loop to test different dimensions
#for (n_dimensions in c(20, 30, 40, 50)) {
n_dimensions <- 40
  
bm <- CreateSeuratObject(counts = BMall.mat, project = "BM")

# I found that increasing the QC thresholds here improves data visualization and clustering. In preprocessing, cells with >20% mitochondrial genes were already removed.
bm <- subset(bm, nCount_RNA > 2000 & nFeature_RNA > 1000)

# Add replicate and tech information to bm object metadata
bm$replicate <- cutf(colnames(bm), d = "-", f = 2)
bm$tech <- ifelse(bm$replicate == "BM6.1", yes = "SW", no = "TenX")

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
bm <- RunHarmony(object = bm, group.by.vars = c("tech", "replicate"), reduction = "pca", plot_convergence = T)

# Run UMAP on Harmony distance
bm <- RunUMAP(bm, reduction = "harmony", dims = 1:n_dimensions)

# Cluster. Recommended resolution 0.4-1.2 for 3,000 cells
bm <- FindNeighbors(bm, reduction = "harmony", dims = 1:n_dimensions)
bm <- FindClusters(bm, resolution = 0.8)

pdf(paste0("All_Seurat_Harmony_", n_dimensions, ".pdf"), width = 6, height = 6)
print( DimPlot(bm, reduction = "umap", group.by = "replicate") + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) + ggtitle("Harmony PCA + Seurat UMAP") )
print( DimPlot(bm, reduction = "umap") + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) + ggtitle("Harmony PCA + Seurat clusters") )
print( FeaturePlot(bm, features = "CD34", reduction = "umap", cols = colItay(c(1:100))) + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) )
print( FeaturePlot(bm, features = "CD14", reduction = "umap", cols = colItay(c(1:100))) + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) )
print( FeaturePlot(bm, features = "FCER1A", reduction = "umap", cols = colItay(c(1:100))) + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) )
print( FeaturePlot(bm, features = "IRF8", reduction = "umap", cols = colItay(c(1:100))) + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) )
print( FeaturePlot(bm, features = "MME", reduction = "umap", cols = colItay(c(1:100))) + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) )
print( FeaturePlot(bm, features = "MS4A1", reduction = "umap", cols = colItay(c(1:100))) + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) )
print( FeaturePlot(bm, features = "CD3D", reduction = "umap", cols = colItay(c(1:100))) + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) )
print( FeaturePlot(bm, features = "CD8A", reduction = "umap", cols = colItay(c(1:100))) + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) )
print( FeaturePlot(bm, features = "GNLY", reduction = "umap", cols = colItay(c(1:100))) + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) )
print( FeaturePlot(bm, features = "HBB", reduction = "umap", cols = colItay(c(1:100))) + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) )
dev.off()
#}

# Add UMAP coordinates to metadata. They are intentionally reversed.
bm$UMAP_1 <- bm@reductions$umap@cell.embeddings[,2]
bm$UMAP_2 <- bm@reductions$umap@cell.embeddings[,1]
# Overwrite the original with inverted coordinates, which looks nicer
bm@reductions$umap@cell.embeddings[,1] <- bm$UMAP_1
bm@reductions$umap@cell.embeddings[,2] <- bm$UMAP_2


# Dimension reduction on subsetted TNK cells ------------------------------------------------------

# Which clusters are TNK cells?
meta_data <- as_tibble(bm@meta.data, rownames = "cell") %>% mutate(TNK = seurat_clusters %in% c(0,2,3,8,13,22))
meta_data %>% ggplot(aes(x = UMAP_1, y = UMAP_2, color = seurat_clusters)) + geom_point()
meta_data %>% ggplot(aes(x = UMAP_1, y = UMAP_2, color = TNK)) + geom_point() + geom_hline(yintercept = -3.5)

# Exclude cells that of which the clustering disagrees with the UMAP coordinates
tnk_ids <- meta_data %>% filter(UMAP_2 < -3.5, TNK == T) %>% .$cell
exclude_ids <-  c(filter(meta_data, UMAP_2 > -3.5, TNK == T)$cell)
meta_data %>% ggplot(aes(x = UMAP_1, y = UMAP_2, color = cell %in% exclude_ids)) + geom_point()
bm$Doublets <- colnames(bm) %in% exclude_ids

# Visualize / cluster TNK cells
tnk <- subset(bm, seurat_clusters %in% c(0,2,3,8,13,22) & Doublets == F)
tnk <- RunPCA(tnk, features = VariableFeatures(object = tnk))
tnk <- RunHarmony(object = tnk, group.by.vars = c("tech", "replicate"), reduction = "pca", plot_convergence = T)
DimPlot(tnk, reduction = "harmony", group.by = "replicate") + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) + ggtitle("Harmony PCA")

# Run UMAP on Harmony distance
tnk <- RunUMAP(tnk, reduction = "harmony", dims = 1:10)

# Cluster. Recommended resolution 0.4-1.2 for 3,000 cells
tnk <- FindNeighbors(tnk, reduction = "harmony", dims = 1:10)
tnk <- FindClusters(tnk, resolution = 0.8)

pdf("TNK_Seurat_Harmony.pdf", width = 6, height = 6)
print( DimPlot(tnk, reduction = "umap", group.by = "replicate") + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) + ggtitle("Harmony PCA + Seurat UMAP") )
print( DimPlot(tnk, reduction = "umap") + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) + ggtitle("Harmony PCA + Seurat clusters") )
dev.off()


# Exclude additional cells -----------------------------------------------------------------------

# Cluster 20 are clearly doublets
FeaturePlot(bm, features = c("HBB", "CD3D")) + theme(aspect.ratio = 1)
as_tibble(bm@meta.data, rownames = "cell") %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = seurat_clusters == 20)) +
  geom_point()

# Cluster 24 also does not belong to any cell type as far as I can tell
as_tibble(bm@meta.data, rownames = "cell") %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = seurat_clusters == 24)) +
  geom_point()

# Exclude clusters 20 and 24 and a few cells closeby
bad_cells <- c(filter(as_tibble(bm@meta.data, rownames = "cell"), seurat_clusters == 20)$cell,
               filter(as_tibble(bm@meta.data, rownames = "cell"), seurat_clusters == 24)$cell,
               filter(as_tibble(bm@meta.data, rownames = "cell"), between(UMAP_1, -7, -2), between(UMAP_2, 0, 3))$cell)
as_tibble(bm@meta.data, rownames = "cell") %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = cell %in% bad_cells)) +
  geom_point()

# Exclude doublets
bm$Doublets <- colnames(bm) %in% c(exclude_ids, bad_cells)


# Clean up & save ---------------------------------------------------------------------------------
bm@commands <- list()
bm@meta.data$RNA_snn_res.0.8 <- NULL
bm <- DietSeurat(bm, assays = "RNA", dimreducs = c("pca", "harmony", "umap"))
saveRDS(bm, file = "BM_Seurat_clusters.rds")

tnk@commands <- list()
tnk@meta.data$RNA_snn_res.0.8 <- NULL
tnk <- DietSeurat(tnk, assays = "RNA", dimreducs = c("pca", "harmony", "umap"))
saveRDS(tnk, file = "TNK_Seurat_clusters.rds")



