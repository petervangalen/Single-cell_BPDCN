# Peter van Galen, 220206
# Project BPDCN cells onto a map of normal BM contour plots

library(tidyverse)
library(Seurat)
library(readxl)
library(KernSmooth)
library(data.table)

setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/Github/4_UMAP_projections")

rm(list=ls())

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")
cell_colors <- popcol.tib$hex[1:17]
names(cell_colors) <- popcol.tib$pop[1:17]

# Load BM Data (Seurat object)
bm <- readRDS("../2_Annotate/BM_Seurat_CellTypes.rds")
bm.umap <- data.frame(bm@reductions$umap@cell.embeddings)

# Plot named clusters
DimPlot(bm, reduction = "umap", label = T, pt.size = 0.5, cols = cell_colors) +
    theme(aspect.ratio = 1, axis.line=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(),
          axis.title.x=element_blank(), axis.title.y=element_blank(), panel.border=element_rect(colour = "black", fill=NA, size=1),
          plot.title = element_text(hjust = 0.5)) +
    ggtitle("Named clusters")


# Visualize ---------------------------------------------------------------------------------------

# List files with scRNA-seq data
seurat_files <- list.files("../3_RandomForest/", pattern = "*.rds")

# Plot every patient
for ( seu_file in seurat_files ) {
#seu_file <- seurat_files[1]
patient_id <- cutf(seu_file, d = "_")

# Load Seurat object with coordinates etc.
seu <- readRDS(paste0("../3_RandomForest/", seu_file))

# Plot
pdf(paste0(patient_id, "_visualizations.pdf"), width = 6, height = 6)
par(mar=c(4,4,4,4))

# Project umap ----------------------------------
# Calculate healthy donors contour density
n <- 100
bins <- bkde2D(bm@reductions$umap@cell.embeddings, bandwidth = c(0.5, 0.5), gridsize = c(n, n))
meltbins <- data.table( reshape2::melt(bins$fhat) )
meltbins$Var1 <- rep(bins$x1, n)
meltbins$Var2 <- rep(bins$x2, each = n)

bpdcn.project.umap <- data.frame(project.umap.x = seu@meta.data$project.umap.x,
                                 project.umap.y = seu@meta.data$project.umap.y,
                                 mycol = cell_colors[seu$CellType])

# Project on bm density
print(
ggplot() +
    geom_contour(data = meltbins, mapping = aes(x = Var1, y = Var2, z = value, color = after_stat(level)), bins = 10, size = 1) +
    scale_color_gradient2(low = "grey", high = "black", guide="none") +
    ggtitle(paste0(patient_id, " Projection (", ncol(seu), " cells)")) +
    geom_point(data = bpdcn.project.umap, mapping = aes(x = project.umap.x, y = project.umap.y), color = bpdcn.project.umap$mycol) +
    theme(aspect.ratio = 1, axis.line=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), panel.border=element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), plot.title = element_text(hjust = 0.5, size = 22))
)

# Seurat UMAP -----------------------------------
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu, features = rownames(seu))
seu <- RunPCA(seu, features = VariableFeatures(object = seu))
seu <- RunUMAP(seu, reduction = "pca", dims = 1:50)
print(
  DimPlot(seu, reduction = "umap", label = F, pt.size = 0.5, cols = cell_colors) +
    theme(aspect.ratio = 1, axis.line=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(),
          axis.ticks=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(),
          panel.border=element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste(patient_id, "Seurat UMAP"))
)

dev.off()

}





