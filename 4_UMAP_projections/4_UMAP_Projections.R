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

# Calculate healthy donors contour density
n <- 100
bins <- bkde2D(bm@reductions$umap@cell.embeddings, bandwidth = c(0.5, 0.5), gridsize = c(n, n))
meltbins <- data.table( reshape2::melt(bins$fhat) )
meltbins$Var1 <- rep(bins$x1, n)
meltbins$Var2 <- rep(bins$x2, each = n)


# Cell Projections --------------------------------------------------------------------------------

# List files with BPDCN scRNA-seq data
seurat_files <- list.files("../3_RandomForest/", pattern = "*.rds")

# Plot every patient
for ( seu_file in seurat_files ) {
#seu_file <- seurat_files[1]

# Load Seurat object with coordinates etc.
seu <- readRDS(paste0("../3_RandomForest/", seu_file))

# Data frame with coordinates
bpdcn.project.umap <- data.frame(project.umap.x = seu@meta.data$project.umap.x,
                                 project.umap.y = seu@meta.data$project.umap.y,
                                 mycol = cell_colors[seu$CellType])

# Plot
patient_id <- cutf(seu_file, d = "_")

pdf(paste0("cell_projections/", patient_id, ".pdf"), width = 6, height = 6)
par(mar=c(4,4,4,4))

# Project on bm density
print(
ggplot() +
    geom_contour(data = meltbins, mapping = aes(x = Var1, y = Var2, z = value, color = after_stat(level)), bins = 10, size = 1) +
    scale_color_gradient2(low = "grey", high = "black", guide="none") +
    ggtitle(paste0(patient_id, " Projection (", ncol(seu), " cells)")) +
    geom_point(data = bpdcn.project.umap, mapping = aes(x = project.umap.x, y = project.umap.y), color = bpdcn.project.umap$mycol) +
    theme(aspect.ratio = 1, axis.line=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), panel.border=element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), plot.title = element_text(hjust = 0.5, size = 22))
)

dev.off()

}


# Mutation Projections ----------------------------------------------------------------------------

# Load genotyping information
genotyping_tables.tib <- read_excel("../7_XV-seq/FilteredCells_files.xlsx")
# Replace different MTAP primers with one, just as in 7_Add_GoT-XV_to_Seurat.R
genotyping_tables.tib$Mut <- gsub("MTAP.rearr.*", "MTAP.rearr", genotyping_tables.tib$Mut)
genotyping_tables.tib <- genotyping_tables.tib %>% select(Sample, Mut) %>% unique

# Load Seurat files with genotyping information
seurat_files <- list.files("../7_XV-seq", pattern = "*.rds", full.names = T)
seu_bpdcn.ls <- lapply(seurat_files, function(x) readRDS(x))
names(seu_bpdcn.ls) <- gsub("_Seurat_Anno.rds", "", cutf(seurat_files, d = "/", f = 3))

for (n in 1:nrow(genotyping_tables.tib)) {
#n <- 1
Sample <- genotyping_tables.tib$Sample[n]
Mut <- genotyping_tables.tib$Mut[n]

print(paste(Sample, Mut))

metadata_df <- seu_bpdcn.ls[[Sample]]@meta.data[,c("project.umap.x", "project.umap.y", Mut)]
colnames(metadata_df) <- c("project.umap.x", "project.umap.y", "mutated")
metadata_df$mutated <- gsub("mutant", "#4B0092", gsub("wildtype", "#1AFF1A", gsub("no call", "grey", metadata_df$mutated)))
metadata_df$mutated <- factor(metadata_df$mutated, levels = c("#4B0092", "#1AFF1A", "grey"))
metadata_df <- metadata_df[rev(order(metadata_df$mutated)),]

pdf(paste0("mut_projections/", n, "_", Sample, "_", gsub("/|:", "-", Mut), ".pdf"), width = 6, height = 6)
par(mar=c(4,4,4,4))
print(
ggplot() +
  geom_contour(data = meltbins, mapping = aes(x = Var1, y = Var2, z = value, color = after_stat(level)), bins = 10, size = 1) +
  scale_color_gradient2(low = "grey", high = "black", guide="none") +
  ggtitle(paste0(Mut, " in ", Sample, " (", nrow(metadata_df), " cells)")) +
  geom_point(data = metadata_df, mapping = aes(x = project.umap.x, y = project.umap.y), color = metadata_df$mutated) +
  theme(aspect.ratio = 1, axis.line=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(),
        axis.ticks=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(),
        panel.border=element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14))
)
dev.off()
}

# RPS24 and PRKDC ---------------------------------------------------------------------------------

Sample <- "Pt10Dx"
Mut <- "RPS24.chr10:79795273:T/C"
Mut <- "PRKDC.chr8:48713410:C/T"
Mut <- "RAB9A.3pUTR"
Mut <- "MTAP.rearr"

metadata_df <- seu_bpdcn.ls[[Sample]]@meta.data[,c("project.umap.x", "project.umap.y", Mut)]
colnames(metadata_df) <- c("project.umap.x", "project.umap.y", "mutated")
metadata_df$mutated <- gsub("mutant", "#4B0092", gsub("wildtype", "#1AFF1A", gsub("no call", "grey", metadata_df$mutated)))
metadata_df$mutated <- factor(metadata_df$mutated, levels = c("#4B0092", "#1AFF1A", "grey"))
#metadata_df <- metadata_df[rev(order(metadata_df$mutated)),]
metadata_df <- metadata_df[sample(nrow(metadata_df)),]
metadata_df <- metadata_df[metadata_df$mutated != "grey",]

pdf(paste0(Sample, "_", gsub("/|:", "-", Mut), ".pdf"), width = 6, height = 6)
par(mar=c(4,4,4,4))
print(
  ggplot() +
    geom_contour(data = meltbins, mapping = aes(x = Var1, y = Var2, z = value, color = after_stat(level)), bins = 10, size = 1) +
    scale_color_gradient2(low = "grey", high = "black", guide="none") +
    ggtitle(paste0(Mut, " in ", Sample, " (", nrow(metadata_df), " cells)")) +
    geom_point(data = metadata_df, mapping = aes(x = project.umap.x, y = project.umap.y), color = metadata_df$mutated) +
    theme(aspect.ratio = 1, axis.line=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(),
          axis.ticks=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(),
          panel.border=element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 14))
)
dev.off()
