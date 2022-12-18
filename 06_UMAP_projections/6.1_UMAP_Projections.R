# Peter van Galen, 220206
# Project BPDCN cells onto a map of normal BM contour plots

library(tidyverse)
library(Seurat)
library(readxl)
library(KernSmooth)
library(data.table)
library(gridExtra)
library(viridis)
library(MASS)

setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/06_UMAP_projections")

rm(list=ls())

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")
cell_colors <- popcol.tib$hex[1:22]
names(cell_colors) <- popcol.tib$pop[1:22]

# Load Seurat objects
seurat_files <- list.files("../04_XV-seq", pattern = "*.rds", full.names = T)
seu.ls <- lapply(seurat_files, function(x) readRDS(x))
names(seu.ls) <- cutf(basename(seurat_files), d = "_")

# Subset BM & plot named clusters
bm <- seu.ls$BM
bm.umap <- data.frame(bm@reductions$umap@cell.embeddings)
DimPlot(bm, reduction = "umap", label = T, pt.size = 0.5, cols = cell_colors) +
    theme(aspect.ratio = 1, axis.line=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(),
          axis.ticks=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(),
          panel.border=element_rect(colour = "black", fill=NA, size=1),
          plot.title = element_text(hjust = 0.5)) +
    ggtitle("Named clusters")

# Calculate healthy donors contour density
n <- 100
bins <- bkde2D(bm@reductions$umap@cell.embeddings, bandwidth = c(0.5, 0.5), gridsize = c(n, n))
meltbins <- data.table( reshape2::melt(bins$fhat) )
meltbins$Var1 <- rep(bins$x1, n)
meltbins$Var2 <- rep(bins$x2, each = n)


# Cell Projections --------------------------------------------------------------------------------

# For every patient, plot the cells on top of a bone marrow density plot. I used this for the initial submission.
for ( patient_id in names(seu.ls)[-1] ) {
#patient_id <- names(seu.ls)[2]

# Load Seurat object with coordinates etc.
seu <- seu.ls[[patient_id]]

# Data frame with coordinates
bpdcn.project.umap <- data.frame(project.umap.x = seu@meta.data$project.umap.x,
                                 project.umap.y = seu@meta.data$project.umap.y,
                                 mycol = cell_colors[as.character(seu$CellType)])

# Plot
pdf(paste0("cell_projections/", patient_id, ".pdf"), width = 6, height = 6)

# Project on bm density
print(
ggplot() +
    geom_contour(data = meltbins, mapping = aes(x = Var1, y = Var2, z = value, color = after_stat(level)), bins = 10, size = 1) +
    scale_color_gradient2(low = "grey", high = "black", guide="none") +
    ggtitle(paste0(patient_id, " Projection (", ncol(seu), " cells)")) +
    geom_point(data = bpdcn.project.umap, mapping = aes(x = project.umap.x, y = project.umap.y),
               color = bpdcn.project.umap$mycol) +
    theme(aspect.ratio = 1, axis.line=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(),
          axis.ticks=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(),
          panel.border=element_rect(colour = "black", fill=NA, size=1),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 22))
)

dev.off()

}


# Density Projections --------------------------------------------------------------------------------

# This function to calculate density is from https://slowkow.com/notes/ggplot2-color-by-density/
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

# Make three Seurat objects, then calculate density on each their metadata tables
bm <- seu.ls[["BM"]]
skin_only <- merge(seu.ls[["Pt1Rem"]], seu.ls[c("Pt5Dx", "Pt9Dx", "Pt10Dx", "Pt12Dx")])
bm_involved <- merge(seu.ls[["Pt1Dx"]], seu.ls[c("Pt10Rel", "Pt12Rel", "Pt14Dx", "Pt15Dx", "Pt16Dx")])

bm_metadata <- tibble(bm@meta.data) %>% dplyr::select(UMAP_1, UMAP_2, bm_involvement) %>%
  rename(x = UMAP_1, y = UMAP_2) %>%
  mutate(Density = get_density(x = x, y = y, n = 100))
skin_only_metadata <- tibble(skin_only@meta.data) %>% dplyr::select(project.umap.x, project.umap.y, bm_involvement) %>%
  rename(x = project.umap.x, y = project.umap.y) %>%
  mutate(Density = get_density(x = x, y = y, n = 100))
bm_involvement_metadata <- tibble(bm_involved@meta.data) %>% dplyr::select(project.umap.x, project.umap.y, bm_involvement) %>%
  rename(x = project.umap.x, y = project.umap.y) %>%
  mutate(Density = get_density(x = x, y = y, n = 100))

# From https://www.schemecolor.com/light-green-to-blue-gradient.php
density_colors <- c("#90EE90", "#67D89A", "#34BEA5", "#1EA9AC", "#1D80AF", "#1C61B1")

# Three plots side-by-side
pdf("6.1_Density_plots.pdf", width = 15, height = 5)

# Using facet_wrap (does not allow separate color scales for each)
anno_tib <- tibble(id = c("HD", "No", "Yes"),
  fraction = c(paste0(round(nrow(filter(bm_metadata, between(x, 6, 8), between(y, 9, 13)))/nrow(bm_metadata)*100, 2), "%"),
               paste0(round(nrow(filter(skin_only_metadata, between(x, 6, 8), between(y, 9, 13)))/nrow(skin_only_metadata)*100, 2), "%"),
               paste0(round(nrow(filter(bm_involvement_metadata, between(x, 6, 8), between(y, 9, 13)))/nrow(bm_involvement_metadata)*100, 2), "%")))
  
bind_rows(HD = bm_metadata, No = skin_only_metadata, Yes = bm_involvement_metadata, .id = "id") %>%
  ggplot(aes(x, y, color = Density)) +
  geom_point(size = 0.8) +
  scale_color_gradientn(colors = density_colors) +
  annotate("rect", xmin = 6, xmax = 8, ymin = 9, ymax = 13, color = "black", fill = NA) +
  geom_text(data = anno_tib, mapping = aes(x = 10, y = 11, label = fraction), color = "black", size = 7) +
  facet_wrap(~ id) +
  theme_bw() +
  theme(aspect.ratio = 1, axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))

# For separate color scales, you need to split up the data. This is not used in the paper.
p1 <- bm_metadata %>%
  ggplot(aes(x, y, color = Density)) +
  geom_point() +
  scale_color_gradientn(colors = density_colors) +
  annotate("rect", xmin = 6, xmax = 8, ymin = 9, ymax = 13, color = "black", fill = NA) +
  annotate("text", label = paste0(round(nrow(filter(bm_metadata, between(x, 6, 8), between(y, 9, 13)))/nrow(bm_metadata)*100, 2), "%"),
           x = 10, y = 11, color = "black", size = 7) +
  theme_bw() +
  theme(aspect.ratio = 1, axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))
p2 <- skin_only_metadata %>%
  ggplot(aes(x, y, color = Density)) +
  geom_point(size = 1) +
  scale_color_gradientn(colors = density_colors) +
  annotate("rect", xmin = 6, xmax = 8, ymin = 9, ymax = 13, color = "black", fill = NA) +
  annotate("text", label = paste0(round(nrow(filter(skin_only_metadata, between(x, 6, 8), between(y, 9, 13)))/nrow(skin_only_metadata)*100, 2), "%"),
           x = 10, y = 11, color = "black", size = 7) +
  theme_bw() +
  theme(aspect.ratio = 1, axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))
p3 <- bm_involvement_metadata %>%
  ggplot(aes(x, y, color = Density)) +
  geom_point(size = 1) +
  scale_color_gradientn(colors = density_colors) +
  annotate("rect", xmin = 6, xmax = 8, ymin = 9, ymax = 13, color = "black", fill = NA) +
  annotate("text", label = paste0(round(nrow(filter(bm_involvement_metadata, between(x, 6, 8), between(y, 9, 13)))/nrow(bm_involvement_metadata)*100, 2), "%"),
           x = 10, y = 11, color = "black", size = 7) +
  theme_bw() +
  theme(aspect.ratio = 1, axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))

grid.arrange(p1, p2, p3, ncol = 3)

dev.off()


# Project detection of each mutation separately ---------------------------------------------------

# Load genotyping information
genotyping_tables.tib <- read_excel("../4_XV-seq/XV-seq_overview.xlsx")
# Replace different MTAP entries with one, just as in 4.1_Add_GoT-XV_to_Seurat.R
genotyping_tables.tib$Mutation <- gsub("MTAP.rearr.*", "MTAP.rearr", genotyping_tables.tib$Mutation)
genotyping_tables.tib <- genotyping_tables.tib %>% dplyr::select(Sample, Mutation) %>% unique

for (n in 1:nrow(genotyping_tables.tib)) {
#n <- 1
Sample <- genotyping_tables.tib$Sample[n]
Mut <- genotyping_tables.tib$Mutation[n]

print(paste(Sample, Mut))

metadata_df <- seu.ls[[Sample]]@meta.data[,c("project.umap.x", "project.umap.y", Mut)]
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





