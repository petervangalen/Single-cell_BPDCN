# Peter van Galen, 220422
# Generate heatmaps of myeloid and erythroid differentiation trajectories

library(tidyverse)
library(Seurat)
library(readxl)
library(ComplexHeatmap)
library(data.table)
library(limma)

setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/07_DEG_Heatmaps")

rm(list=ls())

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")
cell_colors <- popcol.tib$hex[1:22]
names(cell_colors) <- popcol.tib$pop[1:22]
donor_colors <- popcol.tib$hex[24:41]
names(donor_colors) <- popcol.tib$pop[24:41]

# Load Seurat objects, merge
seurat_files <- list.files("../04_XV-seq", pattern = "*.rds", full.names = T)
seu.ls <- lapply(seurat_files, function(x) readRDS(x))
names(seu.ls) <- cutf(basename(seurat_files), d = "_")

# Select cells for further analysis.
celltypes.ch <- c("pDC", "cDC", "ncMono", "Mono", "ProMono", "GMP", "HSPC", "EarlyEry", "LateEry")
seu_filter.ls <- lapply(seu.ls, function(x) subset(x, CellType %in% celltypes.ch))
seu_filter.ls <- lapply(seu_filter.ls, function(x) {
  x$CellType <- factor(x$CellType, levels = celltypes.ch)
  return(x) } )


# Marker genes from BM cells ----------------------------------------------------------------------

# Define new marker genes based on filtered BM populations, similar to 200414_Annotate.R (unless you already did before)
if (file.exists("FilteredMarkerGenes.txt")) {
  markergenes.df <- read.table("FilteredMarkerGenes.txt", header = T)
} else {
  markerGenes <- FindAllMarkers(seu_filter.ls$BM, slot = "data", only.pos = T)
  markergenes.dt.ls <- lapply(split(markerGenes, f = markerGenes$cluster), function(x) data.table(x))
  markergenes.dt.ls <- lapply(markergenes.dt.ls, function(x) setorder(x, -avg_log2FC))
  markergenes.df <- do.call(cbind, lapply(markergenes.dt.ls, function(x) x$gene[1:50]))
  write.table(markergenes.df, file = "FilteredMarkerGenes.txt", sep = "\t", col.names = T, row.names = F, quote = F)
}

CurrentMarkers.ch <- unique( unlist(markergenes.df[1:10,celltypes.ch]) )


# Heatmaps per donor ------------------------------------------------------------------------------

# Generate data frame for heatmap annotations
genotyping_tables.tib <- read_excel("../4_XV-seq/XV-seq_overview.xlsx")

### Plot heatmap for each sample
for (Donor in names(seu_filter.ls)) {
#Donor <- "Pt10Dx"

# Load mutations of current donor
mutations <- genotyping_tables.tib %>% filter(Sample == Donor) %>% .$Mutation
# For MTAP (Pt10Dx and Pt10Rel), edit names
mutations <- unique(gsub("rearr.*", "rearr", mutations))

# Generate metadata table for heatmap annotation
metadata_df <- seu_filter.ls[[Donor]]@meta.data[, c("CellType", mutations), drop = F]
metadata_df <- metadata_df[sample(nrow(metadata_df)),,drop = F]
metadata_df <- metadata_df[order(metadata_df$CellType),,drop = F]

# Prevent error when plotting BM
if (Donor == "BM") {
  mutations <- "placeholder"
  metadata_df$placeholder <- "no call"
  }

# Create matrix to plot gene expression, normalize
plot_mat <- as.matrix(GetAssayData(seu_filter.ls[[Donor]], slot = "data"))[CurrentMarkers.ch,]
plot_mat <- plot_mat - rowMeans(plot_mat)

# Subset metadata and gene expression for current sample
plot_subset_mat <- plot_mat[,rownames(metadata_df)]
z.lim <- c(-2, 4)
plot_subset_mat[plot_subset_mat < z.lim[1]] <- z.lim[1]
plot_subset_mat[plot_subset_mat > z.lim[2]] <- z.lim[2]

# Define annotation objects
top_anno.ha <- HeatmapAnnotation(CellType = metadata_df$CellType,
                                 col = list(CellType = cell_colors),
                                 annotation_name_gp = gpar(fontsize = 10),
                                 border = T)

bottom_anno.ha <- HeatmapAnnotation(Mut = as.matrix(metadata_df[,mutations]),
                                    col = list(Mut = c("no call" = "#FFFFFF", "wildtype" = "#1AFF1A", "mutant" = "#4B0092")),
                                    annotation_name_gp = gpar(fontsize = 10),
                                    border = T)

# Create Heatmap object
hm <- Heatmap(as.matrix(plot_subset_mat),
              col = colItay(c(1:11))[3:11],
              cluster_rows = F,
              cluster_columns = F,
              row_names_gp = gpar(fontsize = 2),
              show_column_names = F,
              top_annotation = top_anno.ha,
              bottom_annotation = bottom_anno.ha,
              name = "Expr",
              column_title = paste0(Donor, " (", ncol(plot_subset_mat), " cells)"),
              column_title_gp = gpar(fontsize = 10),
              border = T,
              use_raster = T,
              raster_quality = 10)

pdf(paste0(Donor, "_Heatmap.pdf"), width = 8, height = 5)
print(hm)
dev.off()

}


# Combined heatmaps for supplement ----------------------------------------------------------------

# For normal
healthy_controls <- seu_filter.ls$BM
metadata_df <- healthy_controls@meta.data[,c("CellType", "replicate")]
metadata_df$Donor <- cutf(metadata_df$replicate, d = "\\.")
metadata_df <- metadata_df[order(metadata_df$CellType),,drop = F]
plot_mat <- as.matrix(GetAssayData(healthy_controls, slot = "data"))[CurrentMarkers.ch,]
pdf_name <- "Combined_Normal"

# Or: for skin-only patients
skin_only  <- merge(seu_filter.ls$Pt1Rem, seu_filter.ls[c("Pt5Dx", "Pt9Dx", "Pt10Dx", "Pt12Dx")])
metadata_df <- skin_only@meta.data[,c("CellType", "orig.ident")]
metadata_df$CellType <- factor(metadata_df$CellType, levels = celltypes.ch)
metadata_df$Donor <- factor(metadata_df$orig.ident, levels = unique(skin_only$orig.ident))
metadata_df <- metadata_df[order(metadata_df$CellType),,drop = F]
plot_mat <- as.matrix(GetAssayData(skin_only, slot = "data"))[CurrentMarkers.ch,]
pdf_name <- "Combined_Skin-only"

# Normalize gene expression matrix
plot_mat <- plot_mat - rowMeans(plot_mat)
  
# Subset metadata and gene expression for current sample
plot_subset_mat <- plot_mat[,rownames(metadata_df)]
z.lim <- c(-2, 4)
plot_subset_mat[plot_subset_mat < z.lim[1]] <- z.lim[1]
plot_subset_mat[plot_subset_mat > z.lim[2]] <- z.lim[2]

# Define annotation objects
top_anno.ha <- HeatmapAnnotation(CellType = metadata_df$CellType,
                                 Donor = metadata_df$Donor,
                                 col = list(CellType = cell_colors, Donor = donor_colors),
                                 annotation_name_gp = gpar(fontsize = 10),
                                 border = T)

# Create Heatmap object
hm <- Heatmap(as.matrix(plot_subset_mat),
              col = colItay(c(1:11))[3:11],
              cluster_rows = F,
              cluster_columns = F,
              row_names_gp = gpar(fontsize = 4),
              show_column_names = F,
              top_annotation = top_anno.ha,
              name = "Expr",
              column_title = pdf_name,
              column_title_gp = gpar(fontsize = 10),
              border = T,
              use_raster = T,
              raster_quality = 10)

pdf(paste0(pdf_name, "_Heatmap.pdf"), width = 8, height = 5)
print(hm)
dev.off()
  




