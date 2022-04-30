# Peter van Galen, 220422
# Generate heatmaps of myeloid and erythroid differentiation trajectories

library(tidyverse)
library(Seurat)
library(readxl)
library(ComplexHeatmap)
library(data.table)

setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/Github/8_Heatmaps")

rm(list=ls())

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")
cell_colors <- popcol.tib$hex[1:17]
names(cell_colors) <- popcol.tib$pop[1:17]

# Load Seurat objects, merge
bm <- readRDS("../2_Annotate/BM_Seurat_CellTypes.rds")
seurat_files <- list.files("../7_XV-seq", pattern = "*.rds", full.names = T)
seu_bpdcn.ls <- lapply(seurat_files, function(x) readRDS(x))
seu_all.ls <- c(list(bm), seu_bpdcn.ls)

# Select cells for further analysis.
celltypes.ch <- c("pDC", "cDC", "ncMono", "Mono", "ProMono", "Prog", "HSC", "EarlyE", "LateE")
seu_filter.ls <- lapply(seu_all.ls, function(x) subset(x, idents = celltypes.ch))
seu_filter.ls <- lapply(seu_filter.ls, function(x) {
  x$CellType <- factor(x$CellType, levels = celltypes.ch)
  return(x) } )
names(seu_filter.ls) <- c("BM", gsub("_Seurat_Anno.rds", "", cutf(seurat_files, "/", f = 3)))


# Marker genes from BM cells ----------------------------------------------------------------------

# Define new marker genes based on filtered BM populations, similar to 200414_Annotate.R (unless you already did before)
if (file.exists("FilteredMarkerGenes.txt")) {
    CurrentMarkers.df <- read.table("FilteredMarkerGenes.txt", header = T)
} else {
CurrentMarkers <- FindAllMarkers(seu_filter.ls[["BM"]], slot = "data", logfc.threshold = 0.25, min.pct = 0.1)
CurrentMarkers.dt.ls <- lapply(split(CurrentMarkers, f = CurrentMarkers$cluster), function(x) data.table(x))
CurrentMarkers.dt.ls <- lapply(CurrentMarkers.dt.ls, function(x) setorder(x, -avg_log2FC)) # sort by fold change
stopifnot(min(unlist(lapply(CurrentMarkers.dt.ls, function(x) x[1:50,avg_log2FC]))) > 0) # check that top 50 genes have +ve fold change
CurrentMarkers.df <- do.call(cbind, lapply(CurrentMarkers.dt.ls, function(x) x$gene[1:50]))
write.table(CurrentMarkers.df[,celltypes.ch], file = "FilteredMarkerGenes.txt", sep = "\t", col.names = T, row.names = F, quote = F)
}

CurrentMarkers.ch <- unique( unlist(CurrentMarkers.df[1:10,celltypes.ch], use.names = F) )


# Heatmap -----------------------------------------------------------------------------------------

# Generate data frame for heatmap annotations
genotyping_tables.tib <- read_excel("../7_XV-seq/FilteredCells_files.xlsx")

### Plot heatmap for each sample
for (Donor in names(seu_filter.ls)) {
#Donor <- "Pt10Dx"

# Load mutations of current donor
mutations <- genotyping_tables.tib %>% filter(Sample == Donor) %>% .$Mut 
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

