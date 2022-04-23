# Peter van Galen, 220422
# Generate heatmaps of myeloid and erythroid differentiation trajectories

library(tidyverse)
library(Seurat)
library(readxl)
library(ComplexHeatmap)
library(data.table)

setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/Github/7_Heatmaps")

rm(list=ls())

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")
cell_colors <- popcol.tib$hex[1:17]
names(cell_colors) <- popcol.tib$pop[1:17]

# Marker genes for the 17 cell types
markerGenes.df <- read.table("../2_Annotate/markerGenes.txt", header = T)

# Load Seurat objects, merge
bm <- readRDS("../2_Annotate/BM_Seurat_CellTypes.rds")
seurat_files <- list.files("../3_RandomForest", pattern = "*.rds", full.names = T)
bpdcn_seu.ls <- lapply(seurat_files, function(x) readRDS(x))
merge_seu <- merge(x = bm, y = bpdcn_seu.ls)
merge_seu$CellType <- factor(merge_seu$CellType, levels = levels(bm$CellType))
merge_seu@active.ident <- merge_seu$CellType
merge_seu$orig.ident <- factor(merge_seu$orig.ident, levels = c("BM", "Pt1Dx", "Pt1Mrd", "Pt5Dx", "Pt9Dx", "Pt10Dx", "Pt10Rel", "Pt12Dx", "Pt12Rel"))

# Functions & colors
#source("~/DropboxPartners/Pipelines/scRNAseq_SeqWell/190601_FunctionsGeneral.R")
#popcol.df <- read.xls("~/DropboxPartners/Pipelines/AuxiliaryFiles/PopCol.xlsx", sheet = 2, row.names = 1)
#heatcol <- read.xls("~/DropboxPartners/Pipelines/AuxiliaryFiles/PopCol.xlsx", header = F, sheet = 4)$V1
#markerGenes.df <- read.table("200414_Annotate/markerGenes.txt", header = T)
# Load Seurat objects, merge and reset cell types as factor
#merge.s <- merge(x = readRDS("200508_RandomForest/BM_Seurat_Predict.rds"),
#                 y = c(readRDS("../AnalysisDaniel/201214_Seurat_GoT/BPDCN628_Seurat_Genotyped.rds"),
#                       readRDS("../AnalysisDaniel/201214_Seurat_GoT/BPDCN712_Seurat_Genotyped.rds"),
#                       readRDS("200508_RandomForest/BPDCN712R_Seurat_Predict.rds")))
#merge.s$CellType <- factor(merge.s$CellType, levels = c("HSC", "Prog", "EarlyE", "LateE", "ProMono", "Mono", "ncMono", "cDC", "pDC", "ProB", "PreB", "B", "Plasma", "T", "CTL", "NK"))


#================================
# Clean up data
#================================

# Plot average expression of top marker genes. Genes that are highly expressed overall may be good contamination markers.
FiveMarkerGenes.ch <- unique( do.call(c, markerGenes.df[1:5,]) )
FiveMarkerGenes.df <- as.matrix(GetAssayData(merge_seu, slot = "data"))[FiveMarkerGenes.ch,]
# The highest expressed genes include HBB
barplot(sort(rowMeans(FiveMarkerGenes.df), decreasing = T), las = 2, ylab = "Mean expression")
abline(h = 1, col = "red")
HighMarkerGenes.num <- rowMeans(FiveMarkerGenes.df[rowMeans(FiveMarkerGenes.df) > 1,])
# Violin plots of highly expressed marker genes shows that HBB is present in some non-erythroid cells 
#VlnPlot(merge_seu, features = names(HighMarkerGenes.num), group.by = "CellType", pt.size = 0.2)

# A proportion of non-erythroid cells have high HBB expression, indicating they are doublets that escaped detection by classification, or cells that have high levels of RNA contamination. For heatmaps in Figure 2, we excluded 11.6% of HSC/Prog/Myeloid/Dendritic cells on the basis of HBB expression.
VlnPlot(merge_seu, features = "HBB", group.by = "CellType", pt.size = 0.2) + geom_hline(yintercept = 4.5)
# Including a lot of promonocytes in Patient 9 (the Seq-Well sample)
VlnPlot(subset(merge_seu, CellType == "ProMono"), features = "HBB", group.by = "orig.ident", pt.size = 0.2) + geom_hline(yintercept = 4.5)

# Remove non-erythroid cells with HBB expression
HPMD_seu <- subset(merge_seu, subset = CellType %in% setdiff(levels(merge_seu$CellType), c("EarlyE", "LateE")))
bad_cells <- subset(HPMD_seu, subset = HBB > 4.5)
ncol(bad_cells) / ncol(HPMD_seu)
# -------------------------------------------------------------------
# ^^^ IS IT REALLY NECESSARY TO EXCLUDE HBB+ NON-ERYTHROID CELLS?
# -------------------------------------------------------------------

# Select cells for further analysis.
filter_seu <- subset(merge_seu, cells = setdiff(colnames(merge_seu), colnames(bad_cells)))
celltypes.ch <- c("pDC", "cDC", "ncMono", "Mono", "ProMono", "Prog", "HSC", "EarlyE", "LateE")
filter_seu <- subset(filter_seu, idents = celltypes.ch)
filter_seu$CellType <- factor(filter_seu$CellType, levels = celltypes.ch)


#================================
# Sort cells
#================================

# Make a list of cell IDs separated by cell type, then randomize order
celltype.ls <- split(colnames(filter_seu), f = filter_seu$CellType)
celltype.ls <- lapply(celltype.ls, function(s) sample(s))

# Add randomized order of cells as meta.data to Seurat object
SortedOrder.num <- match(unlist(celltype.ls, use.names = F), colnames(filter_seu))
filter_seu <- AddMetaData(filter_seu, metadata = SortedOrder.num, col.name = "SortedOrder")


#=============================
# Marker genes from BM cells
#=============================

# Define new marker genes based on filtered BM populations, similar to 200414_Annotate.R (unless you already did before)
if (file.exists("FilteredMarkerGenes.txt")) {
    CurrentMarkers.df <- read.table("FilteredMarkerGenes.txt", header = T)
} else {
CurrentMarkers <- FindAllMarkers(subset(filter_seu, subset = orig.ident ==  "BM"), slot = "data",
                                 logfc.threshold = 0.25, min.pct = 0.1)
CurrentMarkers.dt.ls <- lapply(split(CurrentMarkers, f = CurrentMarkers$cluster), function(x) data.table(x))
CurrentMarkers.dt.ls <- lapply(CurrentMarkers.dt.ls, function(x) setorder(x, -avg_log2FC)) # sort by fold change
stopifnot(min(unlist(lapply(CurrentMarkers.dt.ls, function(x) x[1:50,avg_log2FC]))) > 0) # check that top 50 genes have +ve fold change
CurrentMarkers.df <- do.call(cbind, lapply(CurrentMarkers.dt.ls, function(x) x$gene[1:50]))
write.table(CurrentMarkers.df[,celltypes.ch], file = "FilteredMarkerGenes.txt", sep = "\t", col.names = T, row.names = F, quote = F)
}

CurrentMarkers.ch <- unique( unlist(CurrentMarkers.df[1:10,celltypes.ch], use.names = F) )


#================================
# Heatmap
#================================

# Generate data frame for heatmap annotations
#Mut.ch <- c("TET2.E1437fs*", "TET2.Q1547*", "CUX1.L911fs*", "TET2.S792*", "TET2.Q1034*", "TET2.R1216*", "TET2.H1380Y", "ASXL1.G642fs")
#metadata_df <- filter_seu@meta.data[filter_seu$SortedOrder,c("orig.ident", "CellType", paste0("Predict.", celltypes.ch), Mut.ch),]
#metadata_df[is.na(metadata_df)] <- "no call"
metadata_df <- filter_seu@meta.data[filter_seu$SortedOrder,c("orig.ident", "CellType")]

# Create matrix to plot gene expression, normalize
plot_mat <- as.matrix(GetAssayData(filter_seu, slot = "data"))[CurrentMarkers.ch,]
plot_mat <- plot_mat - rowMeans(plot_mat)

### Plot heatmap for each sample
for (Donor in levels(filter_seu$orig.ident)) {
#Donor <- "BM"

# Subset metadata and gene expression for current sample
metadata_subset_df <- metadata_df[metadata_df$orig.ident == Donor,]

plot_subset_mat <- plot_mat[,rownames(metadata_subset_df)]
z.lim <- c(-2, 4)
plot_subset_mat[plot_subset_mat < z.lim[1]] <- z.lim[1]
plot_subset_mat[plot_subset_mat > z.lim[2]] <- z.lim[2]

# Define annotation objects. For plotting continuous matrix values (e.g. Predict), see archived script "200530_Heatmaps.R".
top_anno.ha <- HeatmapAnnotation(CellType = metadata_subset_df$CellType,
                                 col = list(CellType = cell_colors),
                                 annotation_name_gp = gpar(fontsize = 10),
                                 border = T)

#bottom_anno.ha <- HeatmapAnnotation(Mut = as.matrix(metadata_subset_df[,Mut.ch]),
#                                    col = list(Mut = c("no call" = "#FFFFFF", "wildtype" = "#1AFF1A", "mutant" = "#4B0092")),
#                                    annotation_name_gp = gpar(fontsize = 10),
#                                    border = T)

# Create Heatmap object
hm <- Heatmap(as.matrix(plot_subset_mat),
              col = colItay(c(1:11))[3:11],
              cluster_rows = F,
              cluster_columns = F,
              row_names_gp = gpar(fontsize = 2),
              show_column_names = F,
              top_annotation = top_anno.ha,
#              bottom_annotation = bottom_anno.ha,
              name = "Expr",
              column_title = Donor,
              column_title_gp = gpar(fontsize = 10),
              border = T,
              use_raster = T,
              raster_quality = 10)

pdf(paste0(Donor, "_Heatmap.pdf"), width = 8, height = 5)
print(hm)
dev.off()

}

