# Peter van Galen, 220422
# Generate heatmaps of myeloid and erythroid differentiation trajectories

library(tidyverse)
library(Seurat)
library(readxl)
library(ComplexHeatmap)
#library(gdata)
#library(harmony)
#library(ggplot2)
#library(RColorBrewer)
#library(data.table)

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
merge.seu <- merge(x = bm, y = bpdcn_seu.ls)
merge.seu$CellType <- factor(merge.seu$CellType, levels = levels(bm$CellType))
merge.seu@active.ident <- merge.seu$CellType

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
FiveMarkerGenes.df <- as.matrix(GetAssayData(merge.seu, slot = "data"))[FiveMarkerGenes.ch,]
barplot(sort(rowMeans(FiveMarkerGenes.df), decreasing = T), las = 2, ylab = "Mean expression")
abline(h = 1, col = "red")
HighMarkerGenes.num <- rowMeans(FiveMarkerGenes.df[rowMeans(FiveMarkerGenes.df) > 1,])

# Violin plot separated by cell type shows that HBB is present in some non-erythroid cells
VlnPlot(merge.s, features = names(HighMarkerGenes.num), group.by = "CellType", pt.size = 0.2)

# A proportion of non-erythroid cells have high HBB expression, indicating they are doublets that escaped detection by classification, or cells that have high levels of RNA contamination. For heatmaps in Figure 2, we excluded 11.6% of HSC/Prog/Myeloid/Dendritic cells on the basis of HBB expression.
VlnPlot(merge.seu, features = "HBB", group.by = "CellType", pt.size = 0.2) + geom_hline(yintercept = 4.5)
HPMD.s <- subset(merge.seu, subset = CellType %in% c("HSC", "Prog", "ProMono", "Mono", "ncMono", "cDC", "pDC"))
badcells.s <- subset(HPMD.s, subset = HBB > 4.5)
ncol(badcells.s) / ncol(HPMD.s)

# Select cells for further analysis.
filter.s <- subset(merge.seu, cells = setdiff(colnames(merge.s), colnames(badcells.s)))
celltypes.ch <- c("pDC", "cDC", "ncMono", "Mono", "ProMono", "Prog", "HSC", "EarlyE", "LateE")
filter.s <- subset(filter.s, idents = celltypes.ch)
filter.s$CellType <- factor(filter.s$CellType, levels = celltypes.ch)


VlnPlot(subset(merge.seu, CellType == "ProMono"), features = "HBB", group.by = "orig.ident", pt.size = 0.2) + geom_hline(yintercept = 4.5)

# -------------------------------
# SHOULD I EXCLUDE ALL HBB+ NON-ERYTHROID CELLS? THIS IS WHERE I LEFT OFF 220422
# -------------------------------



#================================
# Sort cells
#================================

# Make a list of cell IDs separated by cell type, then randomize order
celltype.ls <- split(colnames(filter.s), f = filter.s$CellType)
celltype.ls <- lapply(celltype.ls, function(s) sample(s))

# Add randomized order of cells as meta.data to Seurat object
SortedOrder.num <- match(unlist(celltype.ls, use.names = F), colnames(filter.s))
filter.s <- AddMetaData(filter.s, metadata = SortedOrder.num, col.name = "SortedOrder")


#=============================
# Marker genes from BM cells
#=============================

# Define new marker genes based on filtered BM populations, similar to 200414_Annotate.R (unless you already did before)
if (file.exists(paste0(dir.ch, "CurrentMarkers.txt"))) {
    CurrentMarkers.df <- read.table(paste0(dir.ch, "CurrentMarkers.txt"), header = T)
} else {
CurrentMarkers <- FindAllMarkers(subset(filter.s, subset = orig.ident ==  "BM"), slot = "data",
                                 logfc.threshold = 0.25, min.pct = 0.1)
CurrentMarkers.dt.ls <- lapply(split(CurrentMarkers, f = CurrentMarkers$cluster), function(x) data.table(x))
CurrentMarkers.dt.ls <- lapply(CurrentMarkers.dt.ls, function(x) setorder(x, -avg_logFC))
CurrentMarkers.df <- data.frame( do.call(cbind, lapply(CurrentMarkers.dt.ls, function(x) x$gene[1:50])) )
write.table(CurrentMarkers.df[,celltypes.ch], file = paste0(dir.ch, "CurrentMarkers.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
}

CurrentMarkers.ch <- unique( unlist(CurrentMarkers.df[1:10,celltypes.ch], use.names = F) )


#================================
# Heatmap
#================================

# Generate data frame for heatmap annotations
Mut.ch <- c("TET2.E1437fs*", "TET2.Q1547*", "CUX1.L911fs*", "TET2.S792*", "TET2.Q1034*", "TET2.R1216*", "TET2.H1380Y", "ASXL1.G642fs")
MetaData.df <- filter.s@meta.data[filter.s$SortedOrder,c("orig.ident", "CellType", paste0("Predict.", celltypes.ch), Mut.ch),]
MetaData.df[is.na(MetaData.df)] <- "no call"

# Create matrix to plot gene expression, normalize
plot.mat <- as.matrix(GetAssayData(filter.s, slot = "data"))[CurrentMarkers.ch,]
plot.mat <- plot.mat - rowMeans(plot.mat)

### Plot heatmap for each sample
for (Donor in c("BM", "BPDCN628", "BPDCN712", "BPDCN712R")) {
#Donor <- "BM"
#Donor <- "BPDCN628"
#Donor <- "BPDCN712"
#Donor <- "BPDCN712R"

# Subset metadata and gene expression for current sample
MetaData.s.df <- MetaData.df[MetaData.df$orig.ident == Donor,]
plot.s.mat <- plot.mat[,rownames(MetaData.s.df)]
z.lim <- c(-2, 4)
plot.s.mat[plot.s.mat < z.lim[1]] <- z.lim[1]
plot.s.mat[plot.s.mat > z.lim[2]] <- z.lim[2]

# Define annotation objects. For plotting continuous matrix values (e.g. Predict), see archived script "200530_Heatmaps.R".
top_anno.ha <- HeatmapAnnotation(CellType = MetaData.s.df$CellType,
                                 #Predict = as.matrix(MetaData.s.df[,grepl("Predict", colnames(MetaData.s.df))]),
                                 col = list(CellType = setNames(popcol.df[celltypes.ch,"hex"], celltypes.ch)),
                                 annotation_name_gp = gpar(fontsize = 10),
                                 border = T)

bottom_anno.ha <- HeatmapAnnotation(Mut = as.matrix(MetaData.s.df[,Mut.ch]),
                                    col = list(Mut = c("no call" = "#FFFFFF", "wildtype" = "#1AFF1A", "mutant" = "#4B0092")),
                                    annotation_name_gp = gpar(fontsize = 10),
                                    border = T)

# Create Heatmap object
hm <- Heatmap(as.matrix(plot.s.mat),
              col = heatcol[3:11],
              cluster_rows = F,
              cluster_columns = F,
              row_names_gp = gpar(fontsize = 2),
              show_column_names = F,
              top_annotation = top_anno.ha,
              bottom_annotation = bottom_anno.ha,
              name = "Expr",
              column_title = Donor,
              column_title_gp = gpar(fontsize = 10),
              border = T,
              use_raster = T)

pdf(paste0(dir.ch, Donor, "_Heatmap.pdf"), width = 8, height = 5)
print(hm)
dev.off()

}





