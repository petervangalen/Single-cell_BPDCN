# Peter van Galen, 220205
# Now that we have a UMAP coordinates and clusters (1_Seurat_Harmony), we can classify the cell types

library(tidyverse)
library(Seurat)
library(readxl)
library(data.table)

setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/Github/2_Annotate")

rm(list=ls())

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")
cell_colors <- popcol.tib$hex[1:17]
names(cell_colors) <- popcol.tib$pop[1:17]

# Load BM Data
bm <- readRDS("../1_Seurat_Harmony/BM_Seurat_clusters.rds")
# Make Seurat's count matrix easier to access (log, transcript per 10K)
expr.mat <- as.matrix( GetAssayData(object = bm, slot = "data") )
#table(round(colSums(exp(expr.mat)-1)))
#table( cutf(colnames(expr.mat), d = "-", f = 2) )


# Plot signature scores ---------------------------------------------------------------------------

# Average gene expression for scoreSignature function (to save time)
expr.mat.mean <- rowMeans(expr.mat)

# Load signatures from AML paper (Van Galen, Hovestadt ... Bernstein, 2019, Cell)
aml_signs.dt <- fread("markerGenes_AML_project.txt")
# Subset for gene names that are also in the current dataset
signs.ls <- lapply(aml_signs.dt, intersect, rownames(expr.mat))

pdf("1_NormalBM_SignatureScores.pdf", width = 6, height = 6)
# Plot signature scores
for (n in names(signs.ls)) {
    message(n)
    
    # Calculate signature score
    scores <- scoreSignature(CM = expr.mat, signatures = signs.ls[[n]], CM.mean = expr.mat.mean, verbose = TRUE)
    
    # Add as meta.data to Seurat object and plot
    bm <- AddMetaData(bm, metadata = scores, col.name = n)
    print(
    FeaturePlot(bm, features = n, reduction = "umap", cols = colItay(c(1:100))) + theme(aspect.ratio = 1, axis.line=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), panel.border=element_rect(colour = "black", fill=NA, size=1))
    )
    
    # Remove score meta.data
    bm@meta.data[,n] <- NULL
}
dev.off()

# To look at individual markers:
#FeaturePlot(bm, features = "FCGR3A", reduction = "umap", cols = colItay(c(1:100))) + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))


# Name clusters -----------------------------------------------------------------------------------

pdf("2_ClusterNames.pdf", width = 6, height = 6)

# Plot Seurat clusters
DimPlot(bm, reduction = "umap", label = TRUE, pt.size = 0.5) +
  theme(aspect.ratio = 1, axis.line=element_blank(), axis.text.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title.x=element_blank(),
        axis.title.y=element_blank(), panel.border=element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(hjust = 0.5)) + ggtitle("Seurat clusters")

# Name clusters
bm <- RenameIdents(object = bm,
                   "0" = "T",
                   "1" = "EarlyE",
                   "2" = "Mono",
                   "3" = "LateE",
                   "4" = "CTL",
                   "5" = "ProMono",
                   "6" = "LateE",
                   "7" = "PreB",
                   "8" = "Mono", # previously nk
                   "9" = "NK", # previously prog
                   "10" = "Prog", # previously mono
                   "11" = "cDC",
                   "12" = "B",
                   "13" = "HSC",
                   "14" = "EarlyE",
                   "15" = "ProB",
                   "16" = "ncMono",
                   "17" = "Doublets",
                   "18" = "pDC",
                   "19" = "ProB",
                   "20" = "Plasma")

# Add CellType metadata
bm$CellType <- factor(bm@active.ident, levels = c("HSC", "Prog", "EarlyE", "LateE", "ProMono", "Mono", "ncMono", "cDC", "pDC", "ProB", "PreB", "B", "Plasma", "T", "CTL", "NK", "Doublets"))
bm@active.ident <- bm$CellType

# Plot named clusters
DimPlot(bm, reduction = "umap", label = TRUE, pt.size = 0.5, cols = cell_colors) +
  theme(aspect.ratio = 1, axis.line=element_blank(), axis.text.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title.x=element_blank(),
        axis.title.y=element_blank(), panel.border=element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(hjust = 0.5)) + ggtitle("Named clusters")

dev.off()

# Save
bm@commands <- list()
bm@meta.data$seurat_clusters <- NULL
saveRDS(bm, file = "BM_Seurat_CellTypes.rds")


# Frequency barplot -------------------------------------------------------------------------------

CellTypes.ls <- split(bm@active.ident, f = cutf(names(bm@active.ident), d = "-|\\.", f = 2))

# Plot
BMfreq.mat <- do.call(cbind, lapply(CellTypes.ls, table))
BMfreq.mat <- BMfreq.mat[!grepl("Doublets", rownames(BMfreq.mat)),]
BMnorm.mat <- sweep(BMfreq.mat, 2, colSums(BMfreq.mat), "/")*100

pdf("3_BM_CellTypes.pdf", width = 5, height = 6)
par(mar = c(8,4,8,12), xpd = T)

barplot(BMnorm.mat[nrow(BMnorm.mat):1,], col = rev(cell_colors[-17]), xaxt = "n", ylab = "Population frequency (%)", border = NA)
axis(side = 1, at = seq(1,ncol(BMnorm.mat))*1.2-0.5, labels = colnames(BMnorm.mat), las = 2)
legend(x = ncol(BMnorm.mat)*1.2+0.5, y = 100, legend = rownames(BMnorm.mat), fill = cell_colors[-17], bty = "n", border = NA)

dev.off()


# Marker genes ------------------------------------------------------------------------------------

markerGenes <- FindAllMarkers(bm, slot = "data", logfc.threshold = 0.25, min.pct = 0.1)
markergenes.dt.ls <- lapply(split(markerGenes, f = markerGenes$cluster), function(x) data.table(x))
markergenes.dt.ls <- lapply(markergenes.dt.ls, function(x) setorder(x, -avg_log2FC))
# Check that all top 50 genes have a positive fold change
stopifnot(min(unlist(lapply(markergenes.dt.ls, function(x) x[1:50,avg_log2FC]))) > 0)
markergenes.df <- do.call(cbind, lapply(markergenes.dt.ls, function(x) x$gene[1:50]))

markergenesOrder.df <- markergenes.df[,c("HSC", "Prog", "EarlyE", "LateE", "ProMono", "Mono", "ncMono", "cDC", "pDC", "ProB", "PreB", "B", "Plasma", "T", "CTL", "NK", "Doublets")]

write.table(markergenesOrder.df, file = "markerGenes.txt", sep = "\t", col.names = T, row.names = F, quote = F)


# Average expression ------------------------------------------------------------------------------

# Do it in base R
celltypes.ls <- split(rownames(bm@meta.data), bm$CellType)
mean_expr.mat <- do.call(cbind, lapply(celltypes.ls, function(x) rowMeans(expr.mat[,x])))
head(mean_expr.mat)
# Check
mean_expr.mat["HBD",]

write.table(mean_expr.mat, file = "mean_expr.txt", quote = F, sep = "\t")



