# Peter van Galen, 220205
# Now that we have a UMAP coordinates and clusters (1_Seurat_Harmony), we can classify the cell types

library(tidyverse)
library(Seurat)
library(readxl)
library(data.table)
library(limma)

setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/02_Annotate")

rm(list=ls())

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")
my_colors <- popcol.tib$hex
names(my_colors) <- popcol.tib$pop

# Load data
bm <- readRDS("../01_Seurat_Harmony/BM_Seurat_clusters.rds")
tnk <- readRDS("../01_Seurat_Harmony/TNK_Seurat_clusters.rds")
# Make Seurat's count matrix easier to access (log, transcript per 10K)
bm_expr.mat <- as.matrix( GetAssayData(object = bm, slot = "data") )
tnk_expr.mat <- as.matrix( GetAssayData(object = tnk, slot = "data") )
#table(round(colSums(exp(bm_expr.mat)-1)))
#table( cutf(colnames(bm_expr.mat), d = "-", f = 2) )


# All cells: plot signature scores ----------------------------------------------------------------

# Average gene expression for scoreSignature function (to save time)
bm_expr.mat.mean <- rowMeans(bm_expr.mat)

# Load signatures from AML paper (Van Galen, Hovestadt ... Bernstein, 2019, Cell) and Yoke Seng Lee (unpublished). These signature files are not synced to Github to prevent confusion; use markerGenes.txt instead. Email Peter van Galen with questions.
signs_aml.dt <- fread("markerGenes_AML_project.txt")
colnames(signs_aml.dt) <- paste0("AML_", colnames(signs_aml.dt))
signs_yoke.dt <- fread("markerGenes_Yoke_all.txt")
colnames(signs_yoke.dt) <- paste0("Yoke_", colnames(signs_yoke.dt))
signs.dt <- cbind(signs_aml.dt, signs_yoke.dt)
# Subset for gene names in current dataset
signs.ls <- lapply(signs.dt, intersect, rownames(bm_expr.mat))
names(signs.ls) <- gsub(" |-", ".", names(signs.ls))

pdf("1_NormalBM_SignatureScores.pdf", width = 6, height = 6)

DimPlot(bm) + theme(aspect.ratio = 1)

# Plot cell cycle score
bm <-  CellCycleScoring(bm, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
FeaturePlot(bm, features = "G2M.Score", cols = colItay(c(1:100))) + theme(aspect.ratio = 1)

# Plot signature scores
for (n in names(signs.ls)) {
    message(n)
    
    # Calculate signature score. This is quite similar to Seurat's function AddModuleScore.
    scores <- scoreSignature(CM = bm_expr.mat, signatures = signs.ls[[n]], CM.mean = bm_expr.mat.mean, verbose = TRUE)
    
    # Add as meta.data to Seurat object and plot
    bm <- AddMetaData(bm, metadata = scores, col.name = n)
    print(
    FeaturePlot(bm, features = n, reduction = "umap", cols = colItay(c(1:100))) +
      theme(aspect.ratio = 1, axis.line=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(),
            axis.ticks=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(),
            panel.border=element_rect(colour = "black", fill=NA, size=1))
    )
    
    # Remove score meta.data
    bm@meta.data[,n] <- NULL
}
dev.off()

# To look at individual markers:
FeaturePlot(bm, features = "FCGR3A", reduction = "umap", cols = colItay(c(1:100))) + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
# Regarding cluster 23, they are not pDCs, because these are some pDC markers Andy Lane suggested to look at on 220623:
as_tibble(bm@meta.data) %>% ggplot(aes(x = UMAP_1, y = UMAP_2, color = seurat_clusters == 23)) + geom_point() + theme(aspect.ratio = 1)
FeaturePlot(bm, features = "IRF8", reduction = "umap", cols = colItay(c(1:100))) + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) + coord_cartesian(xlim = c(0,14), ylim = c(5,13))
FeaturePlot(bm, features = "TCF4", reduction = "umap", cols = colItay(c(1:100))) + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) + coord_cartesian(xlim = c(0,14), ylim = c(5,13))
FeaturePlot(bm, features = "CLEC4C", reduction = "umap", cols = colItay(c(1:100))) + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) + coord_cartesian(xlim = c(0,14), ylim = c(5,13))
FeaturePlot(bm, features = "NRP1", reduction = "umap", cols = colItay(c(1:100))) + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) + coord_cartesian(xlim = c(0,14), ylim = c(5,13))
FeaturePlot(bm, features = "TLR7", reduction = "umap", cols = colItay(c(1:100))) + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) + coord_cartesian(xlim = c(0,14), ylim = c(5,13))
FeaturePlot(bm, features = "TLR9", reduction = "umap", cols = colItay(c(1:100))) + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) + coord_cartesian(xlim = c(0,14), ylim = c(5,13))
FeaturePlot(bm, features = "IL3RA", reduction = "umap", cols = colItay(c(1:100))) + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) + coord_cartesian(xlim = c(0,14), ylim = c(5,13))


# Name clusters from bone marrow data -------------------------------------------------------------

pdf("2_BM_ClusterNames.pdf", width = 6, height = 6)

# Plot Seurat clusters
DimPlot(bm, reduction = "umap", label = TRUE, pt.size = 0.5) +
  theme(aspect.ratio = 1, axis.line=element_blank(), axis.text.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title.x=element_blank(),
        axis.title.y=element_blank(), panel.border=element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(hjust = 0.5)) + ggtitle("Seurat clusters")

# Name clusters
bm <- RenameIdents(object = bm,
                   "0" = "T",
                   "1" = "Mono",
                   "2" = "T",
                   "3" = "T",
                   "4" = "EarlyEry",
                   "5" = "B",
                   "6" = "ProMono",
                   "7" = "LateEry",
                   "8" = "NK",
                   "9" = "HSPC",
                   "10" = "LateEry",
                   "11" = "cDC",
                   "12" = "ProB",
                   "13" = "NK",
                   "14" = "GMP",
                   "15" = "EarlyEry",
                   "16" = "Plasma",
                   "17" = "ncMono",
                   "18" = "pDC",
                   "19" = "PreB",
                   "21" = "ProB",
                   "22" = "T",
                   "23" = "cDC")

# Annotate doublets
bm@active.ident <- factor(ifelse(bm$Doublets, yes = "Doublets", no = as.character(bm@active.ident)),
                          levels = c("HSPC", "EarlyEry", "LateEry", "GMP", "ProMono", "Mono", "ncMono",
                                     "cDC", "pDC", "ProB", "PreB", "B", "Plasma", "T", "NK", "Doublets"))

# Add CellType metadata
bm$CellType <- bm@active.ident

# Plot named clusters
DimPlot(bm, reduction = "umap", label = TRUE, pt.size = 0.3, cols = my_colors[levels(bm$CellType)]) +
  theme(aspect.ratio = 1, axis.line=element_blank(), axis.text.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title.x=element_blank(),
        axis.title.y=element_blank(), panel.border=element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(hjust = 0.5)) + ggtitle("Named clusters")

dev.off()


# TNK cells: plot signature scores ----------------------------------------------------------------

# Average gene expression for scoreSignature function (to save time)
tnk_expr.mat.mean <- rowMeans(tnk_expr.mat)

# Load T/NK cell signatures from Yoke Seng Lee (unpublished), subset for gene names in current dataset. These signatures are not synced to Github to prevent confusion; use markerGenes_TNK.txt instead. Email Peter van Galen with questions.
tnk_signs.dt <- fread("markerGenes_Yoke_TNK.txt")
signs.ls <- lapply(tnk_signs.dt, intersect, rownames(tnk_expr.mat))
names(signs.ls) <- gsub(" |-", ".", names(signs.ls))

pdf("1_TNK_SignatureScores.pdf", width = 6, height = 6)

DimPlot(tnk) + theme(aspect.ratio = 1)

tnk <-  CellCycleScoring(tnk, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
FeaturePlot(tnk, features = "G2M.Score", cols = colItay(c(1:100))) + theme(aspect.ratio = 1)

# Plot signature scores
for (n in names(signs.ls)) {
  message(n)
  
  # Calculate signature score. This is quite similar to Seurat's function AddModuleScore.
  scores <- scoreSignature(CM = tnk_expr.mat, signatures = signs.ls[[n]], CM.mean = tnk_expr.mat.mean, verbose = TRUE)
  
  # Add as meta.data to Seurat object and plot
  tnk <- AddMetaData(tnk, metadata = scores, col.name = n)
  print(
    FeaturePlot(tnk, features = n, reduction = "umap", cols = colItay(c(1:100))) +
      theme(aspect.ratio = 1, axis.line=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(),
            axis.ticks=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(),
            panel.border=element_rect(colour = "black", fill=NA, size=1))
  )
  
  # Remove score meta.data
  tnk@meta.data[,n] <- NULL
}
dev.off()


# Name TNK clusters -------------------------------------------------------------------------------

pdf("2_TNK_ClusterNames.pdf", width = 6, height = 6)

# Plot Seurat clusters
DimPlot(tnk, reduction = "umap", label = TRUE, pt.size = 0.5) +
  theme(aspect.ratio = 1, axis.line=element_blank(), axis.text.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title.x=element_blank(),
        axis.title.y=element_blank(), panel.border=element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(hjust = 0.5)) + ggtitle("Seurat clusters")

# Name clusters
tnk <- RenameIdents(object = tnk,
                   "0" = "CD4Naive",
                   "1" = "CD4Memory",
                   "2" = "CD8TermExh",
                   "3" = "NK", # "",
                   "4" = "CD8Memory",
                   "5" = "CD8Naive",
                   "6" = "NKT",
                   "7" = "NK",
                   "8" = "GammaDeltaLike")

# Add CellType metadata
tnk$CellType <- factor(tnk@active.ident, levels = c("CD4Naive", "CD4Memory", "CD8Naive", "CD8Memory",
                                                    "CD8TermExh", "GammaDeltaLike", "NKT", "NK"))
tnk@active.ident <- tnk$CellType

# Plot named clusters
DimPlot(tnk, reduction = "umap", label = TRUE, pt.size = 0.3, cols = my_colors[names(my_colors) %in% tnk@active.ident]) +
  theme(aspect.ratio = 1, axis.line=element_blank(), axis.text.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title.x=element_blank(),
        axis.title.y=element_blank(), panel.border=element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(hjust = 0.5)) + ggtitle("Named clusters")

dev.off()


# Combine, plot and save --------------------------------------------------------------------------

# Merge cell type annotations ---------------------------------------------------------------------
stopifnot( all( colnames(tnk) %in% colnames(bm) ) )

# Make metadata tibbles
all_metadata_tib <- as_tibble(bm@meta.data, rownames = "cell")
tnk_metadata_tib <- as_tibble(tnk@meta.data, rownames = "cell")

# Join
merge_metadata_tib <- all_metadata_tib %>% left_join(select(tnk_metadata_tib, cell, CellType), by = "cell")
merge_metadata_tib <- merge_metadata_tib %>% mutate(CellType = ifelse(is.na(CellType.y), yes = as.character(CellType.x), no = as.character(CellType.y)))
merge_metadata_tib$CellType %>% table(useNA = "always")

# Overwrite with more granular clusters
bm$CellType <- NULL
bm@meta.data$CellType <- factor(merge_metadata_tib$CellType, levels = c("HSPC", "EarlyEry", "LateEry",
  "GMP", "ProMono", "Mono", "ncMono", "cDC", "pDC", "ProB", "PreB", "B", "Plasma",
  "CD4Naive", "CD4Memory", "CD8Naive", "CD8Memory", "CD8TermExh", "GammaDeltaLike", "NKT", "NK", "Doublets"))
bm@active.ident <- bm$CellType

# Plot
pdf("3_Combined_ClusterNames.pdf", width = 10, height = 6)
DimPlot(bm, reduction = "umap", label = TRUE, pt.size = 0.3, cols = my_colors[levels(bm@active.ident)]) +
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
BMfreq.mat <- BMfreq.mat[-grep("Doublets", rownames(BMfreq.mat)),]
BMnorm.mat <- sweep(BMfreq.mat, 2, colSums(BMfreq.mat), "/")*100

pdf("4_BM_CellTypes.pdf", width = 6, height = 6)
par(mar = c(8,4,8,12), xpd = T)

barplot(BMnorm.mat[nrow(BMnorm.mat):1,], col = rev(my_colors[rownames(BMnorm.mat)]),
        xaxt = "n", ylab = "Population frequency (%)", border = NA)
axis(side = 1, at = seq(1,ncol(BMnorm.mat))*1.2-0.5, labels = colnames(BMnorm.mat), las = 2)
legend(x = ncol(BMnorm.mat)*1.2+0.5, y = 100, legend = rownames(BMnorm.mat),
       fill = my_colors[rownames(BMnorm.mat)], bty = "n", border = NA)

dev.off()


# Marker genes ------------------------------------------------------------------------------------

# This is used as input for the random forest classifier later
markerGenes <- FindAllMarkers(bm, slot = "data", only.pos = T)
markergenes.dt.ls <- lapply(split(markerGenes, f = markerGenes$cluster), function(x) data.table(x))
markergenes.dt.ls <- lapply(markergenes.dt.ls, function(x) setorder(x, -avg_log2FC))
markergenes.df <- do.call(cbind, lapply(markergenes.dt.ls, function(x) x$gene[1:50]))

markergenesOrder.df <- markergenes.df[,levels(bm$CellType)]

write.table(markergenesOrder.df, file = "markerGenes.txt", sep = "\t", col.names = T, row.names = F, quote = F)

# Marker genes: TNK subset
markerGenes <- FindAllMarkers(tnk, slot = "data", only.pos = T)
markergenes.dt.ls <- lapply(split(markerGenes, f = markerGenes$cluster), function(x) data.table(x))
markergenes.dt.ls <- lapply(markergenes.dt.ls, function(x) setorder(x, -avg_log2FC))
markergenes.df <- do.call(cbind, lapply(markergenes.dt.ls, function(x) x$gene[1:50]))

markergenesOrder.df <- markergenes.df[,levels(tnk$CellType)]

write.table(markergenesOrder.df, file = "markerGenes_TNK.txt", sep = "\t", col.names = T, row.names = F, quote = F)


