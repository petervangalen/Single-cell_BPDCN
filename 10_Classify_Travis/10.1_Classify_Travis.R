# Peter van Galen, 220727
# Use the random forest classifier to evaluate if there are any pDCs in Travis' skin scRNA-seq data (https://doi.org/10.1016/j.immuni.2020.09.015)
# This is similar to 3_RandomForest.R

library(tidyverse)
library(Seurat)
library(readxl)
library(randomForest)
library(gplots)

setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/Github/10_Classify_Travis")

rm(list=ls())

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")
cell_colors <- popcol.tib$hex[1:23]
names(cell_colors) <- popcol.tib$pop[1:23]

# Marker genes for the 23 cell types
markerGenes.df <- read.table("../2_Annotate/markerGenes.txt", header = T)
selected.genes <- unique( unlist(markerGenes.df, use.names = F) )

# Load BM Data (Seurat object)
bm <- readRDS("../2_Annotate/BM_Seurat_CellTypes.rds")
analysis_name <- ""
# OR
bm <- readRDS("../4_XV-seq/BM_Seurat_Final.rds")
analysis_name <- "_NoDoublets"

# Load single-cell count matrix from Travis' paper (prepared by Volker)
load("~/DropboxMGB/Projects/Single-cell_BPDCN/Travis_singlecell_skin/GSE150672_Skin_Expression_counts.Normal_4840cells.RData")

# Get started. The cell annotations are available but not included here.
seu <- CreateSeuratObject(counts = d.skin)
seu <- NormalizeData(seu)
seu@meta.data %>% head

# Subset everything for common genes (b/c Travis and I used different reference genomes)
common.genes <- intersect(rownames(GetAssayData(bm)), rownames(GetAssayData(seu)))
selected.genes.subset <- intersect(selected.genes, common.genes)
bm <- subset(bm, features = common.genes)
seu <- subset(seu, features = common.genes)

# Plot named clusters
DimPlot(bm, reduction = "umap", label = TRUE, pt.size = 0.5, cols = cell_colors) +
  theme(aspect.ratio = 1, axis.line=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(),
        axis.title.x=element_blank(), axis.title.y=element_blank(), panel.border=element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(hjust = 0.5)) + ggtitle("Named clusters")


# Build classifier --------------------------------------------------------------------------------

# Make Seurat's count matrix easier to access (log, transcript per 10K)
bm.expr <- as.matrix( GetAssayData(object = bm, slot = "data") )

# Plant classification trees
set.seed(123)
rf <- randomForest(x = t(bm.expr[selected.genes.subset,]),  # matrix of predictors
                   y = bm@active.ident,              # response vector
                   sampsize = rep(50, length(levels(bm@active.ident))),
                   ntree = 1000,
                   do.trace = 100)

# Plot confusion matrix (based on out-of-bag data), with colors normalized for the total cell number in each population
Conf.mat <- rf$confusion[, rownames(rf$confusion)]
NormConf.mat <- Conf.mat / rowSums(Conf.mat)

pdf(paste0("1_ConfusionMatrix", analysis_name, ".pdf"), width = 12, height = 12)
heatmap.2(NormConf.mat, Rowv = F, Colv = F, dendrogram = "none", scale = "none", zlim = c(0, 1),
          col = colCustom(seq(0, 1, 0.01), color = c("white", "red")), trace = "none", density.info = "none",
          colsep = c(0, ncol(NormConf.mat)), rowsep = c(0, nrow(NormConf.mat)), sepcolor = "black",sepwidth = rep(0.01, 4),
          main = paste0("Confusion matrix, ", round(sum(diag(Conf.mat)) / sum(Conf.mat)*100, 2), "% accurate"),
          add.expr = text(rep(1:ncol(NormConf.mat), each=nrow(NormConf.mat)),
                          rep(ncol(NormConf.mat):1, nrow(NormConf.mat)), Conf.mat))
dev.off()
# For the five-fold cross-validation, see 3_RandomForest.R


# Classify cells from Travis' paper ---------------------------------------------------------------

# Predict cell identities
predictions.mat <- predict(rf, t(as.matrix(GetAssayData(seu, slot = "data"))[selected.genes.subset,]), type = "prob")

# Make the prediction matrix into a factor of cell types
classify.fc <- factor(colnames(predictions.mat)[apply(predictions.mat, 1, which.max)], levels = colnames(predictions.mat))
names(classify.fc) <- rownames(predictions.mat)
table( classify.fc )

as_tibble(predictions.mat, rownames = "cell") %>% write_tsv(file = paste0("Prediction_scores", analysis_name, ".txt"))
as_tibble(data.frame(MaxPrediction = classify.fc), rownames = "cell") %>% write_tsv(file = paste0("Predictions", analysis_name, ".txt"))



# >>>>>>


# <<<<<<



# Plot the clustered cell type frequencies in BM (per donor) and the predicted cell type frequencies in BPDCN
BMfreq.mat <- do.call(cbind, lapply(split(bm@meta.data, f = cutf(bm@meta.data$replicate, d = "\\.")), function(x) table(x$CellType)))
PredictFreq.mat <- do.call(cbind, lapply(CellTypes.ls, table))
# Percent doublets. For healthy donor doublets, see bm$CellType
apply(PredictFreq.mat, 2, function(x) x["Doublets"] / sum(x) * 100)
PlotFreq.mat <- cbind(BMfreq.mat, PredictFreq.mat)
PlotFreq.mat <- PlotFreq.mat[-match("Doublets", rownames(PlotFreq.mat)),]
# Normalize to 100
PlotFreqNorm.mat <- sweep(PlotFreq.mat, 2, colSums(PlotFreq.mat), "/")*100

# Test if there are significant differences in cell type frequencies between BM vs. Dx --> No there aren't.
SubsetFreqNorm.mat <- PlotFreqNorm.mat[,grepl("BM|Dx", colnames(PlotFreqNorm.mat))]
apply(SubsetFreqNorm.mat, 1, function(z) t.test(x = z[1:3], y = z[c(4:8)])$p.value)

pdf("3_CellTypeFrequencies.pdf", width = 10, height = 6)
par(mar = c(8,4,8,12), xpd = T)

barplot(PlotFreqNorm.mat[nrow(PlotFreqNorm.mat):1,], col = rev(cell_colors[rownames(PlotFreqNorm.mat)]),
        xaxt = "n", ylab = "Population frequency (%)", border = NA)
axis(side = 1, at = seq(1,ncol(PlotFreqNorm.mat))*1.2-0.5, labels = colnames(PlotFreqNorm.mat), las = 2)
legend(x = ncol(PlotFreqNorm.mat)*1.2+0.5, y = 100, legend = rownames(PlotFreqNorm.mat), fill = cell_colors[rownames(PlotFreqNorm.mat)], bty = "n", border = NA)

dev.off()


# Project cell types & save -----------------------------------------------------------------------

# Coordinates of normal cells (change UMAP_1 and UMAP_2 to facilitate workflow)
bm.umap <- data.frame(bm@reductions$umap@cell.embeddings)

# Project and save the BPDCN samples
for ( Patient_ID in names(expr.ls) ) {
  #Patient_ID <- names(predictions.mat.ls)[2]
  message(Patient_ID)
  
  # Correlate BPDCN cell prediction scores to BM prediction scores
  cor.mat <- cor(t(bm.predictions.mat), t(predictions.mat.ls[[Patient_ID]]))
  cor.mat[1:3,1:3]
  cor.max <- apply(cor.mat, 2, which.max)
  cor.id <- rownames(bm.predictions.mat)[cor.max]
  
  # Plot single nearest cells
  bpdcn.project.umap <- data.frame(row.names = colnames(cor.mat),
                                   project.umap.x = bm.umap[cor.id, "UMAP_1"],
                                   project.umap.y = bm.umap[cor.id, "UMAP_2"],
                                   CellType = CellTypes.ls[[Patient_ID]],
                                   CellTypeCol = cell_colors[as.character(CellTypes.ls[[Patient_ID]])])
  
  pdf(paste0(Patient_ID, "_cor.predict.pdf"), width = 6, height = 6)
  par(mar=c(4, 4, 4, 4))
  
  plotTSNE(bm.umap, cex = 0.3, col = "#DDDDDD", main = paste(Patient_ID, "on normal BM (grey)"))
  points(bpdcn.project.umap[,1:2], pch = 16, cex = 0.3, col = bpdcn.project.umap$CellTypeCol)
  
  # How does it look for the 2, 5, 10 and 20 most similar cells?
  for (n in c(2, 5, 10, 20)) {
    #n <- 2
    plotTSNE(bm.umap, cex = 0.3, col = "#DDDDDD", main = paste(n, "nearest cells"))
    
    cor.n <- t(apply(cor.mat, 2, function(x) rownames(cor.mat)[ order(x, decreasing = T)[1:n] ] ))
    
    bpdcn.project.umap.n <- data.frame(row.names = rownames(cor.n),
                                       project.umap.x = apply(cor.n, 1, function(x) mean(bm.umap[x,"UMAP_1"])),
                                       project.umap.y = apply(cor.n, 1, function(x) mean(bm.umap[x,"UMAP_2"])),
                                       CellType = CellTypes.ls[[Patient_ID]],
                                       CellTypeCol = cell_colors[as.character(CellTypes.ls[[Patient_ID]])])
    
    points(bpdcn.project.umap.n[,1:2], pch = 16, cex = 0.3, col = bpdcn.project.umap$CellTypeCol)
  }
  
  dev.off()
  
  # Save Seurat object will all relevant data
  s <- seu.ls[[Patient_ID]]
  
  # Check
  stopifnot(all(rownames(bpdcn.project.umap) == colnames(s)))
  stopifnot(all(rownames(predictions.mat.ls[[Patient_ID]]) == colnames(s)))
  
  # Add prediction scores for each cell type
  #for (x in colnames(predictions.mat.ls[[Patient_ID]])) {
  #    s <- AddMetaData(s, predictions.mat.ls[[Patient_ID]][,x], col.name = paste0("Predict.", x))
  #}
  
  # Add other metadata, exclude doublets, save
  s$CellType <- bpdcn.project.umap$CellType
  s@active.ident <- s$CellType
  s$project.umap.x <- bpdcn.project.umap$project.umap.x
  s$project.umap.y <- bpdcn.project.umap$project.umap.y
  s@commands <- list()
  
  # Specify additional doublets based on their UMAP coordinates (same as 1_Seurat_Harmony.R) and exclude them with classified doublets
  more_doublets <- filter(as_tibble(s@meta.data, rownames = "cell"), between(project.umap.x, -7, -2), between(project.umap.y, 0, 3))$cell
  s <- subset(s, cells = setdiff(colnames(s), more_doublets))
  s <- subset(s, subset = CellType != "Doublets")
  
  saveRDS(s, file = paste0(Patient_ID, "_Seurat_Predict.rds"))
  
}




