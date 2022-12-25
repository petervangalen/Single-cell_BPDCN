# Peter van Galen, 220205
# Predict BPDCN cell type identity based on normal BM

library(tidyverse)
library(Seurat)
library(readxl)
library(randomForest)
library(gplots)

setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/03_RandomForest")

rm(list=ls())

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")
cell_colors <- popcol.tib$hex[1:22]
names(cell_colors) <- popcol.tib$pop[1:22]

# Marker genes for the 23 cell types
markerGenes.df <- read.table("../02_Annotate/markerGenes.txt", header = T)
selected.genes <- unique( unlist(markerGenes.df, use.names = F) )

# Load BM Data (Seurat object)
bm <- readRDS("../02_Annotate/BM_Seurat_CellTypes.rds")


# Load all gene expression data -------------------------------------------------------------------

# Make Seurat's count matrix easier to access (log, transcript per 10K)
bm.expr <- as.matrix( GetAssayData(object = bm, slot = "data") )

# Plot named clusters
DimPlot(bm, reduction = "umap", label = TRUE, pt.size = 0.5, cols = cell_colors) +
    theme(aspect.ratio = 1, axis.line=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(),
          axis.title.x=element_blank(), axis.title.y=element_blank(), panel.border=element_rect(colour = "black", fill=NA, size=1),
          plot.title = element_text(hjust = 0.5)) + ggtitle("Named clusters")

# Load and process BPDCN628 SW data similar to the healthy bone marrow controls (see 1_Seurat_Harmony.R)
load("../../ExprStar/191111.190314.BPDCN180628-1.star/191111.190314.BPDCN180628-1.filter.RData", verbose = T)
CM1.df <- CM.df
load("../../ExprStar/191111.190314.BPDCN180628-2.star/191111.190314.BPDCN180628-2.filter.RData", verbose = T)
CM.df <- cbind(CM1.df, CM.df)
# Subset for genes in bm object / exclude some odd-looking genes: setdiff(rownames(CM.df), rownames(bm))
CM.df <- CM.df[rownames(bm),]
# Harmonize cell ids
colnames(CM.df) <- paste0(cutf(colnames(CM.df), d = "_", f = 2), "-", gsub("BPDCN628-", "Pt9Dx.", cutf(colnames(CM.df), d = "_", f = 1)))
# Create Seurat object and add metadata
Pt9Dx.seu <- CreateSeuratObject(counts = CM.df, project = "Pt9Dx")
Pt9Dx.seu <- NormalizeData(Pt9Dx.seu) # log, transcript per 10K
Pt9Dx.seu$replicate <- cutf(colnames(Pt9Dx.seu), d = "-", f = 2)
Pt9Dx.seu$tech <- "SW"

# Function to load and process 10X data in the same way
Process_10X <- function(RData, patient_id) {
    load(RData, verbose = T)
    stopifnot(all(colnames(CM.dgm) == stats.dt$cell))
    CM.df <- as.matrix( CM.dgm[rownames(bm),] )
    colnames(CM.df) <- gsub(1, gsub("-", ".", patient_id), colnames(CM.dgm))
    seu <- CreateSeuratObject(counts = CM.df, project = patient_id)
    seu <- NormalizeData(seu) # log, transcript per 10K
    seu$replicate <- cutf(colnames(seu), d = "-", f = 2)
    seu$tech <- "TenX"
    return(seu)
}

# Create tibble with information on BPDCN samples that were analyzed with 10x
TenX_samples.tib <- tribble(~Patient_ID, ~RData, ~Note,
                            "Pt1Dx-2", "../../CellRanger/Pt1Dx-2/Pt1Dx-2.RData", "",
                            "Pt1Dx-3", "../../CellRanger/Pt1Dx-3/Pt1Dx-3.RData", "",
                            "Pt1Rem-1", "../../CellRanger/Pt1Rem-1/Pt1Rem-1.RData", "",
                            "Pt1Rem-2", "../../CellRanger/Pt1Rem-2/Pt1Rem-2.RData", "",
                            "Pt1Rem-3", "../../CellRanger/Pt1Rem-3/Pt1Rem-3.RData", "",
                            "Pt5Dx-1", "../../CellRanger/Pt5Dx-1/Pt5Dx-1.RData", "",
                            "Pt5Dx-3", "../../CellRanger/Pt5Dx-3/Pt5Dx-3.RData", "",
                            "Pt10Dx-1", "../../CellRanger/Pt10Dx-1/Pt10Dx-1.RData", "Included in submission",
                            "Pt10Dx-2", "../../CellRanger/Pt10Dx-2/Pt10Dx-2.RData", "Included in submission",
                            "Pt10Dx-3", "../../CellRanger/Pt10Dx-3/Pt10Dx-3.RData", "Included in submission",
                            "Pt10Dx-4", "../../CellRanger/Pt10Dx-4/Pt10Dx-4.RData", "Included in submission",
                            "Pt10Rel-1", "../../CellRanger/Pt10Rel-1/Pt10Rel-1.RData", "Included in submission",
                            "Pt10Rel-2", "../../CellRanger/Pt10Rel-2/Pt10Rel-2.RData", "Included in submission",
                            "Pt12Dx-1", "../../CellRanger/Pt12Dx-1/Pt12Dx-1.RData", "",
                            "Pt12Dx-2", "../../CellRanger/Pt12Dx-2/Pt12Dx-2.RData", "",
                            "Pt12Rel-1", "../../CellRanger/Pt12Rel-1/Pt12Rel-1.RData", "",
                            "Pt12Rel-2", "../../CellRanger/Pt12Rel-2/Pt12Rel-2.RData", "",
                            "Pt14Dx-1", "../../CellRanger/Pt14Dx-1/Pt14Dx-1.RData", "",
                            "Pt14Dx-2", "../../CellRanger/Pt14Dx-2/Pt14Dx-2.RData", "",
                            "Pt14Dx-3", "../../CellRanger/Pt14Dx-3/Pt14Dx-3.RData", "",
                            "Pt15Dx-5", "../../CellRanger/Pt15Dx-5/Pt15Dx-5.RData", "",
                            "Pt15Dx-6", "../../CellRanger/Pt15Dx-6/Pt15Dx-6.RData", "",
                            "Pt15Dx-7", "../../CellRanger/Pt15Dx-7/Pt15Dx-7.RData", "",
                            "Pt16Dx-1", "../../CellRanger/Pt16Dx-1/Pt16Dx-1.RData", "",
                            "Pt16Dx-2", "../../CellRanger/Pt16Dx-2/Pt16Dx-2.RData", "",
                            "Pt16Dx-3", "../../CellRanger/Pt16Dx-3/Pt16Dx-3.RData", "")

# Create list of processed Seurat objects
seu_10x.ls <- vector(mode = "list", length = nrow(TenX_samples.tib))
for (x in 1:nrow(TenX_samples.tib)) {
    seu_10x.ls[[x]] <- Process_10X(RData = TenX_samples.tib[[x,"RData"]],
                                   patient_id = TenX_samples.tib[[x,"Patient_ID"]])
}
names(seu_10x.ls) <- TenX_samples.tib$Patient_ID
seu.ls <- list(Pt1Dx = merge(seu_10x.ls$`Pt1Dx-2`, list(seu_10x.ls$`Pt1Dx-3`)),
               Pt1Rem = merge(seu_10x.ls$`Pt1Rem-1`, list(seu_10x.ls$`Pt1Rem-2`, seu_10x.ls$`Pt1Rem-3`)),
               Pt5Dx = merge(seu_10x.ls$`Pt5Dx-1`, list(seu_10x.ls$`Pt5Dx-3`)),
               Pt9Dx = Pt9Dx.seu,
               Pt10Dx = merge(seu_10x.ls$`Pt10Dx-1`, list(seu_10x.ls$`Pt10Dx-2`, seu_10x.ls$`Pt10Dx-3`, seu_10x.ls$`Pt10Dx-4`)),
               Pt10Rel = merge(seu_10x.ls$`Pt10Rel-1`, list(seu_10x.ls$`Pt10Rel-2`)),
               Pt12Dx = merge(seu_10x.ls$`Pt12Dx-1`, list(seu_10x.ls$`Pt12Dx-2`)),
               Pt12Rel = merge(seu_10x.ls$`Pt12Rel-1`, list(seu_10x.ls$`Pt12Rel-2`)),
               Pt14Dx = merge(seu_10x.ls$`Pt14Dx-1`, list(seu_10x.ls$`Pt14Dx-2`, seu_10x.ls$`Pt14Dx-3`)),
               Pt15Dx = merge(seu_10x.ls$`Pt15Dx-5`, list(seu_10x.ls$`Pt15Dx-6`, seu_10x.ls$`Pt15Dx-7`)),
               Pt16Dx = merge(seu_10x.ls$`Pt16Dx-1`, list(seu_10x.ls$`Pt16Dx-2`, seu_10x.ls$`Pt16Dx-3`)))

# Apply the same quality metrics to BPDCN as to healthy controls (see 1_Seurat_Harmony.R). In preprocessing, cells with >20% mitochondrial genes were already removed.
seu.ls <- lapply(seu.ls, subset, nCount_RNA > 2000 & nFeature_RNA > 1000)

# Remove replicate from orig.ident metadata (this is still saved in the replicate column)
for (n in names(seu.ls)) {
  seu.ls[[n]]$orig.ident <- cutf(as.character(seu.ls[[n]]$orig.ident), d = "-")
}

# Get expression matrices
expr.ls <- lapply(seu.ls, function(x) as.matrix(GetAssayData(x), slot = "data"))

# Merge all gene expression data & check if it all makes sense.
all.mat <- cbind(bm.expr, do.call(cbind, expr.ls))
data.frame(table(cutf(colnames(all.mat), d = "-", f = 2)))
table(round(colSums(exp(all.mat)-1)))
all(colnames(bm.expr) == names(bm@active.ident))


# Build classifier --------------------------------------------------------------------------------

# Plant classification trees
set.seed(123)
rf <- randomForest(x = t(bm.expr[selected.genes,]),  # matrix of predictors
                   y = bm@active.ident,              # response vector
                   sampsize = rep(50, length(levels(bm@active.ident))),
                   ntree = 1000,
                   do.trace = 100)

# Plot confusion matrix (based on out-of-bag data), with colors normalized for the total cell number in each population
Conf.mat <- rf$confusion[, rownames(rf$confusion)]
NormConf.mat <- Conf.mat / rowSums(Conf.mat)

pdf("1_ConfusionMatrix.pdf", width = 12, height = 12)
heatmap.2(NormConf.mat, Rowv = F, Colv = F, dendrogram = "none", scale = "none", zlim = c(0, 1),
          col = colCustom(seq(0, 1, 0.01), color = c("white", "red")), trace = "none", density.info = "none",
          colsep = c(0, ncol(NormConf.mat)), rowsep = c(0, nrow(NormConf.mat)), sepcolor = "black",sepwidth = rep(0.01, 4),
          main = paste0("Confusion matrix, ", round(sum(diag(Conf.mat)) / sum(Conf.mat)*100, 2), "% accurate"),
          add.expr = text(rep(1:ncol(NormConf.mat), each=nrow(NormConf.mat)),
                          rep(ncol(NormConf.mat):1, nrow(NormConf.mat)), Conf.mat))
dev.off()


# Five-fold cross-validation ----------------------------------------------------------------------

# Split dataset in five parts
cv <- split(colnames(bm.expr), rep(1:5, 1E6)[1:ncol(bm.expr)])

# Build five forests, each with 4/5 of the data
rf.cv <- lapply(cv, function(n) {
        #n <- cv[[1]]
        set.seed(123)
        randomForest(x = t(bm.expr[selected.genes,setdiff(colnames(bm.expr), n)]),
                     y = bm@active.ident[! colnames(bm.expr) %in% n],
                     sampsize = sapply(table(bm@active.ident[! colnames(bm.expr) %in% n]), min, 50), # check that it's not too low
                     ntree = 1000,
                     do.trace = 100)
})

# Predict the sets that were not used for training
rf.cv.predict.prob <- lapply(rf.cv, function(rf) {
        predict(rf, t(bm.expr[selected.genes, setdiff(colnames(bm.expr), names(rf$y))]), type = "prob")
})

# Maximum prpbability
rf.cv.predict <- lapply(rf.cv.predict.prob, function(x) {
        y <- factor(colnames(x)[apply(x, 1, which.max)], colnames(x))  # maximum sore
        names(y) <- rownames(x)
        y
})
head(rf.cv.predict[[1]])

# Confusion matrix
Conf.cv.mat <- table(bm@active.ident[unlist(lapply(rf.cv.predict, names))], unlist(rf.cv.predict))
NormConf.cv.mat <- Conf.cv.mat / rowSums(Conf.cv.mat)

pdf("2_CrossValidation.pdf", width = 12, height = 12)
heatmap.2(NormConf.cv.mat, Rowv = F, Colv = F, dendrogram = "none", scale = "none", zlim = c(0, 1),
          col = colCustom(seq(0, 1, 0.01), color = c("white", "red")), trace = "none", density.info = "none",
          colsep = c(0, ncol(NormConf.cv.mat)), rowsep = c(0, nrow(NormConf.cv.mat)), sepcolor = "black",sepwidth = rep(0.01, 4),
          main = paste0("Confusion matrix, ", round(sum(diag(Conf.cv.mat)) / sum(Conf.cv.mat)*100, 2), "% accurate"),
          add.expr = text(rep(1:ncol(NormConf.cv.mat), each=nrow(NormConf.cv.mat)),
                          rep(ncol(NormConf.cv.mat):1, nrow(NormConf.cv.mat)), Conf.cv.mat))
dev.off()


# Classify cells from patients --------------------------------------------------------------------

# Predict cell types in each of the BPDCN samples (ties are broken at random)
predictions.mat.ls <- lapply(expr.ls, function(x) predict(rf, t(x[selected.genes,]), type = "prob"))
# Perform similar analysis on bone marrow (for the next section)
bm.predictions.mat <- predict(rf, t(bm.expr[selected.genes,]), type = "prob")

# Maximum probability
CellTypes.ls <- lapply(predictions.mat.ls, function(x) {
        y <- factor(colnames(x)[apply(x, 1, which.max)], colnames(x))
        names(y) <- rownames(x)
        y
})

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

# For 11.2_Sina_or_barplots.R, save pDC prediction scores. These include doublets
saveRDS(predictions.mat.ls, "Prediction_scores.rds")


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




