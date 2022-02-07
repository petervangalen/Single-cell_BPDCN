# Peter van Galen, 220205
# Predict BPDCN cell type identity based on normal BM

library(tidyverse)
library(Seurat)
library(readxl)
library(randomForest)
library(gplots)

setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/Github/3_RandomForest")

rm(list=ls())

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")
cell_colors <- popcol.tib$hex[1:17]
names(cell_colors) <- popcol.tib$pop[1:17]

# Marker genes for the 17 cell types
markerGenes.df <- read.table("../2_Annotate/markerGenes.txt", header = T)
selected.genes <- unique( unlist(markerGenes.df, use.names = F) )

# Load BM Data (Seurat object)
bm <- readRDS("../2_Annotate/BM_Seurat_CellTypes.rds")


# Load all gene expression data -------------------------------------------------------------------

# Make Seurat's count matrix easier to access (log, transcript per 10K)
bm.expr <- as.matrix( GetAssayData(object = bm, slot = "data") )

# Plot named clusters
DimPlot(bm, reduction = "umap", label = TRUE, pt.size = 0.5, cols = cell_colors) +
    theme(aspect.ratio = 1, axis.line=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(),
          axis.title.x=element_blank(), axis.title.y=element_blank(), panel.border=element_rect(colour = "black", fill=NA, size=1),
          plot.title = element_text(hjust = 0.5)) + ggtitle("Named clusters")

# Load and process BPDCN628 SW data similar to the healthy bone marrow controls (see 200413_Seurat_Harmony.R)
load("../../ExprStar/191111.190314.BPDCN628.both/191111.190314.BPDCN628.both.filter.RData", verbose = T)
# Check if you can add columns from CM.stats as features in Seurat
all(CM.stats$cell == colnames(CM.df))
# Subset for genes in bm object / exclude some odd-looking genes: setdiff(rownames(CM.df), rownames(bm))
CM.df <- CM.df[rownames(bm),]
# Harmonize cell ids
colnames(CM.df) <- paste0(cutf(colnames(CM.df), d = "_", f = 2), "-", gsub("BPDCN628-", "Pt9Dx.", cutf(colnames(CM.df), d = "_", f = 1)))
# Create Seurat object and add metadata
Pt9Dx.seu <- CreateSeuratObject(counts = CM.df, project = "Pt9Dx")
Pt9Dx.seu <- NormalizeData(Pt9Dx.seu) # (log, transcript per 10K)
Pt9Dx.seu$replicate <- cutf(colnames(Pt9Dx.seu), d = "-", f = 2)
Pt9Dx.seu$tech <- "SW"
Pt9Dx.seu$my.tSNE.x <- CM.stats$tSNEx
Pt9Dx.seu$my.tSNE.y <- CM.stats$tSNEy
# Create count matrix that is easier to work with for randomForest
Pt9Dx.cm <- as.matrix( GetAssayData(Pt9Dx.seu, slot = "data") )

# Function to load and process 10X data in the same way
Process_10X <- function(RData, patient_id) {
    load(RData, verbose = T)
    stopifnot(all(colnames(CM.dgm) == stats.dt$cell))
    CM.df <- as.matrix( CM.dgm[rownames(bm),] )
    colnames(CM.df) <- colnames(CM.df) %>% str_replace("-1", str_c("-", patient_id, ".1")) %>%
        str_replace("-2", str_c("-", patient_id, ".2")) %>% str_replace("-3", str_c("-", patient_id, ".3")) %>%
        str_replace("-4", str_c("-", patient_id, ".4"))
    seu <- CreateSeuratObject(counts = CM.df, project = patient_id)
    seu <- NormalizeData(seu) # (log, transcript per 10K)
    seu$replicate <- cutf(colnames(seu), d = "-", f = 2)
    seu$tech <- "TenX"
    seu$my.tSNE.x <- stats.dt$tSNEx
    seu$my.tSNE.y <- stats.dt$tSNEy
    return(seu)
}

# Create tibble with information on BPDCN samples that were analyzed with 10x
TenX_samples.tib <- tribble(~Patient_ID, ~RData, ~Note,
                            "Pt1Dx", "../../CellRanger/PT1-DX_agg/PT1-DX_agg.RData", "Added for revision",
                            "Pt1Mrd", "../../CellRanger/PT1-MRD_agg/PT1-MRD_agg.RData", "Added for revision",
                            "Pt5Dx", "../../CellRanger/PT5-DX_agg/PT5-DX_agg.RData", "Added for revision",
                            "Pt10Dx", "../../CellRanger/BPDCN712_agg/BPDCN712_agg.RData", "Included in submission",
                            "Pt10Rel", "../../CellRanger/BPDCN712R_agg/BPDCN712R_agg.RData", "Included in submission",
                            "Pt12Dx", "../../CellRanger/PT12-DX_agg/PT12-DX_agg.RData", "Added for revision",
                            "Pt12Rel", "../../CellRanger/PT12-REL_agg/PT12-REL_agg.RData", "Added for revision") #,
                            # "BPDCN181128", "../../CellRanger/BPDCN181128_agg/BPDCN181128_agg.RData", "Erica's project",
                            # "BPDCN180329", "../../CellRanger/BPDCN180329_agg/BPDCN180329_agg.RData", "Erica's project",
                            # "BPDCN190711", "../../CellRanger/BPDCN190711_agg/BPDCN190711_agg.RData", "Erica's project")

# Create list of processed Seurat objects
seu_10X.ls <- vector(mode = "list", length = nrow(TenX_samples.tib))
for (x in 1:nrow(TenX_samples.tib)) {
    seu_10X.ls[[x]] <- Process_10X(RData = TenX_samples.tib[[x,"RData"]],
                                   patient_id = TenX_samples.tib[[x,"Patient_ID"]])
}
names(seu_10X.ls) <- TenX_samples.tib$Patient_ID
seu.ls <- c(Pt9Dx = list(Pt9Dx.seu), seu_10X.ls)
# Change order
seu.ls <- seu.ls[c("Pt1Dx", "Pt1Mrd", "Pt5Dx", "Pt9Dx", "Pt10Dx", "Pt10Rel", "Pt12Dx", "Pt12Rel")]
expr.ls <- lapply(seu.ls, function(x) as.matrix(GetAssayData(x), slot = "data"))

# Merge all gene expression data & check if it all makes sense.
all.mat <- cbind(bm.expr, do.call(cbind, expr.ls))
table(cutf(colnames(all.mat), d = "-", f = 2))
table(round(colSums(exp(all.mat)-1)))
all(colnames(bm.expr) == names(bm@active.ident))


# Build classifier --------------------------------------------------------------------------------

# Plant classification trees
set.seed(123)
rf <- randomForest(x = t(bm.expr[selected.genes,]),    # matrix of predictors
                   y = bm@active.ident,              # response vector
                   sampsize = rep(50, length(levels(bm@active.ident))),
                   ntree = 1000,
                   do.trace = 100)

# Plot confusion matrix (based on out-of-bag data), with colors normalized for the total cell number in each population
Conf.mat <- rf$confusion[, rownames(rf$confusion)]
NormConf.mat <- Conf.mat / rowSums(Conf.mat)

pdf("1_ConfusionMatrix.pdf", width = 8, height = 8)
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

pdf("2_CrossValidation.pdf", width = 8, height = 8)
heatmap.2(NormConf.cv.mat, Rowv = F, Colv = F, dendrogram = "none", scale = "none", zlim = c(0, 1),
          col = colCustom(seq(0, 1, 0.01), color = c("white", "red")), trace = "none", density.info = "none",
          colsep = c(0, ncol(NormConf.cv.mat)), rowsep = c(0, nrow(NormConf.cv.mat)), sepcolor = "black",sepwidth = rep(0.01, 4),
          main = paste0("Confusion matrix, ", round(sum(diag(Conf.cv.mat)) / sum(Conf.cv.mat)*100, 2), "% accurate"),
          add.expr = text(rep(1:ncol(NormConf.cv.mat), each=nrow(NormConf.cv.mat)),
                          rep(ncol(NormConf.cv.mat):1, nrow(NormConf.cv.mat)), Conf.cv.mat))
dev.off()


# Classify cells from patients --------------------------------------------------------------------

# Predict cell types in each of the BPDCN samples (ties are broken at random)
predictions.mat.ls <- lapply(c(BM = list(bm.expr), expr.ls), function(x) predict(rf, t(x[selected.genes,]), type = "prob"))

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
PlotFreq.mat <- cbind(BMfreq.mat, PredictFreq.mat[,-1])[-match("Doublets", rownames(PredictFreq.mat)),]
# Normalize to 100
PlotFreqNorm.mat <- sweep(PlotFreq.mat, 2, colSums(PlotFreq.mat), "/")*100

# Test if there are significant differences in cell type frequencies between BM vs. Dx
SubsetFreqNorm.mat <- PlotFreqNorm.mat[,grepl("BM|Dx", colnames(PlotFreqNorm.mat))]
apply(SubsetFreqNorm.mat, 1, function(z) t.test(x = z[1:3], y = z[c(4:8)])$p.value)

pdf("3_CellTypeFrequencies.pdf", width = 8, height = 6)
par(mar = c(8,4,8,12), xpd = T)

barplot(PlotFreqNorm.mat[nrow(PlotFreqNorm.mat):1,], col = rev(cell_colors[-17]), xaxt = "n", ylab = "Population frequency (%)", border = NA)
axis(side = 1, at = seq(1,ncol(PlotFreqNorm.mat))*1.2-0.5, labels = colnames(PlotFreqNorm.mat), las = 2)
legend(x = ncol(PlotFreqNorm.mat)*1.2+0.5, y = 100, legend = rownames(PlotFreqNorm.mat), fill = cell_colors[-17], bty = "n", border = NA)

dev.off()


# Project cell types & save -----------------------------------------------------------------------

# Coordinates of normal cells (change UMAP_1 and UMAP_2 to facilitate workflow)
bm.umap <- data.frame(bm@reductions$umap@cell.embeddings)

# Project and save the BPDCN samples
for ( Patient_ID in names(expr.ls) ) {
#Patient_ID <- names(predictions.mat.ls)[2]
message(Patient_ID)

# Correlate BPDCN cell prediction scores to BM prediction scores
cor.mat <- cor(t(predictions.mat.ls[["BM"]]), t(predictions.mat.ls[[Patient_ID]]))
cor.mat[1:3,1:3]
cor.max <- apply(cor.mat, 2, which.max)
cor.id <- rownames(predictions.mat.ls[["BM"]])[cor.max]

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
s$my.tSNE.x <- NULL
s$my.tSNE.y <- NULL
s <- subset(s, subset = CellType != "Doublets")
saveRDS(s, file = paste0(Patient_ID, "_Seurat_Predict.rds"))

}


# Save normal BM object ---------------------------------------------------------------------------

# Add prediction scores to bm & save
#for (x in colnames(predictions.mat.ls[["BM"]])) {
#    bm <- AddMetaData(bm, predictions.mat.ls[["BM"]][,x], col.name = paste0("Predict.", x))
#}

# Some additional wrangling
#bm <- subset(bm, subset = CellType != "Doublets")
#bm <- DietSeurat(bm)

#saveRDS(bm, file = "BM_Seurat_Predict.rds")




