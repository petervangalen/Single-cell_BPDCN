# Peter van Galen, 220703
# Number of cells, UMIs, etc

library(tidyverse)
library(Seurat)

setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/05_Stats")

rm(list=ls())

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")

# Load Seurat objects
seurat_files <- list.files("../04_XV-seq", pattern = "*.rds", full.names = T)
seu.ls <- lapply(seurat_files, function(x) readRDS(x))
names(seu.ls) <- cutf(basename(seurat_files), d = "_")

# Split bm, merge all
seu.ls$BM$orig.ident2 <- cutf(seu.ls$BM$replicate, d = "\\.")
seu_all_ls <- c(SplitObject(seu.ls$BM, split.by = "orig.ident2"), seu.ls[-1])

# Create stats table
stats.df <- data.frame(Sample_ID = names(seu_all_ls),
                       Cell_number = unlist(lapply(seu_all_ls, ncol)),
                       UMIs_mean = unlist(lapply(seu_all_ls, function(x) mean(x$nCount_RNA))),
                       UMIs_min = unlist(lapply(seu_all_ls, function(x) min(x$nCount_RNA))),
                       UMIs_max = unlist(lapply(seu_all_ls, function(x) max(x$nCount_RNA))),
                       Genes_mean = unlist(lapply(seu_all_ls, function(x) mean(x$nFeature_RNA))),
                       Genes_min = unlist(lapply(seu_all_ls, function(x) min(x$nFeature_RNA))),
                       Genes_max = unlist(lapply(seu_all_ls, function(x) max(x$nFeature_RNA))))
stats_order.df <- stats.df[c("BM1", "BM2", "BM3", "BM4", "BM5", "BM6", "Pt1Dx", "Pt1Rem", "Pt5Dx", "Pt9Dx",
                             "Pt10Dx", "Pt10Rel", "Pt12Dx", "Pt12Rel", "Pt14Dx", "Pt15Dx", "Pt16Dx"),]
write.table(stats_order.df, file = "stats.txt", row.names = F, sep = "\t", quote = F)


