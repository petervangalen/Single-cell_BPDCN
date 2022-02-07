# Peter van Galen, 220207
# Generate allowlist of high-quality cell barcodes for BPDCN samples to use for GoT

library(tidyverse)
library(Seurat)

setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/Github/6_Barcode_AllowLists")

rm(list=ls())

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")

# Load list of Seurat files
seurat_files <- list.files("../3_RandomForest", pattern = "*.rds", full.names = T)
bpdcn_seu.ls <- lapply(seurat_files, function(x) readRDS(x))
names(bpdcn_seu.ls) <- cutf(basename(seurat_files), d = "_")


for (n in names(bpdcn_seu.ls)) {
  seu <- bpdcn_seu.ls[[n]]
  allowlist.df <- data.frame(CellBarcode = cutf(colnames(seu), d = "-"))
  write.table(allowlist.df, file = paste0(n, "_CellBarcode_Allowlist.txt"), row.names = F, col.names = F, quote = F)
}