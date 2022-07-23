# Peter van Galen, 220207
# Generate allowlist of high-quality cell barcodes for BPDCN samples to use for GoT

library(tidyverse)
library(Seurat)

setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/Github/5_Stats")

rm(list=ls())

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")

# Load list of Seurat files
seurat_files <- list.files("../4_XV-seq", pattern = "*.rds", full.names = T)
seu.ls <- lapply(seurat_files, function(x) readRDS(x))
names(seu.ls) <- cutf(basename(seurat_files), d = "_")
seu_bpdcn.ls <- seu.ls[-1]

for (n in names(seu_bpdcn.ls)) {
  seu <- seu_bpdcn.ls[[n]]
  allowlist.df <- data.frame(CellBarcode = cutf(colnames(seu), d = "-"))
  write.table(allowlist.df, file = paste0(n, "_CellBarcode_Allowlist.txt"), row.names = F, col.names = F, quote = F)
}

# For bone marrow, do the same, and also make a dataframe with 100 randomly chosen good cells (for another project)
#bm <- seu.ls[["BM"]]

# Also make a file with 100 random cells
#good_cells.df <- data.frame(matrix(ncol = length(unique(bm$replicate)), nrow = 100))
#colnames(good_cells.df) <- unique(bm$replicate)

#for (x in unique(bm$replicate)) {
  #bm_subset <- subset(bm, replicate == x)
  #allowlist.df <- data.frame(CellBarcode = cutf(colnames(bm_subset), d = "-"))
  #print(x); print(nrow(allowlist.df))
  #write.table(allowlist.df, file = paste0(x, "_CellBarcode_Allowlist.txt"), row.names = F, col.names = F, quote = F)
  #good_cells <- paste0(sample(allowlist.df$CellBarcode, size = 100), "-1")
  #good_cells.df[,x] <- good_cells
  #write.table(file = paste0(x, "_good_cells.txt"), data.frame(good_cells), row.names = F, col.names = F, quote = F, sep = "\t")
#}

#write.table(file = "good_cells.txt", good_cells.df, row.names = F, quote = F, sep = "\t")