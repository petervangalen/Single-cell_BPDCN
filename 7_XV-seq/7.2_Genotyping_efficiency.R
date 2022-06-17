# Peter van Galen, 220522
# Determine genotyping efficiency of every mutation

library(tidyverse)
library(Seurat)
library(readxl)
library(data.table)

rm(list=ls())

# Set working directory
setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/Github/7_XV-seq")

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")

# Load Seurat objects
seurat_files <- list.files("../7_XV-seq", pattern = "*.rds", full.names = T)
seu_bpdcn.ls <- lapply(seurat_files, function(x) readRDS(x))
names(seu_bpdcn.ls) <- gsub("_Seurat_Anno.rds", "", cutf(seurat_files, d = "/", f = 3))

# Generate data frame for heatmap annotations
genotyping_tables.tib <- read_excel("../7_XV-seq/FilteredCells_files.xlsx")
# Replace different MTAP primers with one, just as in 7_Add_GoT-XV_to_Seurat.R
genotyping_tables.tib$Mut <- gsub("MTAP.rearr.*", "MTAP.rearr", genotyping_tables.tib$Mut)
genotyping_tables.tib <- genotyping_tables.tib %>% select(Sample, Mut) %>% unique

# Integrate mutation data, plot heatmaps ----------------------------------------------------------

genotyping_tables.tib$Efficiency <- NA

for (n in 1:nrow(genotyping_tables.tib)) {
  #n <- 6
  Sample <- genotyping_tables.tib$Sample[n]
  mut <- genotyping_tables.tib$Mut[n]
  seu <- seu_bpdcn.ls[[Sample]]
  
  # Extract the necessary columns from Seurat metadata
  calls <- seu@meta.data[,mut]
  
  # Extract numbers
  genotyped <- sum(grepl("wildtype|mutant", calls))
  no_call <- sum(grepl("no call", calls))
  total <- length(calls)
  stopifnot(genotyped + no_call == total)
  
  genotyping_tables.tib$Efficiency[n] <- genotyped/total*100
  
}

write_tsv(genotyping_tables.tib, file = "genotyping_efficiency.txt")
