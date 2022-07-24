# Peter van Galen, 220522
# Determine genotyping efficiency of every mutation

library(tidyverse)
library(Seurat)
library(readxl)
library(data.table)

rm(list=ls())

# Set working directory
setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/Github/4_XV-seq")

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")

# Load Seurat objects
seurat_files <- list.files("../4_XV-seq", pattern = "*.rds", full.names = T)
seu_bpdcn.ls <- lapply(seurat_files, function(x) readRDS(x))
names(seu_bpdcn.ls) <- gsub("_Seurat_Final.rds", "", cutf(seurat_files, d = "/", f = 3))

# Generate data frame for genotyping results
genotyping_tables.tib <- read_excel("../4_XV-seq/XV-seq_overview.xlsx")
# Replace different MTAP entries with one, just as in 4.1_Add_GoT-XV_to_Seurat.R
genotyping_tables.tib$Mutation <- gsub("MTAP.rearr.*", "MTAP.rearr", genotyping_tables.tib$Mutation)
genotyping_tables.tib <- genotyping_tables.tib %>% select(Sample, Mutation) %>% unique

# Integrate mutation data, quantify genotyping efficiency -----------------------------------------

genotyping_tables.tib$Efficiency <- NA

for (n in 1:nrow(genotyping_tables.tib)) {
  #n <- 6
  Sample <- genotyping_tables.tib$Sample[n]
  mut <- genotyping_tables.tib$Mutation[n]
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

# Then add the information to the Excel file. Plots are generated in 5_Stats/5.4_GoT_XV-seq.R