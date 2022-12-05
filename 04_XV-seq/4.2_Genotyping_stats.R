# Peter van Galen, 220522
# Determine genotyping efficiency of every mutation

library(tidyverse)
library(Seurat)
library(readxl)
library(data.table)

rm(list=ls())

# Set working directory
setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/04_XV-seq")

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")

# Load Seurat objects
seurat_files <- list.files("../04_XV-seq", pattern = "*.rds", full.names = T)
seu_bpdcn.ls <- lapply(seurat_files, function(x) readRDS(x))
names(seu_bpdcn.ls) <- gsub("_Seurat_Final.rds", "", cutf(seurat_files, d = "/", f = 3))

# Generate data frame for genotyping results
genotyping_tables.tib <- read_excel("../04_XV-seq/XV-seq_overview.xlsx")

# Determine the number of mutated and wildtype UMIs -----------------------------------------------

# Then, add the information to the Excel file "XV-seq_overview.xlsx"
summ_tables_ls <- lapply(genotyping_tables.tib$FilteredCells, read_tsv)
genotyping_tables.tib$wtUMIs <- unlist(lapply(summ_tables_ls, function(x) sum(x$"wtUMIs")))
genotyping_tables.tib$mutUMIs <- unlist(lapply(summ_tables_ls, function(x) sum(x$"mutUMIs")))
write_tsv(genotyping_tables.tib[,c("wtUMIs", "mutUMIs")], file = "4.2_UMI_counts.txt")


# Integrate mutation data, quantify genotyping efficiency and mutated cell fraction ----------------

# Replace different MTAP entries with one, just as in 4.1_Add_GoT-XV_to_Seurat.R
genotyping_tables.tib$Mutation <- gsub("MTAP.rearr.*", "MTAP.rearr", genotyping_tables.tib$Mutation)
genotyping_tables.tib <- genotyping_tables.tib %>% select(Sample, Mutation) %>% unique

# Set up new columns
genotyping_tables.tib$Efficiency <- NA
genotyping_tables.tib$GenotypedCells <- NA
genotyping_tables.tib$MutCellFraction <- NA

for (n in 1:nrow(genotyping_tables.tib)) {
  #n <- 6
  Sample <- genotyping_tables.tib$Sample[n]
  mut <- genotyping_tables.tib$Mutation[n]
  seu <- seu_bpdcn.ls[[Sample]]
  
  # Extract the necessary columns from Seurat metadata
  calls <- seu@meta.data[,mut]
  
  # Extract numbers
  n_wt <- sum(grepl("wildtype", calls))
  n_mut <- sum(grepl("mutant", calls))
  no_call <- sum(grepl("no call", calls))
  total <- length(calls)
  stopifnot(n_wt + n_mut + no_call == total)
  
  genotyping_tables.tib$Efficiency[n] <- (n_wt+n_mut)/total*100
  genotyping_tables.tib$GenotypedCells[n] <- n_wt+n_mut
  genotyping_tables.tib$MutCellFraction[n] <- n_mut/(n_wt+n_mut)*100
  
}

write_tsv(genotyping_tables.tib, file = "4.2_Genotyping_stats.txt")

# Then add the information to the Excel file "XV-seq_overview.xlsx". Plots are generated in 5_Stats/5.4_GoT_XV-seq.R