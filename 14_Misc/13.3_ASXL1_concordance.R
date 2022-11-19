# Peter van Galen, 221119
# Test the concordance between ASXL1 mutations which were genotyped from Pt10Dx using two different primers

library(tidyverse)
library(readxl)
library(Seurat)

# Set working directory
setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/14_Misc")

rm(list=ls())

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")

# Load Seurat objects
seurat_files <- list.files("../04_XV-seq", pattern = "*.rds", full.names = T)
seu_ls <- lapply(seurat_files, function(x) readRDS(x))
names(seu_ls) <- cutf(basename(seurat_files), d = "_")

# Load XV-seq data. Of note, the capture of mutated and wild-type transcripts is mutually exclusive (which makes sense because the capture efficiency was very low so the chance of capturing both loci in one cell is very low)
ASXL1.1886 <- read_table("../04_XV-seq/Pt10Dx/ASXL1.1886/ASXL1.1886.FilteredCells.txt")
ASXL1.1898 <- read_table("../04_XV-seq/Pt10Dx/ASXL1.1898/ASXL1.1898.FilteredCells.txt")
# Filter by cells that pass QC
ASXL1.1886 <- ASXL1.1886 %>% filter(BC %in% cutf(colnames(seu_ls[["Pt10Dx"]]), d = "-"))
ASXL1.1898 <- ASXL1.1898 %>% filter(BC %in% cutf(colnames(seu_ls[["Pt10Dx"]]), d = "-"))

# What are the mutated and wildtype cells?
ASXL1.1886_mutated_cells <- ASXL1.1886 %>% filter(mutUMIs > 0) %>% .$BC
ASXL1.1886_wildtype_cells <- ASXL1.1886 %>% filter(wtUMIs > 0) %>% .$BC
ASXL1.1898_mutated_cells <- ASXL1.1898 %>% filter(mutUMIs > 0) %>% .$BC
ASXL1.1898_wildtype_cells <- ASXL1.1898 %>% filter(wtUMIs > 0) %>% .$BC

# Stats for response letter: ""
length(ASXL1.1886_mutated_cells)
length(ASXL1.1898_mutated_cells)
sum(ASXL1.1886_mutated_cells %in% ASXL1.1898_mutated_cells)
sum(ASXL1.1898_mutated_cells %in% ASXL1.1886_mutated_cells)
ncol(seu_ls[["Pt10Dx"]])

length(ASXL1.1886_wildtype_cells)
length(ASXL1.1898_wildtype_cells)
sum(ASXL1.1886_wildtype_cells %in% ASXL1.1898_wildtype_cells)
sum(ASXL1.1898_wildtype_cells %in% ASXL1.1886_wildtype_cells)








