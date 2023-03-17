# Peter van Galen, 221213
# Quickly visualize expression of a gene

library(tidyverse)
library(Seurat)
library(readxl)
#library(data.table)
#library(ggrepel)
#library(ggforce)
#library(cowplot)
#library(viridis)

rm(list=ls())

# Set working directory
setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/13_Misc")

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")

# Load Seurat objects
seurat_files <- list.files("../04_XV-seq", pattern = "*.rds", full.names = T)
seu_ls <- lapply(seurat_files, function(x) readRDS(x))
names(seu_ls) <- cutf(basename(seurat_files), d = "_")
seu <- merge(seu_ls[[1]], seu_ls[2:length(seu_ls)])
seu$orig.ident2 <- ifelse(grepl("BM", seu$orig.ident), yes = cutf(seu$replicate, d = "\\."), no = seu$orig.ident)

# Add metadata from 11.3_Classify_malignant_BPDCN.R
MalignantCalls_df <- read.table("../11_pDC_expr/11.3_MalignantCalls_Final.txt", header = T, row.names = "cell")
all(rownames(seu@meta.data) == rownames(MalignantCalls_df))
seu$is_malignant <- MalignantCalls_df$is_malignant

# Load genotyping information
genotyping_tables.tib <- read_excel("../04_XV-seq/XV-seq_overview.xlsx")
# Replace different MTAP entries with one, just as in 4.1_Add_GoT-XV_to_Seurat.R
genotyping_tables.tib$Mutation <- gsub("MTAP.rearr.*", "MTAP.rearr", genotyping_tables.tib$Mutation)


GetAssayData(seu, slot = "data")["MALAT1",]
sum(grepl("MALAT1", rownames(seu)))

# Evidently this script is unfinished. To be continued.
subset(seu, orig.ident == "BM")
FeaturePlot(seu_ls[[1]], features = c("CD80", "CD86")) + theme(aspect.ratio = 1)
FeaturePlot(seu_ls[[1]], features = c("SSBP1")) + theme(aspect.ratio = 1)









