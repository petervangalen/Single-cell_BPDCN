# Peter van Galen, 221213
# Load all relevant/final data and visualize gene expression

library(tidyverse)
library(Seurat)
library(readxl)
library(janitor)
library(ggforce)
#library(data.table)
#library(ggrepel)
#library(cowplot)
#library(viridis)

rm(list=ls())

# Set working directory
setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/13_Misc")

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")
cell_colors <- popcol.tib$hex[1:21]
names(cell_colors) <- popcol.tib$pop[1:21]
sample_colors <- popcol.tib$hex[23:40]
names(sample_colors) <- popcol.tib$pop[23:40]
group_colors <- popcol.tib$hex[41:43]
names(group_colors) <- popcol.tib$pop[41:43]
mut_colors <- popcol.tib$hex[44:46]
names(mut_colors) <- popcol.tib$pop[44:46]
malignant_stage <- popcol.tib$hex[47:50]
names(malignant_stage) <- popcol.tib$pop[47:50]

# Load data ---------------------------------------------------------------------------------------

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

# Ensure nice ordering in visualizations
seu$orig.ident <- factor(seu$orig.ident, levels = c("BM", "Pt1Dx", "Pt1Rem", "Pt5Dx", "Pt9Dx", "Pt10Dx", "Pt10Rel",
                                                    "Pt12Dx", "Pt12Rel", "Pt14Dx", "Pt15Dx", "Pt16Dx"))
seu$CellType <- factor(seu$CellType, levels = names(cell_colors))

# Plot expression ---------------------------------------------------------------------------------

# Make metadata tibble
metadata_tib <- as_tibble(seu@meta.data, rownames = "cell")

# Add expression values
gene <- "CXCR4"
metadata_tib[,"expression"] <- GetAssayData(seu, slot = "data")[gene,]

# Sina plot example: expression in pDCs across samples
metadata_tib %>% filter(CellType == "pDC") %>%
  ggplot(aes(x = orig.ident, y = expression, color = CellType)) +
  geom_sina(aes(group = orig.ident)) +
  scale_color_manual(values = cell_colors) +
  ggtitle(label = gene) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Healthy donors only, split by cell type
metadata_tib %>% filter(orig.ident == "BM") %>%
  ggplot(aes(x = CellType, y = expression, color = CellType)) +
  geom_sina(aes(group = CellType)) +
  scale_color_manual(values = cell_colors) +
  ggtitle(label = gene) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

  










