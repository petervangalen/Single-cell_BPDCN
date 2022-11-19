# Peter van Galen, 221002
# What TCR genes are in the gene expression matrix?

library(tidyverse)
library(Seurat)
library(readxl)
library(ggforce)

rm(list=ls())

# Set working directory
setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/13_Misc")

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")
cell_colors <- popcol.tib$hex[1:21]
names(cell_colors) <- popcol.tib$pop[1:21]

# Load Seurat healthy donor BM object
seurat_files <- list.files("../04_XV-seq", pattern = "*.rds", full.names = T)
seu <- readRDS("../04_XV-seq/BM_Seurat_Final.rds")

# This may be too simplistic
tcr_genes <- rownames(GetAssayData(seu))[grepl("^TRAV|^TRAC|^TRBV|^TRBC1", rownames(GetAssayData(seu)))]
write_tsv(tibble(symbol = tcr_genes), file = "TCR_genes.txt")

# Add expression to metadata
seu$GAPDH <- GetAssayData(seu, slot = "data")["GAPDH",]
seu$TCR_union <- colSums(GetAssayData(seu, slot = "data")[tcr_genes,])

pdf("TCR_expr.pdf", height = 4, width = 10)
as_tibble(seu@meta.data) %>% select(CellType, GAPDH, TCR_union) %>%
  pivot_longer(cols = c(GAPDH, TCR_union)) %>%
  ggplot(aes(x = CellType, y = value, color = CellType)) +
  #geom_sina(size = 0.3) +
  geom_boxplot(outlier.size = 0.1) +
  scale_color_manual(values = cell_colors) +
  facet_wrap(~name) +
  ylab("Normalized expression") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none")
dev.off()
