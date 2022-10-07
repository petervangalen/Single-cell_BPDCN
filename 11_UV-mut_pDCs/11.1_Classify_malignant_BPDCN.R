# Peter van Galen, 220923
# Classify malignant BPDCN cells

library(tidyverse)
library(Seurat)
library(readxl)
library(ggforce)
library(cowplot)

rm(list=ls())

# Set working directory
setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/11_UV-mut_pDCs")

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")
cell_colors <- popcol.tib$hex[1:21]
names(cell_colors) <- popcol.tib$pop[1:21]
mut_colors <- popcol.tib$hex[44:46]
names(mut_colors) <- popcol.tib$pop[44:46]
group_colors <- popcol.tib$hex[41:43]
names(group_colors) <- popcol.tib$pop[41:43]

# Load Seurat objects
seurat_files <- list.files("../04_XV-seq", pattern = "*.rds", full.names = T)
seu_ls <- lapply(seurat_files, function(x) readRDS(x))
names(seu_ls) <- cutf(basename(seurat_files), d = "_")
seu <- merge(seu_ls[[1]], seu_ls[2:length(seu_ls)])
seu$CellType <- factor(seu$CellType, levels = levels(seu_ls[[1]]$CellType))

# Add malignant BPDCN cell module score
bpdcn_sign <- read.table("../09_pDC_expr/bpdcn_sign.txt")[,1]
seu <- AddModuleScore(seu, features = list(bpdcn_sign), name = "bpdcn_score")
colnames(seu@meta.data) <- gsub("bpdcn_score1", "bpdcn_score", colnames(seu@meta.data))


# Relate bpdcn_score to bone marrow involvement ---------------------------------------------------

metadata_tib <- as_tibble(seu@meta.data, rownames = "cell") %>%
  dplyr::select(cell, orig.ident, CellType, bm_involvement, bpdcn_score)

pdf(file = "11.1.1_Scores_by_BM_involvement.pdf", width = 6, height = 5)
metadata_tib[sample(nrow(metadata_tib)),] %>%
  ggplot(aes(x = bm_involvement, y = bpdcn_score)) +
  geom_sina(aes(color = CellType, group = bm_involvement), size = 0.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  scale_color_manual(values = cell_colors) +
  theme_bw() +
  theme(aspect.ratio = 1, panel.grid = element_blank(), axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"))
dev.off()


# Classify malignant/normal pDCs ------------------------------------------------------------------

# Sina plots of bpdcn_score across cell types colored by bone marrow involvement
seu@meta.data %>% ggplot(aes(x = CellType, y = bpdcn_score, color = bm_involvement)) +
  geom_sina(aes(group = CellType), size = 0.3) +
  scale_color_manual(values = group_colors) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  theme_bw() +
  theme(aspect.ratio = 0.5, panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(), axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"))

# Based on these data, we decided to classify malignant BPDCN cells as:
#    (1) All cells classified as pDCs in patients with bone marrow involvement
#    (2) All cells in any BPDCN patient  (with or without bone marrow involvement) with a BPDCN signature score exceeding 0.5
seu$Malignant_cell <- ifelse(seu$bm_involvement == "Yes" & seu$CellType == "pDC", yes = "Malignant_cell", no = "Normal_cell")
seu$Malignant_cell <- ifelse(seu$bm_involvement != "HD" & seu$bpdcn_score >= 0.5, yes = "Malignant_cell", no = seu$Malignant_cell)
seu$Malignant_cell <- factor(seu$Malignant_cell, levels = c("Normal_cell", "Malignant_cell"))

pdf(file = "11.1.2_Malignant_cells.pdf", width = 12, height = 5)
seu@meta.data %>% ggplot(aes(x = CellType, y = bpdcn_score, color = Malignant_cell)) +
  geom_sina(aes(group = CellType), size = 1) +
  scale_color_manual(values = c("#bdb76b", "#b22222")) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  ggtitle(paste("Classification of all", ncol(seu), "cells")) +
  theme_bw() +
  theme(aspect.ratio = 0.5, panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(), axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5))
dev.off()


# Check correspondence with surface markers -------------------------------------------------------

all(metadata_tib$cell == colnames(seu))
genes <- c("CD19", "IL3RA", "CD4", "SDC1")
for (g in genes) {
  metadata_tib[,g] <- GetAssayData(seu, slot = "data")[g,]
}

pdf(file = "11.1.3_MarkerGenes.pdf", width = 5, height = 4)

p1 <- metadata_tib %>% filter(CellType == "ProB") %>%
  ggplot(aes(x = bpdcn_score > 0.5, y = CD19)) +
  geom_violin(scale = "width", fill = "#eee8aa") +
  geom_jitter(alpha = 0.4, size = 0.7) +
  ggtitle("ProB cells") +
  theme_bw() +
  theme(aspect.ratio = 2, panel.grid = element_blank(), axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"), plot.title = element_text(hjust = 0.5))

p2 <- metadata_tib %>% filter(CellType == "Plasma") %>%
  ggplot(aes(x = bpdcn_score > 0.5, y = SDC1)) +
  geom_violin(scale = "width", fill = "#eee8aa") +
  geom_jitter(alpha = 0.4, size = 0.7) +
  ggtitle("Plasma cells") +
  theme_bw() +
  theme(aspect.ratio = 2, panel.grid = element_blank(), axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"), plot.title = element_text(hjust = 0.5))

plot_grid(p1, p2)

dev.off()

# Some stats
seu@meta.data %>% filter(CellType == "ProB") %>% mutate(bpdcn_exceed = bpdcn_score >= 0.5) %>%
  group_by(CellType, Malignant_cell, bpdcn_exceed) %>%
  summarize(n = n(), min_s = min(bpdcn_score), max_s = max(bpdcn_score))
seu@meta.data %>% filter(CellType == "Plasma") %>% mutate(bpdcn_exceed = bpdcn_score >= 0.5) %>%
  group_by(CellType, Malignant_cell, bpdcn_exceed) %>%
  summarize(n = n(), min_s = min(bpdcn_score), max_s = max(bpdcn_score))
