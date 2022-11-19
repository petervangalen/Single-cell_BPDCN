# Peter van Galen, 221103
# Classify malignant BPDCN cells

library(tidyverse)
library(Seurat)
library(readxl)
library(ggforce)
library(cowplot)

rm(list=ls())

# Set working directory
setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/09_pDC_BPDCN")

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
bpdcn_sign <- read.table("bpdcn_sign.txt")[,1]
seu <- AddModuleScore(seu, features = list(bpdcn_sign), name = "bpdcn_sign_score")
colnames(seu@meta.data) <- gsub("score1$", "score", colnames(seu@meta.data))


# Relate bpdcn_sign_score to bone marrow involvement ---------------------------------------------------

pdf(file = "9.4.1_Scores_by_BM_involvement.pdf", width = 6, height = 5)
seu@meta.data %>% mutate(my_order = sample(ncol(seu))) %>% arrange(my_order) %>%
  ggplot(aes(x = bm_involvement, y = bpdcn_sign_score)) +
  geom_sina(aes(color = CellType, group = bm_involvement), size = 0.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  scale_color_manual(values = cell_colors) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  theme_bw() +
  theme(aspect.ratio = 1, panel.grid = element_blank(), axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"))
dev.off()


# Classify healthy, premalignant, and malignant pDCs ----------------------------------------------

# Sina plots of bpdcn_sign_score across cell types colored by bone marrow involvement
seu@meta.data %>% mutate(my_order = sample(ncol(seu))) %>% arrange(my_order) %>%
  ggplot(aes(x = CellType, y = bpdcn_sign_score, color = bm_involvement)) +
  geom_sina(aes(group = CellType), size = 0.3) +
  scale_color_manual(values = group_colors) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  theme_bw() +
  theme(aspect.ratio = 0.5, panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(), axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"))

# Classify cells as healthy/premalignant/malignant cells:
# Healthy = cells annotated as pDC in the healthy donors
# Malignant = cells classified as pDC in patients with bone marrow involvement and cells with a BPDCN signature score exceeding 0.5
# Premalignant = cells classified as pDC in patients without bone marrow involvement with a BPDCN signature score lower than 0.5
seu$is_malignant <- ifelse(seu$bm_involvement == "HD" & seu$CellType == "pDC", yes = "Healthy", no = "Other")
seu$is_malignant <- ifelse(seu$bm_involvement == "Yes" & seu$CellType == "pDC", yes = "Malignant", no = seu$is_malignant)
seu$is_malignant <- ifelse(seu$bpdcn_sign_score > 0.5, yes = "Malignant", no = seu$is_malignant)
seu$is_malignant <- ifelse(seu$bm_involvement == "No" & seu$CellType == "pDC" & seu$bpdcn_sign_score < 0.5,
                           yes = "Premalignant", no = seu$is_malignant)
# These classifications are backed up by mutation analysis elsewhere
seu$is_malignant <- factor(seu$is_malignant, levels = c("Healthy", "Premalignant", "Malignant", "Other"))

# Check / numbers for EDF7
as_tibble(seu@meta.data) %>% filter(CellType == "pDC") %>% count(orig.ident, CellType, bm_involvement, is_malignant)
as_tibble(seu@meta.data) %>% .$is_malignant %>% table(useNA = "always")
as_tibble(seu@meta.data) %>% filter(CellType == "pDC", bm_involvement == "Yes") %>% .$is_malignant %>% table
as_tibble(seu@meta.data) %>% filter(CellType == "pDC", bm_involvement == "No") %>% .$is_malignant %>% table
as_tibble(seu@meta.data) %>% filter(CellType != "pDC", bm_involvement == "Yes") %>% .$is_malignant %>% table

pdf(file = "9.4.2_Classify_Premalignant_Malignant.pdf", width = 12, height = 5)

mycol <- c(Healthy = "#bdb76b", Premalignant = "#ff69b4", Malignant = "#b22222", Other = "#d3d3d3")

as_tibble(seu@meta.data) %>%
  mutate(is_malignant = factor(is_malignant, levels = c("Malignant", "Premalignant", "Healthy", "Other"))) %>%
  arrange(is_malignant) %>%
  ggplot(aes(x = CellType, y = bpdcn_sign_score, color = is_malignant)) +
  geom_sina(aes(group = CellType), size = 0.5) +
  scale_color_manual(values = mycol) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  ggtitle(paste("Classification of all", ncol(seu), "cells")) +
  theme_bw() +
  theme(aspect.ratio = 0.5, panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(), axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5))

dev.off()


# Check correspondence with surface markers -------------------------------------------------------

metadata_tib <- as_tibble(seu@meta.data, rownames = "cell")

all(metadata_tib$cell == colnames(seu))
genes <- c("CD19", "IL3RA", "CD4", "SDC1")
for (g in genes) {
  metadata_tib[,g] <- GetAssayData(seu, slot = "data")[g,]
}

pdf(file = "9.4.3_MarkerGenes.pdf", width = 5, height = 4)

p1 <- metadata_tib %>% filter(CellType == "ProB") %>%
  mutate(is_malignant = ifelse(is_malignant == "Malignant", yes = "Reclassified", no = "ProB")) %>%
  ggplot(aes(x = is_malignant, y = CD19)) +
  geom_violin(scale = "width", fill = "#eee8aa") +
  geom_jitter(alpha = 0.4, size = 0.7) +
  ggtitle("ProB cells") +
  theme_bw() +
  theme(aspect.ratio = 2, panel.grid = element_blank(), axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"), plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank())

p2 <- metadata_tib %>% filter(CellType == "Plasma") %>%
  mutate(is_malignant = ifelse(is_malignant == "Malignant", yes = "Reclassified", no = "Plasma")) %>%
  ggplot(aes(x = is_malignant, y = SDC1)) +
  geom_violin(scale = "width", fill = "#eee8aa") +
  geom_jitter(alpha = 0.4, size = 0.7) +
  ggtitle("Plasma cells") +
  theme_bw() +
  theme(aspect.ratio = 2, panel.grid = element_blank(), axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"), plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank())

plot_grid(p1, p2)

dev.off()

# Some stats and save -----------------------------------------------------------------------------

seu@meta.data %>% filter(CellType == "ProB") %>% mutate(bpdcn_exceed = bpdcn_sign_score >= 0.5) %>%
  group_by(CellType, is_malignant, bpdcn_exceed) %>%
  summarize(n = n(), min_s = min(bpdcn_sign_score), max_s = max(bpdcn_sign_score))
seu@meta.data %>% filter(CellType == "Plasma") %>% mutate(bpdcn_exceed = bpdcn_sign_score >= 0.5) %>%
  group_by(CellType, is_malignant, bpdcn_exceed) %>%
  summarize(n = n(), min_s = min(bpdcn_sign_score), max_s = max(bpdcn_sign_score))

# Retrieve random forest prediction scores to plot on x-axis in later analyses
prediction_scores_ls <- readRDS("../03_RandomForest/Prediction_scores.rds")
pDC_scores_tib <- as_tibble(do.call(rbind, lapply(prediction_scores_ls, function(x) data.frame(x[,"pDC"]))), rownames = "cell")
pDC_scores_tib <- pDC_scores_tib %>% rename(RF_pDC_score = "x....pDC..") %>% mutate(cell = gsub("^.*?\\.", "", cell))

# Merge & save
metadata2_tib <- metadata_tib %>% left_join(pDC_scores_tib) %>%
  dplyr::select(cell, orig.ident, CellType, bpdcn_sign_score, RF_pDC_score, is_malignant)
# Quickly visualize (don't include normal cells; they don't have an RF_pDC_score)
metadata2_tib %>% filter(orig.ident != "BM") %>%
  ggplot(aes(x = RF_pDC_score, y = bpdcn_sign_score, color = is_malignant)) +
  geom_point() +
  scale_color_manual(values = mycol) +
  theme_bw() +
  theme(aspect.ratio = 1, panel.grid = element_blank(),
        axis.ticks = element_line(color = "black"), axis.text = element_text(color="black"))
# It looks like a few healthy pDCs from patients with bone marrow involvement are erroneously classified as malignant. This is a limitation of the study that does not affect subsequent analyses or conclusions.

write_tsv(metadata2_tib, file = "9.4_MalignantCalls_Final.txt")

