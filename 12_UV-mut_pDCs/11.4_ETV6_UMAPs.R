# Peter van Galen, 220914
# Generate UMAPs for MTAP in Patient 10 and ETV6 in Patient 14

library(tidyverse)
library(Seurat)
library(readxl)
library(data.table)
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


# Visualize occurrence of UV-associated mutations -------------------------------------------------

# Choose one of these options
  # Note that most non-pDCs are donor-derived
  #seu <- merge(readRDS("../04_XV-seq/Pt10Dx_Seurat_Final.rds"), readRDS("../04_XV-seq/Pt10Rel_Seurat_Final.rds"))
  #mut <- "MTAP.rearr"

  # This doesn't look good because there are too many tumor cells and suspected false-positives by soup RNA
  #seu <- readRDS("../04_XV-seq/Pt15Dx_Seurat_Final.rds")
  #pt <- "Pt15Dx"
  #mut <- "U2AF1.S34F"

  # This one looks good with the exception of two mutated cells outside of the "malignant cluster". There might be a batch effect between Dx and Rel in the myeloid lineage but I have not looked into that (since this figure will not be published anyway)
  seu <- merge(readRDS("../04_XV-seq/Pt12Dx_Seurat_Final.rds"), readRDS("../04_XV-seq/Pt12Rel_Seurat_Final.rds"))
  pt <- "Pt12"
  mut <- "MALAT1.chr11:65270399:G/A"

  # This one will be in the paper
  seu <- readRDS("../04_XV-seq/Pt14Dx_Seurat_Final.rds")
  pt <- "Pt14Dx"
  mut <- "ETV6.R369W"

# Add metadata column for specific mutation
seu$mutation <- seu@meta.data[,mut]
seu$mutation <- factor(seu$mutation, levels = c("no call", "wildtype", "mutant"))
seu$CellType <- factor(seu$CellType, levels = names(cell_colors))

# Add malignant BPDCN cell module score
bpdcn_sign <- read.table("../09_pDC_expr/bpdcn_sign.txt")[,1]
seu <- AddModuleScore(seu, features = list(bpdcn_sign), name = "bpdcn_score")
colnames(seu@meta.data) <- gsub("bpdcn_score1", "bpdcn_score", colnames(seu@meta.data))

# Use standard Seurat functions (skip harmony and cell cycle regression that were used in 1_Seurat_Harmony.R and 9.1_pDC_object.R)
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu, features = VariableFeatures(object = seu))
seu <- RunUMAP(seu, reduction = "pca", dims = 1:20)

# Extract metadata
seu$UMAP_1 <- seu@reductions$umap@cell.embeddings[,1]
seu$UMAP_2 <-seu@reductions$umap@cell.embeddings[,2]
metadata_tib <- as_tibble(seu@meta.data, rownames = "cell")

# Sina plot of BPDCN scores and mutations
metadata_tib %>% arrange(mutation) %>%
  ggplot(aes(x = CellType, y = bpdcn_score, color = mutation)) +
  geom_sina(aes(group = CellType)) +
  scale_color_manual(values = mut_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Pie charts of the cell type proportions in wild-type and mutant cell populations
mut_tib <- metadata_tib %>% filter(get(mut) != "no call") %>%
  group_by(CellType) %>% summarize(wildtype = sum(get(mut) == "wildtype"),
                                   mutant = sum(get(mut) == "mutant"))
mut_tib %>%
  ggplot(aes(x = "", y = wildtype, fill = CellType)) +
  geom_bar(stat="identity", width=1) +
  scale_fill_manual(values = cell_colors) +
  coord_polar("y") +
  theme_bw() +
  theme(axis.title = element_blank())
mut_tib %>%
  ggplot(aes(x = "", y = mutant, fill = CellType)) +
  geom_bar(stat="identity", width=1) +
  scale_fill_manual(values = cell_colors) +
  coord_polar("y") +
  theme_bw() +
  theme(axis.title = element_blank())

# UMAPs
p1 <- metadata_tib %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = CellType)) +
  geom_point(size = 0.5) +
  ggtitle(label = pt) +
  scale_color_manual(values = cell_colors) +
  theme_bw() +
  theme(aspect.ratio = 1, panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))

p2 <- metadata_tib %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = bpdcn_score)) +
  geom_point(size = 0.5) +
  ggtitle(label = pt) +
  theme_bw() +
  theme(aspect.ratio = 1, panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))

p3 <- metadata_tib %>% mutate(mutation = factor(mutation, levels = c("no call", "wildtype", "mutant"))) %>%
  arrange(mutation) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = mutation)) +
  geom_point(size = 0.5) +
  ggtitle(label = mut) +
  scale_color_manual(values = mut_colors) +
  theme_bw() +
  theme(aspect.ratio = 1, panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))

pdf(paste0("11.1_", mut, "_", pt, "_UMAPs.pdf"), height = 4.5, width = 20)
print(
plot_grid(p1, p2, p3, ncol = 3)
)
dev.off()



