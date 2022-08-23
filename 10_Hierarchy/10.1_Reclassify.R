# Peter van Galen, 220718
# Reclassify all cells from BPDCN samples with an additional malignant BPDCN cell classs

library(tidyverse)
library(Seurat)
library(readxl)
library(data.table)
library(randomForest)
library(gplots)
library(ggforce)
library(cowplot)
#library(ggrepel)
#library(viridis)

rm(list=ls())

# Set working directory
setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/10_Hierarchy")

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")
cell_colors <- popcol.tib$hex[1:22]
names(cell_colors) <- popcol.tib$pop[1:22]
donor_colors <- popcol.tib$hex[24:41]
names(donor_colors) <- popcol.tib$pop[24:41]
group_colors <- popcol.tib$hex[42:44]
names(group_colors) <- popcol.tib$pop[42:44]
mut_colors <- popcol.tib$hex[45:47]
names(mut_colors) <- popcol.tib$pop[45:47]

# Load Seurat objects
seurat_files <- list.files("../04_XV-seq", pattern = "*.rds", full.names = T)
seu_ls <- lapply(seurat_files, function(x) readRDS(x))
names(seu_ls) <- cutf(basename(seurat_files), d = "_")
seu <- merge(seu_ls[[1]], seu_ls[2:length(seu_ls)])
seu$orig.ident2 <- ifelse(grepl("BM", seu$orig.ident), yes = cutf(seu$replicate, d = "\\."), no = seu$orig.ident)
seu$CellType <- factor(seu$CellType, levels = levels(seu_ls[[1]]$CellType))

# Load genotyping information
genotyping_tables.tib <- read_excel("../04_XV-seq/XV-seq_overview.xlsx")
# Replace different MTAP entries with one, just as in 4.1_Add_GoT-XV_to_Seurat.R
genotyping_tables.tib$Mutation <- gsub("MTAP.rearr.*", "MTAP.rearr", genotyping_tables.tib$Mutation)

# Add malignant BPDCN cell module score, extract metadata
bpdcn_sign <- read.table("../09_pDC_expr/bpdcn_sign.txt")[,1]
seu <- AddModuleScore(seu, features = list(bpdcn_sign), name = "bpdcn_score")
metadata_tib <- as_tibble(seu@meta.data, rownames = "cell") %>%
  dplyr::select(cell, orig.ident, project.umap.x, project.umap.y, genotyping_tables.tib$Mutation,
                CellType, bm_involvement, bpdcn_score1) %>%
  rename(bpdcn_score = "bpdcn_score1")


# Add genotyping information ----------------------------------------------------------------------

# Add which cells have founder mutations
f_muts <- filter(genotyping_tables.tib, `Founder or progression mutation` == "Founder") %>% .$Mutation
metadata_tib$f_call <- ifelse(apply(dplyr::select(metadata_tib, all_of(f_muts)), 1, 
                                    function(x) sum(grepl("wildtype", x) > 0)), yes = "wildtype", no = "no call")
metadata_tib$f_call <- ifelse(apply(dplyr::select(metadata_tib, all_of(f_muts)), 1,
                                    function(x) sum(grepl("mutant", x) > 0)), yes = "mutant", no = metadata_tib$f_call)
metadata_tib$f_call <- factor(metadata_tib$f_call, levels = c("no call", "wildtype", "mutant"))
# Add which cells have progession mutations
p_muts <- filter(genotyping_tables.tib, `Founder or progression mutation` == "Progression") %>% .$Mutation
metadata_tib$p_call <- ifelse(apply(dplyr::select(metadata_tib, all_of(p_muts)), 1, 
                                    function(x) sum(grepl("wildtype", x) > 0)), yes = "wildtype", no = "no call")
metadata_tib$p_call <- ifelse(apply(dplyr::select(metadata_tib, all_of(p_muts)), 1,
                                    function(x) sum(grepl("mutant", x) > 0)), yes = "mutant", no = metadata_tib$p_call)
metadata_tib$p_call <- factor(metadata_tib$p_call, levels = c("no call", "wildtype", "mutant"))


# Add a new cell type (BPDCN) and reclassify cells accordingly ------------------------------------

# Add new cell type: BPDCN
metadata_tib <- metadata_tib %>% mutate(CellTypeRefined = case_when(
  CellType == "pDC" & bm_involvement == "Yes" ~ "BPDCN",
  bpdcn_score > 0.5 ~ "BPDCN",
  TRUE ~ CellType)) %>%
  mutate(CellTypeRefined = factor(CellTypeRefined, levels = c(levels(metadata_tib$CellType), "BPDCN")))
metadata_tib$CellType %>% table
metadata_tib$CellTypeRefined %>% table

# Visualize BPDCN score, split by patient bone marrow involvement
plot_tib <- metadata_tib[sample(nrow(metadata_tib)),]

pdf(file = "10.1_Scores_by_BM_involvement.pdf", width = 6, height = 5)
plot_tib %>%
  ggplot(aes(x = bm_involvement, y = bpdcn_score)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  geom_sina(aes(color = CellType, group = bm_involvement), size = 0.5) +
  scale_color_manual(values = cell_colors) +
  theme_bw() +
  theme(aspect.ratio = 1, panel.grid = element_blank(), axis.text = element_text(color = "black"))
dev.off()

# Visualize reclassification and mutations
pdf(file = "10.2_Reclassify.pdf", width = 8, heigh = 8)

# Before reclassification
p1 <- metadata_tib %>% arrange(p_call) %>%
  ggplot(aes(x = CellType, y = bpdcn_score)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  geom_sina(aes(color = p_call, group = CellType), size = 0.3) +
  scale_color_manual(values = mut_colors) +
  theme_bw() +
  ggtitle("Before reclassification") +
  theme(aspect.ratio = 0.5, panel.grid = element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank(), axis.text = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5))

# After reclassification
p2 <- metadata_tib %>% arrange(p_call) %>%
  ggplot(aes(x = CellTypeRefined, y = bpdcn_score)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  geom_sina(aes(color = p_call, group = CellTypeRefined), size = 0.3) +
  scale_color_manual(values = mut_colors) +
  theme_bw() +
  ggtitle("After reclassification") +
  theme(aspect.ratio = 0.5, panel.grid = element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank(), axis.text = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5))

plot_grid(p1, p2, ncol = 1)

dev.off()

# Check accuracy of reclassification by evaluating expression of canonical marker genes -----------

# Add gene epxression to metadata
all(metadata_tib$cell == colnames(seu))
genes <- c("CD19", "IL3RA", "CD4", "SDC1")
for (g in genes) {
  metadata_tib[,g] <- GetAssayData(seu, slot = "data")[g,]
}

pdf(file = "10.3_MarkerGenes.pdf", width = 5, height = 4)

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


