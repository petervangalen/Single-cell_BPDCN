# Peter van Galen, 220902
# Relate gene expression in non-malignant pDCs from patients without marrow involvement vs. healthy donors

library(tidyverse)
library(Seurat)
library(readxl)
library(ggforce)
library(cowplot)
library(limma)

rm(list=ls())

# Set working directory
setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/09_pDC_Expr")

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")
donor_colors <- popcol.tib$hex[23:40]
names(donor_colors) <- popcol.tib$pop[23:40]
group_colors <- popcol.tib$hex[41:43]
names(group_colors) <- popcol.tib$pop[41:43]
mut_colors <- popcol.tib$hex[44:46]
names(mut_colors) <- popcol.tib$pop[44:46]

# Load pDC gene expression and metadata and signature
pdcs <- readRDS("pdcs.rds")
bpdcn_sign <- read.table("bpdcn_sign.txt")[,1]

# Make orig.ident2 into levels
ordered_levels <- c(setdiff(unique(pdcs$orig.ident2), unique(pdcs$orig.ident)), levels(pdcs$orig.ident)[-1])
all(ordered_levels %in% pdcs$orig.ident2); all(pdcs$orig.ident2 %in% ordered_levels)
pdcs$orig.ident2 <- factor(pdcs$orig.ident2, levels = ordered_levels)

# Add signatures score
pdcs <- AddModuleScore(pdcs, features = list(bpdcn_sign), name = "bpdcn_sign_score")
colnames(pdcs@meta.data) <- gsub("score1$", "score", colnames(pdcs@meta.data))

# Extract metadata
metadata_tib <- as_tibble(pdcs@meta.data, rownames = "cell")

# Select non-malignant pDCs, calculate P-values for all genes
pdcs_nonmalignant <- subset(pdcs, bm_involvement %in% c("HD", "No") & bpdcn_sign_score < 0)
pdcs_nonmalignant@active.ident <- pdcs_nonmalignant$bm_involvement
markerGenes <- FindAllMarkers(pdcs_nonmalignant, logfc.threshold = 0, min.pct = 0, min.cells.feature = 0, return.thresh = 1.01)
# Extract a table of all genes; those with positive fold change are higher in non-involved bone marrow samples
markerGenes_tib <- as_tibble(markerGenes[markerGenes$cluster == "No",])
# Add and organize columns
means_tib <- tibble(gene = rownames(pdcs_nonmalignant),
                    HD.mean = rowMeans(expm1(GetAssayData(subset(pdcs_nonmalignant, bm_involvement == "HD")))),
                    Uninvolved.mean = rowMeans(expm1(GetAssayData(subset(pdcs_nonmalignant, bm_involvement == "No")))))
markerGenes_tib <- markerGenes_tib %>% left_join(means_tib) %>% rename(Uninvolved.pct = pct.1, HD.pct = pct.2) %>%
  select(gene, HD.mean, Uninvolved.mean, HD.pct, Uninvolved.pct, avg_log2FC, p_val, p_val_adj) %>%
  arrange(-avg_log2FC)
write_tsv(file = "diff_genes_nonmalignant_pDCs.txt", markerGenes_tib)

pdf("9.4_Pre-malignant_pDCs.pdf", width = 10, height = 5)
for (g in c(bpdcn_sign, "CXCR4")) {
  print(g)
  #g <- "CXCR4"
  # Add gene expression to metadata
  stopifnot(all(metadata_tib$cell == colnames(pdcs)))
  metadata_tib$Expression <- GetAssayData(pdcs)[g,]
  
  # Select pDCs that are not malignant
  #metadata_tib %>% filter(bm_involvement %in% c("HD", "No")) %>%
  #  ggplot(aes(x = bm_involvement, y = bpdcn_sign_score)) +
  #  geom_point() + geom_hline(yintercept = 0) + theme(aspect.ratio = 2)
  plot_tib <- metadata_tib %>% filter(bm_involvement %in% c("HD", "No"), bpdcn_sign_score < 0)

  p1 <- plot_tib %>% ggplot(aes(x = orig.ident2, y = Expression, color = orig.ident2)) +
    geom_jitter(height = 0) +
    geom_violin(draw_quantiles = 0.5, fill = NA, color = "black") +
    scale_color_manual(values = donor_colors[as.character(unique(plot_tib$orig.ident2))]) +
    ggtitle(g) +
    theme_bw() +
    theme(aspect.ratio = 1, panel.grid = element_blank(), plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
          axis.text.y = element_text(color = "black"), axis.ticks = element_line(color = "black"))
  p2 <- plot_tib %>% ggplot(aes(x = bm_involvement, y = Expression)) +
    geom_sina(aes(color = orig.ident, group = bm_involvement)) +
    geom_violin(draw_quantiles = 0.5, fill = NA) +
    scale_color_manual(values = donor_colors[as.character(unique(plot_tib$orig.ident))]) +
    annotate("text", y = max(plot_tib$Expression), x = 1.5,
             label = paste0("P = ", round(filter(markerGenes_tib, gene == g)$p_val_adj, 4))) +
    ggtitle(g) +
    theme_bw() +
    theme(aspect.ratio = 1, panel.grid = element_blank(), plot.title = element_text(hjust = 0.5),
          axis.text = element_text(color = "black"), axis.ticks = element_line(color = "black"))
  print(
  plot_grid(p1, p2)
  )
}
dev.off()


# One more thing...what does CXCR4 look like in pDCs of all samples?
stopifnot(all(metadata_tib$cell == colnames(pdcs)))
metadata_tib$CXCR4 <- GetAssayData(pdcs)["CXCR4",]

p3 <- metadata_tib %>% ggplot(aes(x = orig.ident2, y = CXCR4, color = orig.ident2)) +
  geom_jitter(height = 0) +
  geom_violin(draw_quantiles = 0.5, fill = NA, color = "black") +
  scale_color_manual(values = donor_colors) +
  ggtitle("CXCR4") +
  theme_bw() +
  theme(aspect.ratio = 1, panel.grid = element_blank(), plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black"), axis.ticks = element_line(color = "black"))

p4 <- metadata_tib %>% ggplot(aes(x = bm_involvement, y = CXCR4)) +
  geom_sina(aes(color = orig.ident, group = bm_involvement)) +
  geom_violin(draw_quantiles = 0.5, fill = NA) +
  scale_color_manual(values = donor_colors) +
  ggtitle(g) +
  theme_bw() +
  theme(aspect.ratio = 1, panel.grid = element_blank(), plot.title = element_text(hjust = 0.5),
        axis.text = element_text(color = "black"), axis.ticks = element_line(color = "black"))

plot_grid(p3, p4)
