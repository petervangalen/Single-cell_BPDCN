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
setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/11_pDC_expr")

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")
donor_colors <- popcol.tib$hex[23:40]
names(donor_colors) <- popcol.tib$pop[23:40]
group_colors <- popcol.tib$hex[41:43]
names(group_colors) <- popcol.tib$pop[41:43]
mut_colors <- popcol.tib$hex[44:46]
names(mut_colors) <- popcol.tib$pop[44:46]

# Load pDC gene expression and metadata (including UMAP coordinates) from 11.1_pDC_object.R
pdcs <- readRDS("pdcs.rds")
# Add information from previous script
MalignantCalls_df <- read.table("11.3_MalignantCalls_Final.txt", header = T, row.names = "cell")
all(colnames(pdcs) %in% rownames(MalignantCalls_df))
pdcs <- AddMetaData(pdcs, MalignantCalls_df[colnames(pdcs), c("bpdcn_sign_score", "RF_pDC_score", "is_malignant")])
pdcs$is_malignant <- factor(pdcs$is_malignant, levels = c("Healthy", "Premalignant", "Malignant"))


# Calculate differentially expressed genes --------------------------------------------------------

# Select non-malignant pDCs, calculate P-values for all genes
pdcs_nonmalignant <- subset(pdcs, is_malignant != "Malignant")
pdcs_nonmalignant@active.ident <- pdcs_nonmalignant$is_malignant
markerGenes <- FindAllMarkers(pdcs_nonmalignant, logfc.threshold = 0, min.pct = 0, min.cells.feature = 0, return.thresh = 1.01)
# Extract a table of all genes; those with positive fold change are higher in putative pre-malignant pDCs
markerGenes_tib <- as_tibble(markerGenes[markerGenes$cluster == "Premalignant",])
# Add and organize columns
means_tib <- tibble(gene = rownames(pdcs_nonmalignant),
                    healthy_mean = rowMeans(expm1(GetAssayData(subset(pdcs_nonmalignant, bm_involvement == "HD")))),
                    premalignant_mean = rowMeans(expm1(GetAssayData(subset(pdcs_nonmalignant, bm_involvement == "No")))))
markerGenes_tib <- markerGenes_tib %>% left_join(means_tib) %>% rename(premalignant.pct = pct.1, healthy.pct = pct.2) %>%
  select(gene, healthy_mean, premalignant_mean, healthy.pct, premalignant.pct, avg_log2FC, p_val, p_val_adj) %>%
  arrange(-avg_log2FC)
interesting_genes <- markerGenes_tib %>% filter(abs(avg_log2FC) > log2(1.5), p_val_adj < 0.05) %>% .$gene
write_tsv(file = "11.5_Premalignant-vs-Healthy_pDC_DEG.txt", markerGenes_tib)

# Calculate P-values for healthy vs. malignant pDCs (only to annotate the plots below).
pdc_subset <- subset(pdcs, is_malignant != "Premalignant")
pdc_subset@active.ident <- pdc_subset$is_malignant
markerGenes2 <- FindAllMarkers(pdc_subset, logfc.threshold = 0, min.pct = 0, min.cells.feature = 0, return.thresh = 1.01)
# Extract a table of all genes; those with positive fold change are higher in malignant cells
markerGenes2_tib <- as_tibble(markerGenes2[markerGenes2$cluster == "Malignant",])

# Select genes to highlight
bpdcn_sign <- read.table("bpdcn_sign.txt")[,1]
interesting_genes <- c("CXCR3", "CXCR4", "LILRA4", "PLCG2", interesting_genes, bpdcn_sign)

# Plot
metadata_tib <- as_tibble(pdcs@meta.data, rownames = "cell")
metadata_tib$orig.ident2 <- factor(metadata_tib$orig.ident2, levels = c(paste0("BM", 1:6),
  "Pt1Rem", "Pt5Dx", "Pt9Dx", "Pt10Dx", "Pt12Dx", "Pt1Dx", "Pt10Rel", "Pt12Rel", "Pt14Dx", "Pt15Dx", "Pt16Dx"))

#options(warn = 1) # I don't know why but I had to add this on 221112 to make it work

pdf("11.5_Premalignant_pDCs.pdf", width = 12, height = 4)
for (g in interesting_genes) {
  print(g)
  #g <- "CXCR3"
  
  stopifnot(all(metadata_tib$cell == colnames(pdcs)))
  metadata_tib$g <- GetAssayData(pdcs)[g,]
  
  p1 <- metadata_tib %>%
    ggplot(aes(x = orig.ident2, y = g, color = is_malignant)) +
    geom_jitter(height = 0) +
    geom_violin(draw_quantiles = 0.5, scale = "width", fill = NA, color = "black") +
    scale_color_manual(values = c("#00ff7f", "#4682b4", "#ff6347" )) +
    ylab(paste(g, "expression")) +
    theme_bw() +
    theme(aspect.ratio = 1, panel.grid = element_blank(), plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
          axis.text.y = element_text(color = "black"), axis.ticks = element_line(color = "black"),
          axis.title.x = element_blank())
  
  p2 <- metadata_tib %>% mutate(my_order = sample(nrow(metadata_tib))) %>% arrange(my_order) %>%
    ggplot(aes(x = is_malignant, y = g)) +
    geom_sina(aes(color = orig.ident, group = is_malignant), scale = "width") +
    geom_violin(draw_quantiles = 0.5, scale = "width", fill = NA) +
    scale_color_manual(values = donor_colors[levels(metadata_tib$orig.ident)]) +
    ylab(paste(g, "expression")) +
    annotate("text", y = max(metadata_tib$g), x = 2,
             label = paste0("P = ", round(filter(markerGenes_tib, gene == g)$p_val_adj, 4))) +
    annotate("text", y = max(metadata_tib$g), x = 3,
             label = paste0("P = ", round(filter(markerGenes2_tib, gene == g)$p_val_adj, 4))) +
    theme_bw() +
    theme(aspect.ratio = 1, panel.grid = element_blank(), plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
          axis.text.y = element_text(color = "black"), axis.ticks = element_line(color = "black"),
          axis.title.x = element_blank())
  
  print(
  plot_grid(p1, p2)
  )
}
dev.off()

