# Peter van Galen, 220902
# Relate gene expression in non-malignant pDCs from patients without marrow involvement vs. healthy donors

library(tidyverse)
library(Seurat)
library(readxl)
library(ggforce)
library(cowplot)
library(limma)
library(ggrepel)

rm(list=ls())

# Set working directory
setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/11_pDC_expr")

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")
donor_colors <- popcol.tib$hex[23:40]
names(donor_colors) <- popcol.tib$pop[23:40]

# Load all data
seurat_files <- list.files("../04_XV-seq", pattern = "*.rds", full.names = T)
seu_ls <- lapply(seurat_files, function(x) readRDS(x))
names(seu_ls) <- cutf(basename(seurat_files), d = "_")
seu <- merge(seu_ls[[1]], seu_ls[2:length(seu_ls)])
seu$orig.ident2 <- ifelse(grepl("BM", seu$orig.ident), yes = cutf(seu$replicate, d = "\\."), no = seu$orig.ident)

# Add metadata from 11.3_Classify_malignant_BPDCN.R
MalignantCalls_df <- read.table("11.3_MalignantCalls_Final.txt", header = T, row.names = "cell")
all(rownames(seu@meta.data) == rownames(MalignantCalls_df))
seu$is_malignant <- MalignantCalls_df$is_malignant


# Subset cells ------------------------------------------------------------------------------------

# Similar to 11.1_pDC_object.R/11.2_pDC_expression.R, subset for cells of interest and exclude RPS and RPL genes
seu_subset <- subset(seu, is_malignant != "Other")
rp.ch <- rownames(seu)[grepl("^RPS|^RPL", rownames(seu))]
seu_subset <- subset(seu_subset, features = setdiff(rownames(seu_subset), rp.ch))
# Renormalize so the columns add up to 10,000 counts (check using sum(expm1(GetAssayData(seu_subset)[,1])))
seu_subset <- NormalizeData(seu_subset)
seu_subset$is_malignant <- factor(seu_subset$is_malignant, levels = c("Healthy", "Premalignant", "Malignant"))

# The number of pDCs is vastly different per sample, which could skew the analysis
as_tibble(seu_subset@meta.data) %>% dplyr::select(orig.ident, orig.ident2, bm_involvement) %>% group_by_all() %>% count()

# This could be resolved by subsetting for up to 50 cells per sample. Ultimately, I decided not do do this.
#H_cells <- as_tibble(subset(seu, is_malignant == "Healthy")@meta.data, rownames = "cell") %>%
#  group_by(orig.ident2) %>% slice_sample(n = 50) %>% .$cell
#P_cells <- as_tibble(subset(seu, is_malignant == "Premalignant")@meta.data, rownames = "cell") %>%
#  group_by(orig.ident2) %>% slice_sample(n = 50) %>% .$cell
#M_cells <- as_tibble(subset(seu, is_malignant == "Malignant")@meta.data, rownames = "cell") %>%
#  group_by(orig.ident2) %>% slice_sample(n = 50) %>% .$cell
#seu_subset <- subset(seu_subset, cells = c(H_cells, P_cells, M_cells))


# Calculate differentially expressed genes --------------------------------------------------------

# You can skip to the visualizaton section by loading the summary table
#dge_tib <- read_tsv("11.5_Healthy.pDC-Premalignant-Malignant_DEG.txt")

# Calculate log2FC and P-values for Healthy vs. Premalignant (HP) pDCs
nonmalignant <- subset(seu_subset, is_malignant != "Malignant")
nonmalignant@active.ident <- nonmalignant$is_malignant
markerGenes <- FindAllMarkers(nonmalignant, logfc.threshold = 0, min.pct = 0, min.cells.feature = 0, return.thresh = 1.01)
# Extract a table of all genes; those with positive fold change are higher in putative pre-malignant pDCs
markerGenes_tib <- as_tibble(markerGenes[markerGenes$cluster == "Premalignant",])
# Add and organize columns
means_tib <- tibble(gene = rownames(nonmalignant),
                    H_mean = rowMeans(expm1(GetAssayData(subset(nonmalignant, bm_involvement == "HD")))),
                    P_mean = rowMeans(expm1(GetAssayData(subset(nonmalignant, bm_involvement == "No")))))
# Combine relevant information into one table
dge_tib <- markerGenes_tib %>% left_join(means_tib) %>%
  rename(H_pct = pct.2, P_pct = pct.1, HP_avg_log2FC = avg_log2FC, HP_p_val = p_val, HP_p_val_adj = p_val_adj) %>%
  select(gene, H_mean, P_mean, H_pct, P_pct, HP_avg_log2FC, HP_p_val, HP_p_val_adj) %>%
  arrange(-HP_avg_log2FC)

# Calculate log2FC and P-values for Healthy vs. Malignant (HM) pDCs
healthy_malignant <- subset(seu_subset, is_malignant != "Premalignant")
healthy_malignant@active.ident <- healthy_malignant$is_malignant
markerGenes2 <- FindAllMarkers(healthy_malignant, logfc.threshold = 0, min.pct = 0, min.cells.feature = 0, return.thresh = 1.01)
# Extract a table of all genes; those with positive fold change are higher in malignant cells
markerGenes2_tib <- as_tibble(markerGenes2[markerGenes2$cluster == "Malignant",])
# Add these data malignant pDC numbers to the summary table
means2_tib <- tibble(gene = rownames(healthy_malignant),
                     H_mean = rowMeans(expm1(GetAssayData(subset(healthy_malignant, bm_involvement == "HD")))),
                     M_mean = rowMeans(expm1(GetAssayData(subset(healthy_malignant, bm_involvement == "Yes")))))
# Check that this is redundant
stopifnot(means2_tib$H_mean == means_tib$H_mean)
means2_tib$H_mean <- NULL
# Combine relevant information into one table
dge_tib <- dge_tib %>% left_join(select(markerGenes2_tib, gene, pct.1, avg_log2FC, p_val, p_val_adj)) %>%
  left_join(means2_tib) %>%
  rename(M_pct = pct.1, HM_avg_log2FC = avg_log2FC, HM_p_val = p_val, HM_p_val_adj = p_val_adj) %>%
  relocate(gene, H_mean, P_mean, M_mean, H_pct, P_pct, M_pct, HP_avg_log2FC, HP_p_val, HP_p_val_adj,
           HM_avg_log2FC, HM_p_val, HM_p_val_adj)
# Save
write_tsv(file = "11.5_Healthy.pDC-Premalignant-Malignant_DEG.txt", dge_tib)


# Visualize ---------------------------------------------------------------------------------------

# Select genes to highlight
bpdcn_sign <- read.table("bpdcn_sign.txt")[,1]
highlight_genes <- c(slice_max(dge_tib, order_by = HP_avg_log2FC, n = 5)$gene,
                     slice_min(dge_tib, order_by = HP_avg_log2FC, n = 5)$gene,
                     slice_min(dge_tib, order_by = HP_p_val_adj, n = 5)$gene,
                     slice_max(dge_tib, order_by = HM_avg_log2FC, n = 5)$gene,
                     slice_min(dge_tib, order_by = HM_avg_log2FC, n = 5)$gene,
                     slice_min(dge_tib, order_by = HM_p_val_adj, n = 5)$gene,
                     "CXCR4", "BCL2", "TCL1A") %>% sort %>% unique
highlight_genes <- setdiff(highlight_genes, c("HBA1", "HBA2", "HBB", "LYZ", "RGS18", "HLA-DRB1", "F13A1",
                                              "DPPA4", "DACH1", "TRPS1", "VCAN", "S100A8", "S100A9"))
axis_lim <- max(abs(dge_tib$HM_avg_log2FC), abs(dge_tib$HP_avg_log2FC))

# Correlation between premalignant and malignant fold changes
pdf("11.5.1_FoldChange_correlation.pdf")
dge_tib %>%
  ggplot(aes(x = HP_avg_log2FC, y = HM_avg_log2FC, color = gene %in% highlight_genes)) +
  geom_hline(yintercept = 0, color = "grey") +
  geom_vline(xintercept = 0, color = "grey") +
  geom_point() +
  geom_text_repel(data = filter(dge_tib, gene %in% highlight_genes),
                  aes(x = HP_avg_log2FC, y = HM_avg_log2FC, label = gene),
                  color = "black", size = 3, max.overlaps = 30) +
  scale_color_manual(values = c("#708090", "#9acd32")) +
  coord_cartesian(xlim = c(-axis_lim, axis_lim), ylim = c(-axis_lim, axis_lim)) +
  theme_bw() +
  theme(aspect.ratio = 1, axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"), axis.line = element_blank(),
        panel.grid = element_blank())
dev.off()

cor(dge_tib$HP_avg_log2FC, dge_tib$HM_avg_log2FC)

# Volcano for premalignant changes
dge_tib %>% ggplot(aes(x = HP_avg_log2FC, y = -log10(HP_p_val_adj), color = gene %in% highlight_genes)) +
  geom_point() +
  geom_text_repel(data = filter(dge_tib, gene %in% highlight_genes, HP_p_val_adj < 1),
                  aes(x = HP_avg_log2FC, y = -log10(HP_p_val_adj), label = gene),
                  color = "black", size = 3, max.overlaps = 30) +
  scale_color_manual(values = c("#708090", "#9acd32")) +
  coord_cartesian(xlim = c(-max(abs(dge_tib$HP_avg_log2FC)), max(abs(dge_tib$HP_avg_log2FC)))) +
  theme_bw() +
  theme(aspect.ratio = 1, axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"), axis.line = element_blank(),
        panel.grid = element_blank())

# Volcano for malignant changes
dge_tib %>% ggplot(aes(x = HM_avg_log2FC, y = -log10(HM_p_val_adj), color = gene %in% highlight_genes)) +
  geom_point() +
  geom_text_repel(data = filter(dge_tib, gene %in% highlight_genes, HM_p_val_adj < 1),
                  aes(x = HM_avg_log2FC, y = -log10(HM_p_val_adj), label = gene),
                  color = "black", size = 3, max.overlaps = 30) +
  scale_color_manual(values = c("#708090", "#9acd32")) +
  coord_cartesian(xlim = c(-max(abs(dge_tib$HM_avg_log2FC)), max(abs(dge_tib$HM_avg_log2FC)))) +
  theme_bw() +
  theme(aspect.ratio = 1, axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"), axis.line = element_blank(),
        panel.grid = element_blank())

# Sina plots
metadata_tib <- as_tibble(seu_subset@meta.data, rownames = "cell")
metadata_tib$orig.ident2 <- factor(metadata_tib$orig.ident2, levels = c(paste0("BM", 1:6),
  "Pt1Rem", "Pt5Dx", "Pt9Dx", "Pt10Dx", "Pt12Dx", "Pt1Dx", "Pt10Rel", "Pt12Rel", "Pt14Dx", "Pt15Dx", "Pt16Dx"))

#options(warn = 1) # I had to add this on 221112 to make it work

pdf("11.5.2_Premalignant_pDCs.pdf", width = 12, height = 4)
for (g in highlight_genes) {
  print(g)
  #g <- "CXCR3"
  
  stopifnot(all(metadata_tib$cell == colnames(seu_subset)))
  metadata_tib$g <- GetAssayData(seu_subset)[g,]
  
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
    scale_color_manual(values = donor_colors[unique(metadata_tib$orig.ident)]) +
    ylab(paste(g, "expression")) +
    annotate("text", y = max(metadata_tib$g), x = 2,
             label = paste0("P = ", round(filter(dge_tib, gene == g)$HP_p_val_adj, 4))) +
    annotate("text", y = max(metadata_tib$g), x = 3,
             label = paste0("P = ", round(filter(dge_tib, gene == g)$HM_p_val_adj, 4))) +
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

