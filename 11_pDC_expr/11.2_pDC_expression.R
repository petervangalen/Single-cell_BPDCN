# Peter van Galen, 220726
# Differential gene expression between pDCs and malignant BPDCN cells

library(tidyverse)
library(Seurat)
library(readxl)
library(data.table)
library(ggrepel)

rm(list=ls())

# Set working directory
setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/11_pDC_expr")

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")
group_colors <- popcol.tib$hex[41:43]
names(group_colors) <- popcol.tib$pop[41:43]

# Load pDC gene expression and metadata
pdcs <- readRDS("pdcs.rds")

# Load genotyping information
genotyping_tables.tib <- read_excel("../04_XV-seq/XV-seq_overview.xlsx")
# Replace different MTAP entries with one, just as in 4.1_Add_GoT-XV_to_Seurat.R
genotyping_tables.tib$Mutation <- gsub("MTAP.rearr.*", "MTAP.rearr", genotyping_tables.tib$Mutation)


# Select cells to compare -------------------------------------------------------------------------

# Define two groups of cells: (1) pDCs from healthy donors and pDCs without progression mutations from skin-only patients, and (2) pDCs from patients with bone marrow involvement
normal_pdcs <- subset(pdcs, cells = colnames(subset(pdcs, bm_involvement != "Yes" & pdcs$progression != "mutant")))
malignant_pdcs <- subset(pdcs, cells = colnames(subset(pdcs, bm_involvement == "Yes")))

# There are widely different numbers of cells for different samples
as_tibble(normal_pdcs@meta.data) %>% dplyr::select(orig.ident, orig.ident2, bm_involvement) %>% group_by_all() %>% count()
as_tibble(malignant_pdcs@meta.data) %>% dplyr::select(orig.ident, orig.ident2, bm_involvement) %>% group_by_all() %>% count()

# Take at most 50 cells from each sample. This may yield an error if you don't use dplyr_1.0.99.9000
normal_pdcs_subset_ids <- as_tibble(normal_pdcs@meta.data, rownames = "cell") %>% group_by(orig.ident2) %>% slice_sample(n = 50) %>% .$cell
malignant_pdcs_subset_ids <- as_tibble(malignant_pdcs@meta.data, rownames = "cell") %>% group_by(orig.ident2) %>% slice_sample(n = 50) %>% .$cell
normal_pdcs_subset <- subset(normal_pdcs, cells = normal_pdcs_subset_ids)
malignant_pdcs_subset <- subset(malignant_pdcs, cells = malignant_pdcs_subset_ids)

# ...that's better
as_tibble(normal_pdcs_subset@meta.data) %>% dplyr::select(orig.ident, orig.ident2, bm_involvement) %>% group_by_all() %>% count()
as_tibble(malignant_pdcs_subset@meta.data) %>% dplyr::select(orig.ident, orig.ident2, bm_involvement) %>% group_by_all() %>% count()


# Differentially expressed genes ------------------------------------------------------------------

# Generate signature
markerGenes <- FindMarkers(pdcs, ident.1 = malignant_pdcs_subset_ids, ident.2 = normal_pdcs_subset_ids,
                           logfc.threshold = 0, min.pct = 0, min.cells.feature = 0, return.thresh = 1.01)
markerGenes_tib <- as_tibble(markerGenes, rownames = "gene") %>% arrange(-avg_log2FC)
write_tsv(markerGenes_tib, file = "11.2_Malignant-vs-Healthy_pDC_DEG.txt")

# Define signatures
bpdcn_sign <- markerGenes_tib %>% filter(p_val_adj < 1E-30, avg_log2FC > 1) %>% .$gene
highlight_genes <- c("IGLL1", "HES6", "SERPINB2", "ARPP21", "PDLIM1", "TCL1A", "SCN3A", "LAMP5",
                     "FCER1A", "MYBPH", "CLEC11A", "SOX4", "GSTP1", "BCL2", "S100A4")
highlight_genes_down <- c("GZMB", "LILRA4", "SERPINF1", "PSAP", "APP", "IRF7", "TGFBI", "RUNX1")

# Volcano plot
pdf("11.2.1_Volcano.pdf")
markerGenes_tib %>%
  ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj), color = gene %in% bpdcn_sign)) +
  geom_point() +
  geom_text_repel(data = filter(markerGenes_tib, gene %in% c(highlight_genes, highlight_genes_down)),
                  aes(x = avg_log2FC, y = -log10(p_val_adj), label = gene),
                  color = "#006400", size = 3, max.overlaps = 30) +
  scale_color_manual(values = c("#708090", "#9acd32")) +
  coord_cartesian(xlim = c(-max(abs(markerGenes_tib$avg_log2FC)), max(abs(markerGenes_tib$avg_log2FC)))) +
  theme_bw() +
  theme(aspect.ratio = 1, axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"), axis.line = element_blank(),
        panel.grid = element_blank())
dev.off()

# How does this compare to Chloe Villani's signature?
pdcs <- AddModuleScore(pdcs, features = list(bpdcn_sign), name = "bpdcn_sign_score")
colnames(pdcs@meta.data) <- gsub("score1$", "score", colnames(pdcs@meta.data))

pdf("11.2.2_Signature_correlation.pdf", width = 5, height = 5)
print(
pdcs@meta.data %>% arrange(fct_rev(bm_involvement)) %>%
  ggplot(aes(x = VILLANI_BPDCN_UP_Score, y = bpdcn_sign_score, color = bm_involvement)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = group_colors) +
  theme_bw() +
  annotate("text", x = -0.5, y = 1, color = "black",
           label = paste0("r = ", round(cor(pdcs$VILLANI_BPDCN_UP_Score, pdcs$bpdcn_sign_score), 2))) +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.text = element_text(color = "black"))
)
dev.off()

# Save
write.table(bpdcn_sign, file = "bpdcn_sign.txt", quote = F, sep = "\t", row.names = F, col.names = F)






# Differentially expressed genes - alternative approach -------------------------------------------

# The following yields very similar results; this fold change ranked gene list can be usedfor GSEA, but
# we ended up just using the section above for the paper

# Extract gene expression matrices
normal_mat <- GetAssayData(normal_pdcs_subset, slot = "data")
malignant_mat <- GetAssayData(malignant_pdcs_subset, slot = "data")

# Make table with gene info
diffgenes.tib <- tibble(gene = rownames(normal_mat),
                        normal_pdcs = rowMeans(normal_mat),
                        malignant_pdcs = rowMeans(malignant_mat))
diffgenes.tib <- mutate(diffgenes.tib, logFC = malignant_pdcs - normal_pdcs)

# Add P-values
diffgenes.tib$pval <- p.adjust(sapply(diffgenes.tib$gene, function(z)
  wilcox.test(x = normal_mat[z,], y = malignant_mat[z,])$p.value), method = "fdr")

# Order, check & save
diffgenes.tib <- arrange(diffgenes.tib, desc(logFC))

pdf("11.2.3_Volcano_alternative.pdf")
diffgenes.tib %>%
  ggplot(aes(x = logFC, y = -log10(pval), color = gene %in% bpdcn_sign)) +
  geom_point() +
  geom_text_repel(data = filter(diffgenes.tib, gene %in% c(bpdcn_sign, tail(diffgenes.tib, 20)$gene)),
                  aes(x = logFC, y = -log10(pval), label = gene),
                  color = "#006400", size = 3, max.overlaps = 30) +
  scale_color_manual(values = c("#708090", "#9acd32")) +
  coord_cartesian(xlim = c(-max(abs(diffgenes.tib$logFC)), max(abs(diffgenes.tib$logFC)))) +
  theme_bw() +
  theme(aspect.ratio = 1, axis.line = element_line(color = "black"), axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"), panel.grid = element_blank())
dev.off()

# Save
#write.table(diffgenes.tib, file = "DiffGenes.txt", quote = F, sep = "\t", row.names = F)
#write.table(dplyr::select(diffgenes.tib, gene, logFC), file = "DiffGenes.rnk", quote = F, sep = "\t", row.names = F, col.names = F)
#bpdcn_sign <- diffgenes.tib %>% filter(pval < 1E-30, logFC > log(2)) %>% .$gene
#write.table(bpdcn_sign, file = "bpdcn_sign.txt", quote = F, sep = "\t", row.names = F, col.names = F)


# Run GSEA (outside of R) -------------------------------------------------------------------------

# Use the gmt file "201204_CombinedSigns.gmt" and the rnk file DiffGenes.rnk
# Pre-Ranked GSEA was run using default settings
# This shows significant enrichment of BPDCN tumor signatures, and also mitochondrial / OXPHOS signatures






# Files for Volker's analysis, 220708 -------------------------------------------------------------
# As of 221020, I don't remember what this is for
#expr.mat <- matrix(nrow = nrow(seu_ls$BM), ncol=length(levels(seu_ls$BM$CellType)), dimnames = list(rownames(seu_ls$BM), levels(seu_ls$BM$CellType)))
#for (n in levels(seu_ls$BM$CellType)) {
#  expr.mat[,n] <- rowMeans( GetAssayData(subset(seu_ls$BM, CellType == n), slot = "data") )
#}
#expr.mat[c("HBD", "CD3D", "CD14"),]
#write.table(expr.mat, file = "Mean_BM_expr.txt", quote = F, sep = "\t")
#seu_ls$BM$CellType %>% table %>% data.frame %>% write.table(file = "CellNumbers.txt", sep = "\t", quote = F, row.names = F, col.names = F)



