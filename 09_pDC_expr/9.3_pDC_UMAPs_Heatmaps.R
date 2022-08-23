# Peter van Galen, 220722
# Relate gene expression in pDCs and malignant BPDCN cells to their mutational status

library(tidyverse)
library(Seurat)
library(readxl)
library(ggforce)
library(gridExtra)
library(cowplot)
library(viridis)
library(ComplexHeatmap)
library(limma)
library(data.table)
library(ggrepel)
library(circlize)

rm(list=ls())

# Set working directory
setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/09_pDC_Expr")

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")
donor_colors <- popcol.tib$hex[24:41]
names(donor_colors) <- popcol.tib$pop[24:41]
group_colors <- popcol.tib$hex[42:44]
names(group_colors) <- popcol.tib$pop[42:44]
mut_colors <- popcol.tib$hex[45:47]
names(mut_colors) <- popcol.tib$pop[45:47]

# Load pDC gene expression and metadata and signature
pdcs <- readRDS("pdcs.rds")
bpdcn_sign <- read.table("bpdcn_sign.txt")[,1]

# Load genotyping information
genotyping_tables.tib <- read_excel("../04_XV-seq/XV-seq_overview.xlsx")
# Replace different MTAP entries with one, just as in 4.1_Add_GoT-XV_to_Seurat.R
genotyping_tables.tib$Mutation <- gsub("MTAP.rearr.*", "MTAP.rearr", genotyping_tables.tib$Mutation)


# Make violin plots with mutation calls vs. BPDCN sign scores -------------------------------------

# Add signatures score
pdcs <- AddModuleScore(pdcs, features = list(bpdcn_sign), name = "bpdcn_sign_score")
colnames(pdcs@meta.data) <- gsub("score1$", "score", colnames(pdcs@meta.data))

# Extract metadata, look at BPDCN signature score vs. sample vs. mutations
metadata_tib <- as_tibble(pdcs@meta.data, rownames = "cell")
pdf("9.3.1_Non-involved_pDC_scores.pdf")
metadata_tib %>% filter(bm_involvement != "Yes") %>%
  ggplot(aes(x = orig.ident, y = bpdcn_sign_score)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_sina(color = gsub("wildtype", mut_colors[1], gsub("mutant", mut_colors[2], gsub("no call", mut_colors[3],
    filter(arrange(metadata_tib, orig.ident), bm_involvement != "Yes")$progression)))) +
  theme_bw() +
  theme(panel.grid = element_blank())
dev.off()

# Wrangle
Pt10_Muts <- filter(genotyping_tables.tib, Sample == "Pt10Dx") %>% arrange(`Founder or progression mutation`) %>% .$Mutation
plot_tib <- metadata_tib %>% dplyr::select(orig.ident, all_of(Pt10_Muts), bpdcn_sign_score) %>% filter(orig.ident == "Pt10Dx") %>%
  pivot_longer(cols = Pt10_Muts, names_to = "Mutation", values_to = "XVseq") %>%
  mutate(XVseq = factor(XVseq, levels = c("no call", "wildtype", "mutant")),
         Mutation = factor(Mutation, levels = unique(Pt10_Muts)))
plot_tib <- plot_tib %>% arrange(Mutation, XVseq) %>%
  mutate(XVseq = gsub("wildtype", "#32cd32", gsub("mutant", "#dc143c", gsub("no call", "#dcdcdc", as.character(XVseq)))))

# Visualize violin plots with mut/wt calls for every mutation
pdf("9.3.2_Pt10Dx_pDC_scores_XVseq.pdf", width = 10, height = 6)
plot_tib %>%
  ggplot(aes(x = Mutation, y = bpdcn_sign_score)) +
  geom_sina(color = plot_tib$XVseq) +
  geom_violin(fill = NA) +
  theme_bw() +
  theme(aspect.ratio = 0.6, panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"), axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# Wrangle to combine collapse violins from the previous visualization
plot2_tib <- metadata_tib %>% dplyr::select(orig.ident, founder, progression, bpdcn_sign_score) %>%
  filter(orig.ident == "Pt10Dx") %>%
  pivot_longer(cols = c("founder", "progression"), names_to = "Mutation", values_to = "XVseq") %>%
  mutate(XVseq = factor(XVseq, levels = c("no call", "wildtype", "mutant")))
plot2_tib <- plot2_tib %>% arrange(Mutation, XVseq) %>%
  mutate(XVseq = gsub("wildtype", "#32cd32", gsub("mutant", "#dc143c", gsub("no call", "#dcdcdc", as.character(XVseq)))))

# Visualize
pdf("9.3.3_Pt10Dx_pDC_scores_combined_XVseq.pdf", width = 3, height = 4)
plot2_tib %>%
  ggplot(aes(x = Mutation, y = bpdcn_sign_score)) +
  geom_sina(color = plot2_tib$XVseq) +
  geom_violin(fill = NA) +
  theme_bw() +
  theme(aspect.ratio = 2, panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"), axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


# Skin-only pDC UMAPs -----------------------------------------------------------------------------

# Subset pDC object for skin-only patients
seu_skin_only <- subset(pdcs, bm_involvement == "No")
seu_skin_only$orig.ident <- factor(seu_skin_only$orig.ident, levels = c("Pt1Rem", "Pt5Dx", "Pt9Dx", "Pt10Dx", "Pt12Dx"))

# Calculate UMAP coordinates for pDCs using malignant BPDCN cell signature genes
seu_skin_only <- NormalizeData(seu_skin_only)
seu_skin_only <- ScaleData(seu_skin_only, features = rownames(seu_skin_only))
seu_skin_only <- RunPCA(seu_skin_only, features = bpdcn_sign)
seu_skin_only <- RunUMAP(seu_skin_only, reduction = "pca", dims = 1:10)

# Plot
p1 <- DimPlot(seu_skin_only, group.by = "orig.ident", cols = donor_colors[unique(as.character(seu_skin_only$orig.ident))]) +
  theme_void() + theme(aspect.ratio = 1, panel.border = element_rect(color = "black", fill = NA))
p2 <- DimPlot(seu_skin_only, group.by = "progression", cols = mut_colors, order = "call") + theme(aspect.ratio = 1) +
  theme_void() + theme(aspect.ratio = 1, panel.border = element_rect(color = "black", fill = NA))
p3 <- FeaturePlot(seu_skin_only, features = "VILLANI_BPDCN_UP_Score") + scale_color_viridis() + theme(aspect.ratio = 1) +
  theme_void() + theme(aspect.ratio = 1, panel.border = element_rect(color = "black", fill = NA))
p4 <- FeaturePlot(seu_skin_only, features = "bpdcn_sign_score") + scale_color_viridis() + theme(aspect.ratio = 1) +
  theme_void() + theme(aspect.ratio = 1, panel.border = element_rect(color = "black", fill = NA))

pdf("9.3.4_Skin_only_pDCs_UMAP.pdf", width = 7, height = 7)
plot_grid(p1, p2, p3, p4, align = "v", rel_heights = c(1,1,1,1))
dev.off()


# Subset data to make two heatmaps ----------------------------------------------------------------

# First, subset for a heatmap showing pDCs from healthy donors and those with malignant infiltration
healthy_ids <- subset(pdcs, subset = bm_involvement == "HD") %>% colnames() %>% sample(6*30)
malignant_ids <- as_tibble(pdcs@meta.data, rownames = "cell") %>%
  filter(bm_involvement == "Yes", any_mut != "no call") %>%
  group_by(orig.ident) %>% slice_sample(n = 30) %>% .$cell
seu_healthy_malignant <- subset(pdcs, cells = c(healthy_ids, malignant_ids))

# Then, subset data for for a heatmap showing pDCs of skin-only patients
remove_cells <- as_tibble(seu_skin_only@meta.data, rownames = "cell") %>%
  filter(bpdcn_sign_score < 0, any_mut == "no call") %>% .$cell
keep_cells <- as_tibble(seu_skin_only@meta.data, rownames = "cell") %>%
  filter(! cell %in% remove_cells) %>% arrange(orig.ident) %>% .$cell
seu_skin_subset <- subset(seu_skin_only, cells = keep_cells)

# What mutations to show (for both heatmaps)?
mutation_detection <- as_tibble(merge(seu_healthy_malignant, seu_skin_subset)@meta.data) %>%
  dplyr::select(unique(genotyping_tables.tib$Mutation)) %>%
  apply(., 2, function(x) sum(grepl("mutant|wildtype", x)))
show_f_mut <- intersect(names(mutation_detection)[mutation_detection > 3],
                        filter(genotyping_tables.tib, `Founder or progression mutation` == "Founder")$Mutation)
show_p_mut <- intersect(names(mutation_detection)[mutation_detection > 3],
                        filter(genotyping_tables.tib, `Founder or progression mutation` == "Progression")$Mutation)

# Load gene expression data, then normalize (like in 8_Heatmaps.R)
expr_mat <- as.matrix(GetAssayData(seu_healthy_malignant, slot = "data"))[bpdcn_sign,]
expr_mat <- expr_mat - rowMeans(expr_mat)
z.lim <- c(-2, 4)
expr_mat[expr_mat < z.lim[1]] <- z.lim[1]
expr_mat[expr_mat > z.lim[2]] <- z.lim[2]

# Adjust some colors
mut_colors[3] <- "white"
col_bw <- colorRamp2(breaks = c(min(pdcs$bpdcn_sign_score), max(pdcs$bpdcn_sign_score)), colors = c("white", "black"))

# Define annotation objects
hm_anno_df <- seu_healthy_malignant@meta.data[,c("orig.ident", "bpdcn_sign_score", show_f_mut, show_p_mut)]
top_anno.ha <- HeatmapAnnotation(Donor = as.character(hm_anno_df$orig.ident),
                                 Score = hm_anno_df$bpdcn_sign_score,
                                 col = list(Donor = donor_colors, Score = col_bw),
                                 annotation_name_gp = gpar(fontsize = 10),
                                 border = T, simple_anno_size = unit(3, "mm"))
bottom_anno.ha <- HeatmapAnnotation(Founder = as.matrix(hm_anno_df[,show_f_mut]),
                                    Prog = as.matrix(hm_anno_df[,show_p_mut]),
                                    col = list(Founder = mut_colors, Prog = mut_colors),
                                    annotation_name_gp = gpar(fontsize = 6),
                                    border = T, na_col = "white", simple_anno_size = unit(1.8, "mm"))

# Create Heatmap object
hm <- Heatmap(as.matrix(expr_mat),
              col = colItay(c(1:11))[3:11],
              cluster_rows = F,
              cluster_columns = F,
              row_names_gp = gpar(fontsize = 6),
              show_column_names = F,
              column_split = factor(ifelse(grepl("BM", hm_anno_df$orig.ident), yes = "BM", no = "BPDCN"), levels = c("BM", "BPDCN")),
              top_annotation = top_anno.ha,
              bottom_annotation = bottom_anno.ha,
              name = "Expr",
              column_title = "Cells classified as pDC in healthy BM\nand samples with BM involvement",
              column_title_gp = gpar(fontsize = 10),
              border = T,
              use_raster = T,
              raster_quality = 10)

pdf("9.3.5_Healthy_malignant_pDC_heatmap.pdf", width = 5, height = 6)
print(hm)
dev.off()


# Skin-only pDC heatmap ---------------------------------------------------------------------------

# Load gene expression data, then normalize (like in 8_Heatmaps.R)
expr_mat <- as.matrix(GetAssayData(seu_skin_subset, slot = "data"))[bpdcn_sign,]
expr_mat <- expr_mat - rowMeans(expr_mat)
z.lim <- c(-2, 4)
expr_mat[expr_mat < z.lim[1]] <- z.lim[1]
expr_mat[expr_mat > z.lim[2]] <- z.lim[2]

# Define annotation objects
hm_anno_df <- seu_skin_subset@meta.data[,c("orig.ident", "bpdcn_sign_score", show_f_mut, show_p_mut)]
top_anno.ha <- HeatmapAnnotation(Donor = as.character(hm_anno_df$orig.ident),
                                 Score = hm_anno_df$bpdcn_sign_score,
                                 col = list(Donor = donor_colors, Score = col_bw),
                                 annotation_name_gp = gpar(fontsize = 10),
                                 border = T, simple_anno_size = unit(3, "mm"))
bottom_anno.ha <- HeatmapAnnotation(Founder = as.matrix(hm_anno_df[,show_f_mut]),
                                    Prog = as.matrix(hm_anno_df[,show_p_mut]),
                                    col = list(Founder = mut_colors, Prog = mut_colors),
                                    annotation_name_gp = gpar(fontsize = 6),
                                    border = T, na_col = "white", simple_anno_size = unit(1.8, "mm"))

# Create Heatmap object
hm <- Heatmap(as.matrix(expr_mat),
              col = colItay(c(1:11))[3:11],
              cluster_rows = F,
              cluster_columns = F,
              row_names_gp = gpar(fontsize = 6),
              show_column_names = F,
              column_split = factor(ifelse(hm_anno_df$bpdcn_sign_score > 0, yes = "Putative malignant", no = "Normal"),
                                    levels = c("Normal", "Putative malignant")),
              top_annotation = top_anno.ha,
              bottom_annotation = bottom_anno.ha,
              name = "Expr",
              column_title = "Cells classified as pDC in samples\nwithout bone marrow involvement",
              column_title_gp = gpar(fontsize = 10),
              border = T,
              use_raster = T,
              raster_quality = 10)

pdf("9.3.6_Skin_only_pDC_heatmap.pdf", width = 5, height = 6)
print(hm)
dev.off()


