# Peter van Galen, 220722
# Relate gene expression in pDCs and malignant BPDCN cells to their mutational status

library(tidyverse)
library(Seurat)
library(readxl)
library(viridis)
library(ComplexHeatmap)
library(circlize)
library(ggforce)

rm(list=ls())

# Set working directory
setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/09_pDC_BPDCN")

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")
donor_colors <- popcol.tib$hex[23:40]
names(donor_colors) <- popcol.tib$pop[23:40]
group_colors <- popcol.tib$hex[41:43]
names(group_colors) <- popcol.tib$pop[41:43]
mut_colors <- popcol.tib$hex[44:46]
names(mut_colors) <- popcol.tib$pop[44:46]

# Load all data
seurat_files <- list.files("../04_XV-seq", pattern = "*.rds", full.names = T)
seu_ls <- lapply(seurat_files, function(x) readRDS(x))
names(seu_ls) <- cutf(basename(seurat_files), d = "_")
seu <- merge(seu_ls[[1]], seu_ls[2:length(seu_ls)])

# Score BPDCN signature that was made in the previous script
bpdcn_sign <- read.table("bpdcn_sign.txt")[,1]
seu <- AddModuleScore(seu, features = list(bpdcn_sign), name = "bpdcn_sign_score")
colnames(seu@meta.data) <- gsub("score1$", "score", colnames(seu@meta.data))

# Load pDC gene expression and metadata (including UMAP coordinates) from 9.1_pDC_object.R
pdcs <- readRDS("pdcs.rds")
# Add signature score to pdcs
pdcs$bpdcn_sign_score <- seu$bpdcn_sign_score[colnames(pdcs)]

# Load genotyping information
genotyping_tables.tib <- read_excel("../04_XV-seq/XV-seq_overview.xlsx")
# Replace different MTAP entries with one, just as in 4.1_Add_GoT-XV_to_Seurat.R
genotyping_tables.tib$Mutation <- gsub("MTAP.rearr.*", "MTAP.rearr", genotyping_tables.tib$Mutation)


# UMAP of bpdcn_sign_score in pDCs ----------------------------------------------------------------
pdf("9.3.1_ModuleScores.pdf", width = 6, height = 6)
FeaturePlot(pdcs, features = "bpdcn_sign_score", pt.size = 0.6, order = T) +
  scale_color_viridis() +
  theme_void() +
  theme(aspect.ratio = 1, panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(hjust = 0.5))
dev.off()


# Make violin plots with mutation calls vs. BPDCN sign scores -------------------------------------

# Extract metadata
metadata_tib <- as_tibble(pdcs@meta.data, rownames = "cell")

# Sina plot of scores in all donors without involvement - color by progression mutations
pdf("9.3.2_Non-involved_pDC_scores.pdf")
metadata_tib %>% filter(bm_involvement != "Yes") %>%
  ggplot(aes(x = orig.ident, y = bpdcn_sign_score)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  geom_sina(color = gsub("wildtype", mut_colors[1], gsub("mutant", mut_colors[2], gsub("no call", mut_colors[3],
    filter(arrange(metadata_tib, orig.ident), bm_involvement != "Yes")$progression)))) +
  theme_bw() +
  theme(panel.grid = element_blank())
dev.off()

# Wrangle
Pt10_Muts <- filter(genotyping_tables.tib, Sample == "Pt10Dx") %>% arrange(`Founder or progression mutation`) %>% .$Mutation
plot_tib <- metadata_tib %>% dplyr::select(orig.ident, all_of(Pt10_Muts), bpdcn_sign_score) %>% filter(orig.ident == "Pt10Dx") %>%
  pivot_longer(cols = all_of(Pt10_Muts), names_to = "Mutation", values_to = "XVseq") %>%
  mutate(XVseq = factor(XVseq, levels = c("no call", "wildtype", "mutant")),
         Mutation = factor(Mutation, levels = unique(Pt10_Muts)))
plot_tib <- plot_tib %>% arrange(Mutation, XVseq) %>%
  mutate(XVseq = gsub("wildtype", "#32cd32", gsub("mutant", "#dc143c", gsub("no call", "#dcdcdc", as.character(XVseq)))))

# Visualize Pt10Dx sina plots with mut/wt calls for every mutation
pdf("9.3.3_Pt10Dx_pDC_scores_XVseq.pdf", width = 10, height = 6)
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
pdf("9.3.4_Pt10Dx_pDC_scores_combined_XVseq.pdf", width = 3, height = 4)
plot2_tib %>%
  ggplot(aes(x = Mutation, y = bpdcn_sign_score)) +
  geom_sina(color = plot2_tib$XVseq) +
  geom_violin(fill = NA) +
  theme_bw() +
  theme(aspect.ratio = 2, panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"), axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


# Subset data to make two heatmaps ----------------------------------------------------------------

# First, subset for a heatmap showing pDCs from healthy donors and those with malignant infiltration
healthy_ids <- subset(pdcs, subset = bm_involvement == "HD") %>% colnames() %>% sample(6*30)
malignant_ids <- as_tibble(pdcs@meta.data, rownames = "cell") %>%
  filter(bm_involvement == "Yes", any_mut != "no call") %>%
  group_by(orig.ident) %>% slice_sample(n = 30) %>% .$cell
pdcs_healthy_malignant <- subset(pdcs, cells = c(healthy_ids, malignant_ids))

# Then, subset data for for a heatmap showing pDCs of skin-only patients
pdcs_skin_only <- subset(pdcs, bm_involvement == "No")
#pdcs_skin_only$orig.ident <- factor(seu_skin_only$orig.ident, levels = c("Pt1Rem", "Pt5Dx", "Pt9Dx", "Pt10Dx", "Pt12Dx"))
keep_cells_left <- as_tibble(pdcs_skin_only@meta.data, rownames = "cell") %>%
  filter(bpdcn_sign_score < 0.5, any_mut != "no call") %>% arrange(orig.ident) %>% .$cell
keep_cells_right <- as_tibble(pdcs_skin_only@meta.data, rownames = "cell") %>%
  filter(bpdcn_sign_score > 0.5) %>% arrange(orig.ident) %>% .$cell
pdcs_skin_subset <- subset(pdcs_skin_only, cells = c(keep_cells_left, keep_cells_right))


# What mutations to show (for both heatmaps)?
mutation_detection <- as_tibble(merge(pdcs_healthy_malignant, pdcs_skin_subset)@meta.data) %>%
  dplyr::select(unique(na.omit(genotyping_tables.tib$Mutation))) %>%
  apply(., 2, function(x) sum(grepl("mutant|wildtype", x)))
show_f_mut <- intersect(names(mutation_detection)[mutation_detection > 3],
                        filter(genotyping_tables.tib, `Founder or progression mutation` == "Founder")$Mutation)
show_p_mut <- intersect(names(mutation_detection)[mutation_detection > 3],
                        filter(genotyping_tables.tib, `Founder or progression mutation` == "Progression")$Mutation)

# Load gene expression data, then normalize (like in 8_Heatmaps.R)
expr_mat <- as.matrix(GetAssayData(pdcs_healthy_malignant, slot = "data"))[bpdcn_sign,]
expr_mat <- expr_mat - rowMeans(expr_mat)
z.lim <- c(-2, 4)
expr_mat[expr_mat < z.lim[1]] <- z.lim[1]
expr_mat[expr_mat > z.lim[2]] <- z.lim[2]

# Adjust some colors
mut_colors[3] <- "white"
col_bw <- colorRamp2(breaks = c(min(pdcs$bpdcn_sign_score), max(pdcs$bpdcn_sign_score)), colors = c("white", "black"))

# Define annotation objects
hm_anno_df <- pdcs_healthy_malignant@meta.data[,c("orig.ident", "bpdcn_sign_score", show_f_mut, show_p_mut)]
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
expr_mat <- as.matrix(GetAssayData(pdcs_skin_subset, slot = "data"))[bpdcn_sign,]
expr_mat <- expr_mat - rowMeans(expr_mat)
z.lim <- c(-2, 4)
expr_mat[expr_mat < z.lim[1]] <- z.lim[1]
expr_mat[expr_mat > z.lim[2]] <- z.lim[2]

# Define annotation objects
hm_anno_df <- pdcs_skin_subset@meta.data[,c("orig.ident", "bpdcn_sign_score", show_f_mut, show_p_mut)]
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
              column_split = factor(ifelse(hm_anno_df$bpdcn_sign_score > 0.5, yes = "Putative malignant", no = "Normal"),
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


# Simplified skin-only pDC heatmap for Keynote ----------------------------------------------------

# Collapse mut/wt calls
hm_anno_df$founder_summary <- ifelse(apply(hm_anno_df[,show_f_mut], 1, function(x) sum(grepl("wildtype", x) > 0)), yes = "wildtype", no = "no call")
hm_anno_df$founder_summary <- ifelse(apply(hm_anno_df[,show_f_mut], 1, function(x) sum(grepl("mutant", x) > 0)), yes = "mutant", no = hm_anno_df$founder_summary)
hm_anno_df$prog_summary <- ifelse(apply(hm_anno_df[,show_p_mut], 1, function(x) sum(grepl("wildtype", x) > 0)), yes = "wildtype", no = "no call")
hm_anno_df$prog_summary <- ifelse(apply(hm_anno_df[,show_p_mut], 1, function(x) sum(grepl("mutant", x) > 0)), yes = "mutant", no = hm_anno_df$prog_summary)

# Define annotation objects
hm_anno_df <- pdcs_skin_subset@meta.data[,c("orig.ident", "bpdcn_sign_score", show_f_mut, show_p_mut)]
top_anno.ha <- HeatmapAnnotation(Donor = as.character(hm_anno_df$orig.ident),
                                 Score = hm_anno_df$bpdcn_sign_score,
                                 col = list(Donor = donor_colors, Score = col_bw),
                                 annotation_name_gp = gpar(fontsize = 10),
                                 border = T, simple_anno_size = unit(3, "mm"))
bottom_anno.ha <- HeatmapAnnotation(Founder = as.matrix(hm_anno_df$founder_summary),
                                    Prog = as.matrix(hm_anno_df$prog_summary),
                                    col = list(Founder = mut_colors, Prog = mut_colors),
                                    annotation_name_gp = gpar(fontsize = 6),
                                    border = T, na_col = "white", simple_anno_size = unit(3, "mm"))

# Create Heatmap object
hm <- Heatmap(as.matrix(expr_mat),
              col = colItay(c(1:11))[3:11],
              cluster_rows = F,
              cluster_columns = F,
              row_names_gp = gpar(fontsize = 6),
              show_column_names = F,
              column_split = factor(ifelse(hm_anno_df$bpdcn_sign_score > 0.5, yes = "Putative malignant", no = "Normal"),
                                    levels = c("Normal", "Putative malignant")),
              top_annotation = top_anno.ha,
              bottom_annotation = bottom_anno.ha,
              name = "Expr",
              column_title = "Cells classified as pDC in samples\nwithout bone marrow involvement",
              column_title_gp = gpar(fontsize = 10),
              border = T,
              use_raster = T,
              raster_quality = 10)

pdf("9.3.7_Skin_only_pDC_heatmap_simple.pdf", width = 5, height = 6)
print(hm)
dev.off()



