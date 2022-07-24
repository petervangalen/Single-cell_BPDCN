# Peter van Galen, 220722
# Relate gene expression in pDCs and malignant BPDCN cells to their mutational status

library(tidyverse)
library(Seurat)
library(readxl)
library(ggforce)
library(gridExtra)
library(viridis)
library(ComplexHeatmap)
library(limma)
library(data.table)
library(ggrepel)
library(circlize)

rm(list=ls())

# Set working directory
setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/Github/9_pDC_Expr")

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")
donor_colors <- popcol.tib$hex[24:41]
names(donor_colors) <- popcol.tib$pop[24:41]
group_colors <- popcol.tib$hex[42:44]
names(group_colors) <- popcol.tib$pop[42:44]
mut_colors <- popcol.tib$hex[45:47]
names(mut_colors) <- popcol.tib$pop[45:47]

# Or pDCs alone
pdcs <- readRDS("pDCs.rds")

# Load genotyping information
genotyping_tables.tib <- read_excel("../4_XV-seq/XV-seq_overview.xlsx")
# Replace different MTAP entries with one, just as in 4.1_Add_GoT-XV_to_Seurat.R
genotyping_tables.tib$Mutation <- gsub("MTAP.rearr.*", "MTAP.rearr", genotyping_tables.tib$Mutation)


# Add mutation calls to Seurat object -------------------------------------------------------------

# Make vectors of different mutation types
a_mut <- genotyping_tables.tib$Mutation %>% unique
f_mut <- filter(genotyping_tables.tib, `Founder or progression mutation` == "Founder") %>% .$Mutation %>% unique
p_mut <- filter(genotyping_tables.tib, `Founder or progression mutation` == "Progression") %>% .$Mutation %>% unique

# Add columns with mutation calls to metadata
pdcs$any_mut <- ifelse(apply(pdcs@meta.data[,a_mut], 1, function(x) sum(grepl("wildtype", x) > 0)), yes = "wildtype", no = "no call")
pdcs$any_mut <- ifelse(apply(pdcs@meta.data[,a_mut], 1, function(x) sum(grepl("mutant", x) > 0)), yes = "mutant", no = pdcs$any_mut)
pdcs$founder <- ifelse(apply(pdcs@meta.data[,f_mut], 1, function(x) sum(grepl("wildtype", x) > 0)), yes = "wildtype", no = "no call")
pdcs$founder <- ifelse(apply(pdcs@meta.data[,f_mut], 1, function(x) sum(grepl("mutant", x) > 0)), yes = "mutant", no = pdcs$founder)
pdcs$progression <- ifelse(apply(pdcs@meta.data[,p_mut], 1, function(x) sum(grepl("wildtype", x) > 0)), yes = "wildtype", no = "no call")
pdcs$progression <- ifelse(apply(pdcs@meta.data[,p_mut], 1, function(x) sum(grepl("mutant", x) > 0)), yes = "mutant", no = pdcs$progression)


# Generate malignant BPDCN gene signature ---------------------------------------------------------

# Identify genes that correlate with the Villani score
#expr_mat <- GetAssayData(pdcs, slot = "data")
#villani_score <- pdcs$VILLANI_BPDCN_UP_Score
#options(warn=0)
#my_cor <- apply(expr_mat, 1, function(x) cor(x, villani_score))
#my_cor <- sort(my_cor, decreasing = T) # BCL2 is #189
#plot(my_cor, pch = ".")
#abline(v = 50, col = "red")

# Select different signatures; eventually, I settled on ...
#bpdcn_sign <- names(my_cor[1:20]); pdf_name <- "_Cor_Top20"
#bpdcn_sign <- names(my_cor[1:30]); pdf_name <- "_Cor_Top30" # better than villani
#bpdcn_sign <- names(my_cor[1:40]); pdf_name <- "_Cor_Top40" # better than villani
#bpdcn_sign <- names(my_cor[1:50]); pdf_name <- "_Cor_Top50" # better than villani
#bpdcn_sign <- names(my_cor[1:60]); pdf_name <- "_Cor_Top60"
#bpdcn_sign <- names(my_cor[1:70]); pdf_name <- "_Cor_Top70"

# Or just use my data (don't correlate with Villani)
#pdcs_expr <- rowMeans(GetAssayData(subset(pdcs, donor_group != "bm_involvement" & progression != "mutant")))
#bpdcn_expr <- rowMeans(GetAssayData(subset(pdcs, donor_group == "bm_involvement")))
#my_sig <- sort(bpdcn_expr - pdcs_expr, decreasing = T) # BCL2 is #34
#bpdcn_sign <- names( my_sig[my_sig > log(2)] ); pdf_name <- "_MySig1"

# Using Chloe Villani's existing signature
bpdcn_sign <- intersect(read_excel("../9_pDC_Expr/Signatures.xlsx")[-1,]$VILLANI_BPDCN_UP, rownames(pdcs)); pdf_name <- "_Villani"

# Or define my own signature (to the end of this section)
normal_ids <- colnames(subset(pdcs, donor_group != "bm_involvement" & pdcs$progression != "mutant"))
malignant_ids <- colnames(subset(pdcs, donor_group == "bm_involvement"))
markerGenes <- FindMarkers(pdcs, ident.1 = malignant_ids, ident.2 = normal_ids)
markerGenes_tib <- as_tibble(markerGenes, rownames = "gene") %>% arrange(-avg_log2FC)
bpdcn_sign <- markerGenes_tib %>% filter(p_val_adj < 1E-100, avg_log2FC > 1) %>% .$gene
highlight_genes <- c("IGLL1", "HES6", "PDLIM1", "SERPINB2", "ARPP21", "TCL1A", "LAMP5", "FCER1A", "MYBPH", "GSTP1",
                     "CLEC11A", "MS4A4E", "CD7", "S100A4", "SOX4", "HLA-DQA2", "NCAM1", "BCL2", "ASXL1")
pdf_name <- "_MySigFinal"
# Volcano plot
pdf("9.2.0_MySigFinal_Volcano.pdf")
markerGenes_tib %>%
  ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj), color = gene %in% bpdcn_sign)) +
  geom_point() +
  geom_text_repel(data = filter(markerGenes_tib, gene %in% highlight_genes),
                  aes(x = avg_log2FC, y = -log10(p_val_adj), label = gene),
                  color = "#006400", size = 3, max.overlaps = 15) +
  scale_color_manual(values = c("#708090", "#9acd32")) +
  theme_bw() +
  theme(aspect.ratio = 1, axis.line = element_line(color = "black"), axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"))
dev.off()

# How does this compare to Chloe Villani's signature?
pdcs$bpdcn_sign_score <- NULL
pdcs <- AddModuleScore(pdcs, features = list(bpdcn_sign), name = "bpdcn_sign_score")
colnames(pdcs@meta.data) <- gsub("score1$", "score", colnames(pdcs@meta.data))
pdcs@meta.data %>% 
  ggplot(aes(x = VILLANI_BPDCN_UP_Score, y = bpdcn_sign_score)) +
  geom_point() +
  theme(aspect.ratio = 1)


# Make violin plots with mutation calls vs. BPDCN sign scores -------------------------------------

# What does the data look like?
p1 <- DimPlot(pdcs, group.by = "orig.ident", cols = donor_colors[unique(pdcs$orig.ident)]) + theme(aspect.ratio = 1)
p2 <- DimPlot(pdcs, group.by = "donor_group", cols = group_colors) + theme(aspect.ratio = 1)
p3 <- FeaturePlot(pdcs, features = "bpdcn_sign_score") + scale_color_viridis() + theme(aspect.ratio = 1)
grid.arrange(p1, p2, p3, ncol = 2)

# Extract metadata, look at BPDCN signature score vs. sample vs. mutations
metadata_tib <- as_tibble(pdcs@meta.data, rownames = "cell") %>% mutate(UMAP1 = pdcs@reductions$umap@cell.embeddings[,1],
                                                                        UMAP2 = pdcs@reductions$umap@cell.embeddings[,2])
pdf("9.2.0_Non-involved_pDC_scores.pdf")
metadata_tib %>% filter(donor_group != "bm_involvement") %>% 
  ggplot(aes(x = orig.ident, y = bpdcn_sign_score)) +
  geom_hline(yintercept = 0) +
  geom_sina(color = gsub("wildtype", "#32cd32", gsub("mutant", "#dc143c", gsub("no call", "#dcdcdc",
    filter(metadata_tib, donor_group != "bm_involvement")$progression)))) +
  theme_bw() +
  theme(panel.grid = element_blank())
dev.off()

# Of the cells with a malignant BPDCN expression profile, mutations were only captured in Pt10Dx
metadata_tib %>% filter(donor_group != "bm_involvement", bpdcn_sign_score > 0.4) %>%
  filter(if_any(everything(), ~ grepl("mutant", .))) %>% .$orig.ident

# Wrangle
Pt10_Muts <- filter(genotyping_tables.tib, Sample == "Pt10Dx") %>% arrange(`Founder or progression mutation`) %>% .$Mutation
plot_tib <- metadata_tib %>% dplyr::select(orig.ident, Pt10_Muts, bpdcn_sign_score) %>% filter(orig.ident == "Pt10Dx") %>%
  pivot_longer(cols = Pt10_Muts, names_to = "Mutation", values_to = "XVseq") %>%
  mutate(XVseq = factor(XVseq, levels = c("no call", "wildtype", "mutant")),
         Mutation = factor(Mutation, levels = unique(Pt10_Muts)))
plot_tib <- plot_tib %>% arrange(Mutation, XVseq) %>%
  mutate(XVseq = gsub("wildtype", "#32cd32", gsub("mutant", "#dc143c", gsub("no call", "#dcdcdc", as.character(XVseq)))))

# Visualize violin plots with mut/wt calls for every mutation
pdf(paste0("9.2.1_Pt10Dx_pDCs_XVseq", pdf_name, ".pdf"), width = 10, height = 6)
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
pdf(paste0("9.2.1_Pt10Dx_pDCs_XVseq2", pdf_name, ".pdf"), width = 3, height = 4)
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
seu_skin_only <- subset(pdcs, donor_group == "skin_only")
seu_skin_only$orig.ident <- factor(seu_skin_only$orig.ident, levels = c("Pt1Mrd", "Pt5Dx", "Pt9Dx", "Pt10Dx", "Pt12Dx"))

# Calculate UMAP coordinates for pDCs using malignant BPDCN cell signature genes
seu_skin_only <- NormalizeData(seu_skin_only)
seu_skin_only <- ScaleData(seu_skin_only, features = rownames(seu_skin_only))
seu_skin_only <- RunPCA(seu_skin_only, features = bpdcn_sign)
seu_skin_only <- RunUMAP(seu_skin_only, reduction = "pca", dims = 1:10)

# Plot
p1 <- DimPlot(seu_skin_only, group.by = "orig.ident", cols = donor_colors[unique(as.character(seu_skin_only$orig.ident))]) +
  theme(aspect.ratio = 1)
p2 <- DimPlot(seu_skin_only, group.by = "progression", cols = mut_colors, order = "call") + theme(aspect.ratio = 1)

pdf(paste0("9.2.2_Skin_only_pDCs_UMAP", pdf_name, ".pdf"), width = 10, height = 5)
grid.arrange(p1, p2, ncol = 2)
dev.off()

FeaturePlot(seu_skin_only, features = "TCL1A") + scale_color_viridis() + theme(aspect.ratio = 1)
FeaturePlot(seu_skin_only, features = "bpdcn_sign_score") + scale_color_viridis() + theme(aspect.ratio = 1)


# Normal vs. maligannt heatmap --------------------------------------------------------------------

# First, make a heatmap showing pDCs from healthy donors and those with malignant infiltration
healthy_seu <- subset(pdcs, subset = donor_group == "healthy_bm")
malignant_seu <- subset(pdcs, subset = donor_group == "bm_involvement")
malignant_subset <- subset(malignant_seu, cells = sample(colnames(malignant_seu), size = ncol(healthy_seu)))
seu_healthy_malignant <- merge(healthy_seu, malignant_subset)

# Load relevant metadata for heatmap annotation & remove mutations that don't have any mutant cells
hm_anno_df <- seu_healthy_malignant@meta.data[,c("orig.ident", "orig.ident2", "bpdcn_sign_score")]

# Load gene expression data
expr_mat <- as.matrix(GetAssayData(seu_healthy_malignant, slot = "data"))[bpdcn_sign,]
expr_mat[1:3,1:3]

# Normalize (like in 8_Heatmaps.R)
expr_mat <- expr_mat - rowMeans(expr_mat)
z.lim <- c(-2, 4)
expr_mat[expr_mat < z.lim[1]] <- z.lim[1]
expr_mat[expr_mat > z.lim[2]] <- z.lim[2]

# Define annotation objects
col_bw <- colorRamp2(breaks = c(min(hm_anno_df$bpdcn_sign_score), max(hm_anno_df$bpdcn_sign_score)), colors = c("white", "black"))
top_anno.ha <- HeatmapAnnotation(Donor = hm_anno_df$orig.ident,
                                 Donor2 = hm_anno_df$orig.ident2,
                                 Score = hm_anno_df$bpdcn_sign_score,
                                 col = list(Donor = donor_colors, Donor2 = donor_colors, Score = col_bw),
                                 annotation_name_gp = gpar(fontsize = 10),
                                 border = T)

# Create Heatmap object
hm <- Heatmap(as.matrix(expr_mat),
              col = colItay(c(1:11))[3:11],
              cluster_rows = F,
              cluster_columns = T,
              row_names_gp = gpar(fontsize = 6),
              show_column_names = F,
              column_split = 2,
              top_annotation = top_anno.ha,
#              bottom_annotation = bottom_anno.ha,
              name = "Expr",
              column_title = "Cells classified as pDC in healthy donors (BM) and samples with bone marrow involvement",
              column_title_gp = gpar(fontsize = 12),
              border = T,
              use_raster = T,
              raster_quality = 10)

pdf(paste0("9.2.3_Healthy_malignant_pDC_heatmap", pdf_name, ".pdf"), width = 10, height = 5)
print(hm)
dev.off()


# Skin-only pDC heatmap ---------------------------------------------------------------------------

# Remove some cells to improve visualization
seu_skin_only$any_p_mut <- apply(seu_skin_only@meta.data[,p_mut], 1, function(x) sum(grepl("mutant|wildtype", x) > 0))
remove_cells <- as_tibble(seu_skin_only@meta.data, rownames = "cell") %>%
  filter(bpdcn_sign_score < 0, any_p_mut == F) %>% .$cell
remove_cells <- sample(remove_cells, size = length(remove_cells)*4/5)
seu_skin_subset <- subset(seu_skin_only, cells = setdiff(colnames(seu_skin_only), remove_cells))

# Load relevant metadata for heatmap annotation & remove mutations that don't have any mutant cells
hm_anno_df <- seu_skin_subset@meta.data[,c("orig.ident", p_mut, "bpdcn_sign_score")]
detected_mutations <- colnames(hm_anno_df)[apply(hm_anno_df, 2, function(x) sum(grepl("mutant", x))) > 0]

# Load gene expression data
expr_mat <- as.matrix(GetAssayData(seu_skin_subset, slot = "data"))[bpdcn_sign,]
expr_mat[1:3,1:3]

# Normalize (like in 8_Heatmaps.R)
expr_mat <- expr_mat - rowMeans(expr_mat)
z.lim <- c(-2, 4)
expr_mat[expr_mat < z.lim[1]] <- z.lim[1]
expr_mat[expr_mat > z.lim[2]] <- z.lim[2]


# Define annotation objects
col_bw <- colorRamp2(breaks = c(min(hm_anno_df$bpdcn_sign_score), max(hm_anno_df$bpdcn_sign_score)), colors = c("white", "black"))
top_anno.ha <- HeatmapAnnotation(Donor = hm_anno_df$orig.ident,
                                 Score = hm_anno_df$bpdcn_sign_score,
                                 col = list(Donor = donor_colors, Score = col_bw),
                                 annotation_name_gp = gpar(fontsize = 10),
                                 border = T)
bottom_anno.ha <- HeatmapAnnotation(Mut = as.matrix(hm_anno_df[,detected_mutations]),
                                    col = list(Mut = mut_colors),
                                    annotation_name_gp = gpar(fontsize = 10),
                                    border = T)

# Create Heatmap object
hm <- Heatmap(as.matrix(expr_mat),
              col = colItay(c(1:11))[3:11],
              cluster_rows = F,
              cluster_columns = T,
              row_names_gp = gpar(fontsize = 6),
              show_column_names = F,
              column_split = 2,
              top_annotation = top_anno.ha,
              bottom_annotation = bottom_anno.ha,
              name = "Expr",
              column_title = "Cells classified as pDC in samples without bone marrow involvement",
              column_title_gp = gpar(fontsize = 12),
              border = T,
              use_raster = T,
              raster_quality = 10)

pdf(paste0("9.2.3_Skin_only_pDC_heatmap", pdf_name, ".pdf"), width = 10, height = 6)
print(hm)
dev.off()






