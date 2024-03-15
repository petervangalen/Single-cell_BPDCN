# Peter van Galen, 221213
# Load all relevant/final data and visualize gene expression

library(tidyverse)
library(Seurat)
library(readxl)
library(janitor)
library(ggforce)
library(ggrastr)
library(cowplot)
library(ggrepel)
#library(data.table)
#library(viridis)

rm(list=ls())

# Set working directory
setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/13_Misc")

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")
cell_colors <- popcol.tib$hex[1:21]
names(cell_colors) <- popcol.tib$pop[1:21]
sample_colors <- popcol.tib$hex[23:40]
names(sample_colors) <- popcol.tib$pop[23:40]
group_colors <- popcol.tib$hex[41:43]
names(group_colors) <- popcol.tib$pop[41:43]
mut_colors <- popcol.tib$hex[44:46]
names(mut_colors) <- popcol.tib$pop[44:46]
malignant_stage <- popcol.tib$hex[47:50]
names(malignant_stage) <- popcol.tib$pop[47:50]


# Load data ---------------------------------------------------------------------------------------

# Load Seurat objects. These can be downloaded from https://vangalenlab.bwh.harvard.edu/resources/bpdcn-uv/
seurat_files <- list.files("../04_XV-seq", pattern = "*.rds", full.names = T)
seu_ls <- lapply(seurat_files, function(x) readRDS(x))
names(seu_ls) <- cutf(basename(seurat_files), d = "_")
seu <- merge(seu_ls[[1]], seu_ls[2:length(seu_ls)])
seu$orig.ident2 <- ifelse(grepl("BM", seu$orig.ident), yes = cutf(seu$replicate, d = "\\."), no = seu$orig.ident)

# Add metadata from 11.3_Classify_malignant_BPDCN.R
MalignantCalls_df <- read.table("../11_pDC_expr/11.3_MalignantCalls_Final.txt", header = T, row.names = "cell")
all(rownames(seu@meta.data) == rownames(MalignantCalls_df))
seu$is_malignant <- MalignantCalls_df$is_malignant

# Load genotyping information
genotyping_tables.tib <- read_excel("../04_XV-seq/XV-seq_overview.xlsx")
# Replace different MTAP entries with one, just as in 4.1_Add_GoT-XV_to_Seurat.R
genotyping_tables.tib$Mutation <- gsub("MTAP.rearr.*", "MTAP.rearr", genotyping_tables.tib$Mutation)

# Ensure nice ordering in visualizations
seu$orig.ident <- factor(seu$orig.ident, levels = c("BM", "Pt1Dx", "Pt1Rem", "Pt5Dx", "Pt9Dx", "Pt10Dx", "Pt10Rel",
                                                    "Pt12Dx", "Pt12Rel", "Pt14Dx", "Pt15Dx", "Pt16Dx"))
seu$CellType <- factor(seu$CellType, levels = names(cell_colors))

# UMAP of normal cell  types
p1 <- as_tibble(seu@meta.data) %>% filter(orig.ident == "BM") %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = CellType)) +
  geom_point_rast(size = 0.1) +
  scale_color_manual(values = cell_colors) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.text = element_text(color = "black"))
p1

# Plot expression in healthy donors ---------------------------------------------------------------

# Subset Seurat object
seu_healthy <- subset(seu, orig.ident == "BM")

# Make metadata tibble
metadata_healthy_tib <- as_tibble(seu_healthy@meta.data, rownames = "cell")

# Add expression values
gene <- "GAPDH"
gene <- "DNMT3A"
gene <- "HOXA9"
gene <- "PIK3R5"
gene <- "SPI1"
metadata_healthy_tib[,"expression"] <- LayerData(seu_healthy, layer = "data")[gene,]

# Bar plot of mean expression split by cell type
p2 <- metadata_healthy_tib %>% group_by(CellType) %>%
  summarize(n = n(), mean_expr = mean(expression)) %>%
  ggplot(aes(x = CellType, y = mean_expr, fill = CellType)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = paste0(cell_colors, "CC")) +
  ylab("Mean normalized expression log(TP10K+1)") +
  ggtitle(gene) +
  theme_bw() +
  theme(aspect.ratio = 0.6,
        panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())

# Sina plot of expression split by cell type
p3 <- metadata_healthy_tib %>%
  ggplot(aes(x = CellType, y = expression, color = CellType)) +
  geom_violin(scale = "width") +
  geom_sina(aes(group = CellType), scale = "width", size = 0.5) +
  scale_color_manual(values = paste0(cell_colors, "CC")) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  ylab("Normalized expression log(TP10K+1)") +
  ggtitle(label = gene) +
  theme_bw() +
  theme(aspect.ratio = 0.6,
        panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())

# UMAP with gene expression
p4 <- metadata_healthy_tib %>% filter(orig.ident == "BM") %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = expression)) +
  geom_point_rast(size = 0.1) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.text = element_text(color = "black"))

pdf(paste0(gene, "_expr_healthy.pdf"), width = 16, height = 8)
plot_grid(p2, p3, p1, p4, ncol = 2)
dev.off()


# Plot gene expression in healthy donors compared to BPDCN patients -------------------------------

# Make metadata tibble
metadata_tib <- as_tibble(seu@meta.data, rownames = "cell") %>%
  mutate(is_malignant = factor(is_malignant, levels = c("Healthy", "Premalignant", "Malignant", "Other")))

# Add expression values
gene <- "BCL2"
gene <- "CXCR4"
gene <- "PIK3R5"
gene <- "SPI1"
metadata_tib[,"expression"] <- LayerData(seu, layer = "data")[gene,]

# Gene expression in pDCs according to malignant state
p5 <- metadata_tib %>% filter(CellType == "pDC") %>%
  group_by(orig.ident, is_malignant) %>%
  summarize(n = n(), mean_expr = mean(expression)) %>%
  ggplot(aes(x = orig.ident, y = mean_expr, fill = is_malignant)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = n), position = position_dodge(width = 1), vjust = -0.25) +
  scale_fill_manual(values = paste0(malignant_stage, "CC")) +
  ylab("Mean normalized expression log(TP10K+1)") +
  ggtitle(paste0(gene, " in pDCs")) + # make sure this is correct
  theme_bw() +
  theme(aspect.ratio = 0.6,
        panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())

# Sina/violin plot version of this
p6 <- metadata_tib %>% filter(CellType == "pDC") %>%
  ggplot(aes(x = orig.ident, y = expression, color = is_malignant)) +
  geom_violin(scale = "width", fill = NA) +
  geom_sina(aes(color = is_malignant), scale = "width") +
  scale_color_manual(values = malignant_stage) +
  ggtitle(paste0(gene, " in pDCs")) + # make sure this is correct
  ylab("Normalized expression log(TP10K+1)") +
  theme_bw() +
  theme(aspect.ratio = 0.6,
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())

pdf(paste0(gene, "_expr_withBPDCN.pdf"), width = 16, height = 4)
plot_grid(p5, p6, ncol = 2)
dev.off()

  

# Correlate
#gene1 <- "PIK3R5"
#gene2 <- "SPI1"
#metadata_tib[,"expression1"] <- LayerData(seu, layer = "data")[gene1,]
#metadata_tib[,"expression2"] <- LayerData(seu, layer = "data")[gene2,]

#metadata_tib %>%
#  ggplot(aes(x = expression1, y = expression2, color = CellType)) +
#  geom_point() +
#  geom_smooth(method = "lm", aes(color = NULL)) +
#  scale_color_manual(values = cell_colors) +
#  xlab(paste(gene1, "expression")) +
#  ylab(paste(gene2, "expression")) +
#  theme_bw() +
#  theme(aspect.ratio = 1,
#        panel.grid = element_blank())


# NOTE: THE SCRIPTS BELOW ARE UNDERDEVELOPED



# Age correlations --------------------------------------------------------------------------------

# The scripts below describe a rudimentary attempt to correlate patient age to HSPC gene expression
# Later, I wrote ~/DropboxMGB/GitHub/blood-aging/ (private repository as of 230910)

# Add healthy donor age
metadata_healthy_tib <- metadata_healthy_tib %>%
  mutate(age = case_when(orig.ident2 == "BM1" ~ 31,
                         orig.ident2 == "BM2" ~ 37,
                         orig.ident2 == "BM3" ~ 45,
                         orig.ident2 == "BM4" ~ 75,
                         orig.ident2 == "BM5" ~ 74,
                         orig.ident2 == "BM6" ~ 52))
metadata_healthy_hspc <- filter(metadata_healthy_tib, CellType == "HSPC")
seu_healthy_hspc <- subset(seu_healthy, CellType == "HSPC")
all(metadata_healthy_hspc$cell == colnames(seu_healthy_hspc))
tabyl(metadata_healthy_hspc$orig.ident2)

# Method 1: Correlate age to gene expression (using all cells)
expr_mat <- as.matrix(LayerData(seu_healthy_hspc, layer = "data"))
age_expr_cor <- tibble(gene = rownames(seu_healthy_hspc),
                       cor = apply(expr_mat, 1, cor, metadata_healthy_hspc$age),
                       pval = apply(expr_mat, 1, function(x) cor.test(x, metadata_healthy_hspc$age)$p.val))
age_expr_cor <- na.omit(age_expr_cor)
age_expr_cor <- age_expr_cor %>% arrange(cor) %>%
  mutate(index = 1:nrow(age_expr_cor), .before = 1)

# Method 2: Correlate age to gene expression (collapse by donor)
ave_sue <- AverageExpression(seu_healthy_hspc, group.by = "orig.ident2", return.seurat = T)
expr_norm <- as.matrix( LayerData(ave_sue, layer = "data") )
age_expr_cor <- tibble(gene = rownames(seu_healthy_hspc),
                       cor = apply(expr_norm, 1, cor, c(31, 37, 45, 75, 74, 52)),
                       pval = apply(expr_norm, 1, function(x) cor.test(x, c(31, 37, 45, 75, 74, 52))$p.val))
age_expr_cor <- na.omit(age_expr_cor)
age_expr_cor <- age_expr_cor %>% arrange(cor) %>%
  mutate(index = 1:nrow(age_expr_cor), .before = 1)
plot(x = c(31, 37, 45, 75, 74, 52), y = expr_norm["DDIT3",])

# Select genes to visualize from pang,weissman2011 supplement
highlight_genes <- c("FIGN", "PRR16", "CDH7", "SLC1A6", "TRIM13", "MYC", "EGR1", "LOX")
# Select genes to visualize from personal interest
highlight_genes <- c("ID2", "DNAJB9", "FOS", "JUN", "DDIT3", "IL1B", "EGR1", "CEBPA")

age_expr_cor %>%
  ggplot(aes(x = cor, y = -log10(pval), color = gene %in% highlight_genes)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05)) +
  geom_text_repel(data = filter(age_expr_cor, gene %in% highlight_genes),
          aes(x = cor, y = -log10(pval), label = gene), min.segment.length = 0) +
  theme_bw() +
  theme(aspect.ratio = 1)

# Look at genes
age_expr_cor %>% arrange(-index) %>% print(n = 50)
# Based on biological understanding, I think Method 1 yields better results (e.g. EGR1). However, in later scripts (GitHub/blood-aging/) and following an email exchange with Nathan Salomonis, I decided to collapse values by donor, i.e. comparing pseudobulks. While reducing power, this prevents the results from being dominated by a single donor.

# To perform GSEA
#write_tsv(age_expr_cor[,2:3], file = "age_cor.rnk", col_names = F)


