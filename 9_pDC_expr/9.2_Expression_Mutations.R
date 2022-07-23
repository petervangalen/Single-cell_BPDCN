# Peter van Galen, 220722
# Relate gene expression in pDCs and malignant BPDCN cells to their mutational status

library(tidyverse)
library(Seurat)
library(readxl)
library(ggforce)
library(gridExtra)
library(viridis)
#library(data.table)
#library(ggrepel)

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

# Load Seurat object
seurat_files <- list.files("../4_XV-seq", pattern = "*.rds", full.names = T)
seu_ls <- lapply(seurat_files, function(x) readRDS(x))
names(seu_ls) <- cutf(basename(seurat_files), d = "_")
seu <- merge(seu_ls[[1]], seu_ls[2:length(seu_ls)])
seu$orig.ident2 <- ifelse(grepl("BM", seu$orig.ident), yes = cutf(seu$replicate, d = "\\."), no = seu$orig.ident)

# Or pDCs alone
pdcs <- readRDS("pDCs.rds")

# Load genotyping information
genotyping_tables.tib <- read_excel("../4_XV-seq/XV-seq_overview.xlsx")
# Replace different MTAP primers with one, just as in 4.1_Add_GoT-XV_to_Seurat.R
genotyping_tables.tib$Mutation <- gsub("MTAP.rearr.*", "MTAP.rearr", genotyping_tables.tib$Mutation)


# Approach 1: overlay mutation status on UMAPs based on signature genes ---------------------------

# What does the data look like?
p1 <- DimPlot(pdcs, group.by = "orig.ident", cols = donor_colors[unique(pdcs$orig.ident)]) + theme(aspect.ratio = 1)
p2 <- DimPlot(pdcs, group.by = "donor_group", cols = group_colors) + theme(aspect.ratio = 1)
p3 <- FeaturePlot(pdcs, features = "VILLANI_BPDCN_UP_Score") + scale_color_viridis() + theme(aspect.ratio = 1)
grid.arrange(p1, p2, p3, ncol = 2)






metadata_tib <- as_tibble(pdcs@meta.data, rownames = "cell")





metadata_tib %>% filter(donor_group != "bm_involvement") %>% 
  ggplot(aes(x = orig.ident, y = VILLANI_BPDCN_UP_Score)) +
  geom_hline(yintercept = 0.1) +
  geom_sina()

metadata_tib %>% filter(donor_group != "bm_involvement", VILLANI_BPDCN_UP_Score > 0.1) %>% view
metadata_tib %>% filter(orig.ident == "Pt9Dx", VILLANI_BPDCN_UP_Score > 0.1) %>% view # no call for all three genotyped mutations
metadata_tib %>% filter(orig.ident == "Pt12Dx", VILLANI_BPDCN_UP_Score > 0.1) %>% view # no call for all four genotyped mutations
metadata_tib %>% filter(orig.ident == "Pt10Dx", VILLANI_BPDCN_UP_Score > 0.1) %>% view

metadata_tib %>% filter(orig.ident == "Pt10Dx") %>% mutate(ScorePos = VILLANI_BPDCN_UP_Score > 0.1) %>% .$ScorePos %>% table

seu_ls$Pt9Dx
seu_ls$Pt10Dx
seu_ls$Pt12Dx



# ------------------------------------------------------------------------------------------------

# Wrangle
seu_skin_only <- subset(pdcs, donor_group == "skin_only")
skin_only_metadata <- as_tibble(seu_skin_only@meta.data, rownames = "cell") %>%
  dplyr::select(cell, orig.ident, all_of(prog_mut), VILLANI_BPDCN_UP_Score)

# Which mutations to assess? Let's start with all progression mutations, which may identify malignant cells
prog_mut <- genotyping_tables.tib %>% filter(`Founder or progression mutation` == "Progression") %>% .$Mutation
# However, only Pt10Dx has any  mutations. 
#skin_only_metadata$n_mut <- apply(skin_only_metadata, 1, function(x) sum(grepl("mutant", x)))
#skin_only_metadata %>% filter(n_mut > 0) %>% .$orig.ident
# In Pt10Dx, MTAP.rearr and RAB9A.3pUTR are particularly useful, because they only have the mutated allele. So let's just use those.
#prog_mut <- c("MTAP.rearr", "RAB9A.3pUTR")

# Add mut/wt column to metadata
skin_only_metadata$call <- ifelse(apply(dplyr::select(skin_only_metadata, all_of(prog_mut)), 1,
                                        function(x) sum(grepl("wildtype", x) > 0)), yes = "wildtype", no = "no_call")
skin_only_metadata$call <- ifelse(apply(dplyr::select(skin_only_metadata, all_of(prog_mut)), 1,
                                        function(x) sum(grepl("mutant", x) > 0)), yes = "mutant", no = skin_only_metadata$call)

# Refine Villani signature genes ------------------------------------------------------------------

#signs.tib <- read_excel("../9_pDC_Expr/Signatures.xlsx")[-1,]
#villani_genes <- signs.tib$VILLANI_BPDCN_UP
#villani_genes <- intersect(villani_genes, rownames(seu_skin_only))

# Identify genes that correlate with the Villani score
expr_mat <- GetAssayData(pdcs, slot = "data")
villani_score <- pdcs$VILLANI_BPDCN_UP_Score
my_cor <- apply(expr_mat, 1, function(x) cor(x, villani_score))
my_cor <- sort(my_cor, decreasing = T)
plot(my_cor, pch = ".")
abline(v = 100, col = "red")

# Select the top 100 as a new signature
bpdcn_sign <- names(my_cor[1:100])


# Skin-only pDC UMAPs -----------------------------------------------------------------------------

# Calculate UMAP coordinates for pDCs using this signature
seu_skin_only <- NormalizeData(seu_skin_only)
#seu_skin_only <- FindVariableFeatures(seu_skin_only)
seu_skin_only <- ScaleData(seu_skin_only, features = rownames(seu_skin_only))
seu_skin_only <- RunPCA(seu_skin_only, features = bpdcn_sign)
seu_skin_only <- RunUMAP(seu_skin_only, reduction = "pca", dims = 1:10)

# Add some more metadata
seu_skin_only <- CellCycleScoring(seu_skin_only, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
if (all(colnames(seu_skin_only) == skin_only_metadata$cell)) { seu_skin_only$call <- skin_only_metadata$call }

# Plot
p1 <- DimPlot(seu_skin_only, group.by = "orig.ident", cols = donor_colors[unique(seu_skin_only$orig.ident)]) + theme(aspect.ratio = 1)
p2 <- DimPlot(seu_skin_only, group.by = "call", cols = mut_colors) + theme(aspect.ratio = 1)

pdf("9.2.1_Skin_only_pDCs.pdf", width = 10, height = 5)
grid.arrange(p1, p2, ncol = 2)
dev.off()

FeaturePlot(seu_skin_only, features = "TCL1A") + scale_color_viridis() + theme(aspect.ratio = 1)
FeaturePlot(seu_skin_only, features = "S.Score") + scale_color_viridis() + theme(aspect.ratio = 1)
FeaturePlot(seu_skin_only, features = "G2M.Score") + scale_color_viridis() + theme(aspect.ratio = 1)
FeaturePlot(seu_skin_only, features = "VILLANI_BPDCN_UP_Score") + scale_color_viridis() + theme(aspect.ratio = 1)


