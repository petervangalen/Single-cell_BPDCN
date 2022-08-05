# Peter van Galen, 220718
# Assess gene expression in pDCs and BPDCN tumor cells

library(tidyverse)
library(Seurat)
library(readxl)
library(data.table)
library(ggrepel)
library(ggforce)
library(cowplot)
library(viridis)

rm(list=ls())

# Set working directory
setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/9_pDC_Expr")

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")
donor_colors <- popcol.tib$hex[24:41]
names(donor_colors) <- popcol.tib$pop[24:41]
group_colors <- popcol.tib$hex[42:44]
names(group_colors) <- popcol.tib$pop[42:44]

# Load Seurat objects
seurat_files <- list.files("../4_XV-seq", pattern = "*.rds", full.names = T)
seu_ls <- lapply(seurat_files, function(x) readRDS(x))
names(seu_ls) <- cutf(basename(seurat_files), d = "_")
seu <- merge(seu_ls[[1]], seu_ls[2:length(seu_ls)])
seu$orig.ident2 <- ifelse(grepl("BM", seu$orig.ident), yes = cutf(seu$replicate, d = "\\."), no = seu$orig.ident)

# Define sample groups
metadata_tib <- as_tibble(seu@meta.data, rownames = "cell")
skin_only_samples <- unique(filter(metadata_tib, bm_involvement == "No")$orig.ident) %>% .[c(3,4,5,1,2)]
bm_involvement_samples <- unique(filter(metadata_tib, bm_involvement == "Yes")$orig.ident) %>% .[c(6,1,2,3,4,5)]

# Load genotyping information
genotyping_tables.tib <- read_excel("../4_XV-seq/XV-seq_overview.xlsx")
# Replace different MTAP entries with one, just as in 4.1_Add_GoT-XV_to_Seurat.R
genotyping_tables.tib$Mutation <- gsub("MTAP.rearr.*", "MTAP.rearr", genotyping_tables.tib$Mutation)


# pDC proportion ----------------------------------------------------------------------------------

# Calculate proportion of pDCs over HSPC+Erythroid+Myeloid cells
celltypes.ch <- c("pDC", "cDC", "ncMono", "Mono", "ProMono", "Prog", "HSC", "EarlyE", "LateE")
seu_split_ls <- SplitObject(seu, split = "orig.ident2")
pdcs_num <- unlist( lapply(seu_split_ls, function(x)
  ncol(subset(x, CellType == "pDC")) / ncol(subset(x, CellType %in% celltypes.ch))) )
# Make tibble with results
pdc_proportion_tib <- tibble(Sample = names(pdcs_num), pDCs = pdcs_num) %>%
  mutate(Group = case_when(grepl("^BM", Sample) ~ "healthy_bm",
                           Sample %in% skin_only_samples ~ "skin_only",
                           Sample %in% bm_involvement_samples ~ "bm_involvement")) %>%
  mutate(Group = factor(Group, levels = c("healthy_bm", "skin_only", "bm_involvement")))

# Plot
pdf("9.1.1_pDC_Proportion.pdf", width = 6, height = 6)
pdc_proportion_tib %>%
  ggplot(aes(x = Group, y = pDCs, label = Sample, color = Sample)) +
  geom_point(size = 2) +
  geom_label_repel(max.overlaps = 15, show.legend  = F) +
  scale_color_manual(values = donor_colors, na.value = donor_colors["BM"]) +
  ylab("pDCs (% of HSPC+erythroid+myeloid)") +
  theme_bw() +
  theme(aspect.ratio = 1,
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10, color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        panel.grid = element_blank())
dev.off()

# The proportion of pDCs is higher in patients with involvement
t.test(filter(pdc_proportion_tib, Group == "healthy_bm")$pDCs,
       filter(pdc_proportion_tib, Group == "skin_only")$pDCs, var.equal = T)$p.value
t.test(filter(pdc_proportion_tib, Group == "healthy_bm")$pDCs,
       filter(pdc_proportion_tib, Group == "bm_involvement")$pDCs, var.equal = T)$p.value
t.test(filter(pdc_proportion_tib, Group == "skin_only")$pDCs,
       filter(pdc_proportion_tib, Group == "bm_involvement")$pDCs, var.equal = T)$p.value


# Check for technical artifacts -------------------------------------------------------------------

# Subset for pDCs
pdcs <- subset(seu, CellType == "pDC")

# Plot data complexity
VlnPlot(pdcs, features = c("nFeature_RNA", "nCount_RNA"), group.by = "orig.ident2")

# There are substantial differences in ribosomal protein genes in Seq-Well data (BM6, Pt9Dx) compared to 10x data
rp.ch <- rownames(seu)[grepl("^RPS|^RPL", rownames(seu))]
# Commented out for now b/c it takes a few minutes
#pdf("9.1.2_RP_genes.pdf", width = 6, height = 6)
#sapply(rp.ch, function(x) print( VlnPlot(seu, features = x, group.by = "orig.ident2") ))
#dev.off()

# This technology-driven artifact causes some issues later on; remove RPS and RPL genes.
pdcs <- subset(pdcs, features = setdiff(rownames(pdcs), rp.ch))


# pDC UMAP ----------------------------------------------------------------------------------------

# Calculate UMAP coordinates for pDCs See 1_Seurat_Harmony for more information.
pdcs <- CellCycleScoring(pdcs, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
pdcs <- NormalizeData(pdcs)
pdcs <- FindVariableFeatures(pdcs)
pdcs <- ScaleData(pdcs, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(pdcs))
pdcs <- RunPCA(pdcs, features = VariableFeatures(object = pdcs))

# I'm skipping these Harmony lines
#pdcs <- RunHarmony(object = pdcs, group.by.vars = "tech", reduction = "pca", plot_convergence = T)
#pdcs <- RunUMAP(pdcs, reduction = "harmony", dims = 1:50)
pdcs <- RunUMAP(pdcs, reduction = "pca", dims = 1:20)

# Reverse coordinates for visualization purposes
pdcs$UMAP_1 <- pdcs@reductions$umap@cell.embeddings[,2]
pdcs$UMAP_2 <- pdcs@reductions$umap@cell.embeddings[,1]
pdcs@reductions$umap@cell.embeddings[,1] <- pdcs$UMAP_1
pdcs@reductions$umap@cell.embeddings[,2] <- pdcs$UMAP_2
# Start like this if you didn't change anything above & want to save time
#pdcs <- readRDS("pdcs.rds")

# Make relevant metadata into factors
pdcs$orig.ident <- factor(pdcs$orig.ident, levels = c("BM", skin_only_samples, bm_involvement_samples))
pdcs$bm_involvement <- factor(pdcs$bm_involvement, levels = c("HD", "No", "Yes"))

# Plot
pdf(paste0("9.1.3_pDC_UMAPs.pdf"), width = 7, height = 7)
p1 <- DimPlot(pdcs, group.by = "orig.ident", shuffle = T, cols = donor_colors[levels(pdcs$orig.ident)]) +
  theme_void() + theme(aspect.ratio = 1, panel.border = element_rect(color = "black", fill = NA))
p2 <- DimPlot(pdcs, group.by = "bm_involvement", shuffle = T, cols = group_colors[levels(pdcs$bm_involvement)]) +
  theme_void() + theme(aspect.ratio = 1, panel.border = element_rect(color = "black", fill = NA))
p3 <- FeaturePlot(pdcs, features = "IL3RA") + scale_color_viridis() + theme(aspect.ratio = 1) +
  theme_void() + theme(aspect.ratio = 1, panel.border = element_rect(color = "black", fill = NA))
p4 <- FeaturePlot(pdcs, features = "IGLL1") + scale_color_viridis() + theme(aspect.ratio = 1) +
  theme_void() + theme(aspect.ratio = 1, panel.border = element_rect(color = "black", fill = NA))
plot_grid(p1, p2, p3, p4, align = "v", rel_heights = c(1,1,1,1))
dev.off()


# Project signatures ------------------------------------------------------------------------------

signs.tib <- read_excel("Signatures.xlsx")[-1,]
signs.ls <- lapply(signs.tib, intersect, rownames(GetAssayData(pdcs)))

# Add module scores & plot on UMAP
pdf("9.1.4_ModuleScores.pdf", width = 6, height = 6)
for (n in names(signs.ls)) {
  # n <- names(signs.ls)[1]
  column_name <- paste0(n, "_Score")
  if (! column_name %in% colnames(pdcs@meta.data)) {
    pdcs <- AddModuleScore(object = pdcs, features = signs.ls[n], name = n)
    colnames(pdcs@meta.data) <- gsub(str_c(n, "1$"), column_name, colnames(pdcs@meta.data))
  }
  print(
    FeaturePlot(pdcs, features = column_name, pt.size = 0.6, order = T) +
      scale_color_viridis() +
      ggtitle(n) +
      theme_void() +
      theme(aspect.ratio = 1, panel.border = element_rect(color = "black", fill = NA))
  )
}
dev.off()

# Plot as bar graph
skin_only_samples <- c("Pt1Rem", "Pt5Dx", "Pt9Dx", "Pt10Dx", "Pt12Dx")
bm_involvement_samples <- c("Pt1Dx", "Pt10Rel", "Pt12Rel", "Pt14Dx", "Pt15Dx", "Pt16Dx")

pdf("9.1.5_Modules_sina.pdf", width = 6, height = 5)
for (n in names(signs.ls)) {
  #n <- "VILLANI_BPDCN_UP" 
  print(
    as_tibble(pdcs@meta.data, rownames = "cell") %>% 
      mutate(orig.ident2 = factor(orig.ident2, levels = c(paste0("BM", 1:6), skin_only_samples, bm_involvement_samples))) %>%
      ggplot(aes_string(y = paste0(n, "_Score"), x = "orig.ident2", color = "bm_involvement")) +
      geom_sina(size = 0.3) +
      scale_color_manual(values = group_colors[unique(pdcs$bm_involvement)]) +
      theme_bw() +
      theme(aspect.ratio = 2/3,
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
            axis.text.y = element_text(color = "black"),
            panel.grid = element_blank())
  )
}
dev.off()


# Add mutation calls and save Seurat object -------------------------------------------------------------

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

# Save
saveRDS(pdcs, file = "pDCs.rds")


