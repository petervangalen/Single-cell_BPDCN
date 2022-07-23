# Peter van Galen, 220718
# Assess gene expression in pDCs and BPDCN tumor cells

library(tidyverse)
library(Seurat)
library(readxl)
library(data.table)
library(ggrepel)
library(ggforce)
library(gridExtra)
library(viridis)

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

# Load Seurat objects
seurat_files <- list.files("../4_XV-seq", pattern = "*.rds", full.names = T)
seu_ls <- lapply(seurat_files, function(x) readRDS(x))
names(seu_ls) <- cutf(basename(seurat_files), d = "_")
seu <- merge(seu_ls[[1]], seu_ls[2:length(seu_ls)])
seu$orig.ident2 <- ifelse(grepl("BM", seu$orig.ident), yes = cutf(seu$replicate, d = "\\."), no = seu$orig.ident)

# Define sample groups
metadata_tib <- as_tibble(seu@meta.data, rownames = "cell")
skin_only_samples <- unique(filter(metadata_tib, donor_group == "skin_only")$orig.ident)
bm_involvement_samples <- unique(filter(metadata_tib, donor_group == "bm_involvement")$orig.ident)


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
pdcs <- NormalizeData(pdcs)
pdcs <- FindVariableFeatures(pdcs)
pdcs <- ScaleData(pdcs, features = rownames(pdcs))
pdcs <- RunPCA(pdcs, features = VariableFeatures(object = pdcs))
# I'm skipping these Harmony lines
#pdcs <- RunHarmony(object = pdcs, group.by.vars = "tech", reduction = "pca", plot_convergence = T)
#pdcs <- RunUMAP(pdcs, reduction = "harmony", dims = 1:50)
pdcs <- RunUMAP(pdcs, reduction = "pca", dims = 1:50)

# Make donor group a factor & score cell cycle
pdcs$donor_group <- factor(pdcs$donor_group, levels = c("healthy_bm", "skin_only", "bm_involvement"))
pdcs <- CellCycleScoring(pdcs, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)

# Plot
pdf("9.1.3_pDC_UMAPs.pdf", width = 20, height = 8)
p1 <- DimPlot(pdcs, group.by = "orig.ident", cols = donor_colors[unique(pdcs$orig.ident)]) + theme(aspect.ratio = 1)
p2 <- DimPlot(pdcs, group.by = "donor_group", cols = group_colors[unique(pdcs$donor_group)]) + theme(aspect.ratio = 1) 
p3 <- FeaturePlot(pdcs, features = "IL3RA") + scale_color_viridis() + theme(aspect.ratio = 1)
p4 <- FeaturePlot(pdcs, features = "TCL1A") + scale_color_viridis() + theme(aspect.ratio = 1)
p5 <- FeaturePlot(pdcs, features = "S.Score") + scale_color_viridis() + theme(aspect.ratio = 1)
p6 <- FeaturePlot(pdcs, features = "G2M.Score") + scale_color_viridis() + theme(aspect.ratio = 1)
grid.arrange(p1, p2, p3, p4, p5, p6, layout_matrix = matrix(1:6, ncol = 3, byrow = T))
dev.off()


# Differential gene expression --------------------------------------------------------------------

# There are widely different numbers of cells for different samples
as_tibble(pdcs@meta.data) %>% dplyr::select(orig.ident, orig.ident2, donor_group) %>% group_by_all() %>% count()

# Make matrices for gene expression comparison. Use 20 pDCs max per sample (to correct for the different numbers)
healthy_bm_pdcs_ls <- subset(pdcs, donor_group == "healthy_bm") %>% SplitObject(split.by = "orig.ident2")
healthy_bm_pdcs20 <- lapply(healthy_bm_pdcs_ls, function(x) x[,sample(colnames(x), size = min(20, ncol(x)))])
healthy_bm_pdcs_merge <- merge(healthy_bm_pdcs20[[1]], healthy_bm_pdcs20[2:length(healthy_bm_pdcs20)])
healthy_bm_mat <- GetAssayData(healthy_bm_pdcs_merge, slot = "data")

skin_only_pdcs_ls <- subset(pdcs, donor_group == "skin_only") %>% SplitObject(split.by = "orig.ident")
skin_only_pdcs20 <- lapply(skin_only_pdcs_ls, function(x) x[,sample(colnames(x), size = min(20, ncol(x)))])
skin_only_pdcs_merge <- merge(skin_only_pdcs20[[1]], skin_only_pdcs20[2:length(skin_only_pdcs20)])
skin_only_mat <- GetAssayData(skin_only_pdcs_merge, slot = "data")

bm_involvement_pdcs_ls <- subset(pdcs, donor_group == "bm_involvement") %>% SplitObject(split.by = "orig.ident")
bm_involvement_pdcs20 <- lapply(bm_involvement_pdcs_ls, function(x) x[,sample(colnames(x), size = min(20, ncol(x)))])
bm_involvement_pdcs_merge <- merge(bm_involvement_pdcs20[[1]], bm_involvement_pdcs20[2:length(bm_involvement_pdcs20)])
bm_involvement_mat <- GetAssayData(bm_involvement_pdcs_merge, slot = "data")

# Make table with gene info
diffgenes.tib <- tibble(gene = rownames(healthy_bm_mat),
                        healthy_bm = rowMeans(healthy_bm_mat),
                        skin_only = rowMeans(skin_only_mat),
                        bm_involvement = rowMeans(bm_involvement_mat))
diffgenes.tib <- mutate(diffgenes.tib, logFC.bm_involvement.healthy_bm = bm_involvement - healthy_bm)
diffgenes.tib <- mutate(diffgenes.tib, logFC.bm_involvement.skin_only = bm_involvement - skin_only)

# Are pDC gene expression changes correlated between when comparing bm_involvement to the other two conditions? Yes.
diffgenes.tib %>% ggplot(aes(x = logFC.bm_involvement.skin_only, y = logFC.bm_involvement.healthy_bm)) +
  coord_cartesian(xlim = c(-2.5,2.5), ylim = c(-2.5,2.5)) +
  geom_point() +
  theme(aspect.ratio = 1)

# Add P-values
diffgenes.tib$pval <- p.adjust(sapply(diffgenes.tib$gene, function(z)
  wilcox.test(x = bm_involvement_mat[z,], y = healthy_bm_mat[z,])$p.value), method = "fdr")

# Order, check & save
diffgenes.tib <- arrange(diffgenes.tib, desc(logFC.bm_involvement.healthy_bm))

write.table(diffgenes.tib, file = "DiffGenes.txt", quote = F, sep = "\t", row.names = F)
write.table(dplyr::select(diffgenes.tib, gene, logFC.bm_involvement.healthy_bm),
            file = "DiffGenes.rnk", quote = F, sep = "\t", row.names = F, col.names = F)


# Run GSEA (outside of R) -------------------------------------------------------------------------

# Use the gmt file "201204_CombinedSigns.gmt" and the rnk file DiffGenes.rnk
# Pre-Ranked GSEA was run using default settings
# This GSEA shows significant enrichment of BPDCN tumor signatures, and also some mitochondrial
# signatures, and downregulation of BIOCARTA_INFLAM_PATHWAY and PHONG_TNF_TARGETS_UP
# Some interesting genes that came up are BCL2 (up) and CXCR4 (down)

FeaturePlot(pdcs, features = c("BCL2", "CXCR4"), cols = viridis(100)) + theme(aspect.ratio = 1)


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
      ggtitle(n) +
      theme(aspect.ratio = 1)
  )
}
dev.off()

# Plot as bar graph
skin_only_samples <- c("Pt1Mrd", "Pt5Dx", "Pt9Dx", "Pt10Dx", "Pt12Dx")
bm_involvement_samples <- c("Pt1Dx", "Pt10Rel", "Pt12Rel", "Pt14Dx", "Pt15Dx", "Pt16Dx")

pdf("9.1.5_Modules_sina.pdf", width = 6, height = 5)
for (n in names(signs.ls)) {
  #n <- "VILLANI_BPDCN_UP" 
  print(
    as_tibble(pdcs@meta.data, rownames = "cell") %>% 
      mutate(orig.ident2 = factor(orig.ident2, levels = c(paste0("BM", 1:6), skin_only_samples, bm_involvement_samples))) %>%
      ggplot(aes_string(y = paste0(n, "_Score"), x = "orig.ident2", color = "donor_group")) +
      geom_sina(size = 0.3) +
      scale_color_manual(values = group_colors[unique(pdcs$donor_group)]) +
      theme_bw() +
      theme(aspect.ratio = 2/3,
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
            axis.text.y = element_text(color = "black"),
            panel.grid = element_blank())
  )
}
dev.off()

# Save
saveRDS(pdcs, file = "pDCs.rds")






#### OTHER THINGS I TRIED THAT WEREN'T VERY HELPFUL

# Quantify pDC priming in Progenitors with or without TET2 ----------------------------------------

prog <- subset(seu, CellType %in% c("HSPC", "GMP"))
prog <- NormalizeData(prog)
pdc_feats <- read.table("../2_Annotate/markerGenes.txt", header = T)[,"pDC"]
prog <- AddModuleScore(prog, features = list(pdc_feats), name = "pDC_Score")
prog@meta.data %>% head

prog_metdata.tib <- as_tibble(prog@meta.data, rownames = "cell")
apply(dplyr::select(prog_metdata.tib, contains("TET2")), 2, unique)
tet2_mut_prog <- prog_metdata.tib %>% filter(if_any(.cols = contains("TET2"), .fns = ~ grepl("mutant", .x))) %>% .$cell

prog_metdata.tib %>% mutate(TET2.mut = cell %in% tet2_mut_prog) %>%
  ggplot(aes_string(y = "pDC_Score1", x = "orig.ident2", color = "TET2.mut")) + # color = "donor_group"
  geom_sina() +
  theme_bw() +
  ylab("pDC signature score") +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# This analysis does not show that TET2 mutated progenitors are primed to pDC differentiation.


# Files for Volker's analysis, 220708 -------------------------------------------------------------
expr.mat <- matrix(nrow = nrow(seu.ls$BM), ncol=length(levels(seu.ls$BM$CellType)), dimnames = list(rownames(seu.ls$BM), levels(seu.ls$BM$CellType)))
for (n in levels(seu.ls$BM$CellType)) {
  expr.mat[,n] <- rowMeans( GetAssayData(subset(seu.ls$BM, CellType == n), slot = "data") )
}
expr.mat[c("HBD", "CD3D", "CD14"),]
write.table(expr.mat, file = "Mean_BM_expr.txt", quote = F, sep = "\t")
seu.ls$BM$CellType %>% table %>% data.frame %>% write.table(file = "CellNumbers.txt", sep = "\t", quote = F, row.names = F, col.names = F)



