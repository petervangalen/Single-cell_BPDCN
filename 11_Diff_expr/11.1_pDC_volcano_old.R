# Peter van Galen, 220617
# Compare expression between pDCs in healthy individuals, individuals with skin-only BPDCN, and individuals with clinically defined bone marrow involvement

library(tidyverse)
library(Seurat)
library(readxl)
library(ggrepel)
library(ggforce)
#library(data.table)

rm(list=ls())

# Set working directory
setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/Github/11_Diff_expr")

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")
mycol <- popcol.tib$hex
names(mycol) <- popcol.tib$pop

# Load Seurat objects
bm <- readRDS("../2_Annotate/BM_Seurat_CellTypes.rds")
seurat_files <- list.files("../7_XV-seq", pattern = "*.rds", full.names = T)
seu_bpdcn_ls <- lapply(seurat_files, function(x) readRDS(x))
names(seu_bpdcn_ls) <- gsub("_Seurat_Anno.rds", "", cutf(seurat_files, "/", f = 3))

# Merge
seu_all_ls <- c(BM = list(bm), seu_bpdcn_ls)

# Add more metadata; split into three groups based on patient / clinical information
skin_only_samples <- c("Pt1Mrd", "Pt9Dx", "Pt10Dx", "Pt12Dx")
bm_involvement_samples <- c("Pt1Dx", "Pt5Dx", "Pt10Rel", "Pt12Rel")

seu_all_ls <- lapply(names(seu_all_ls), function(n) {
  x <- seu_all_ls[[n]]
  if (n == "BM") {
    x$orig.ident2 <- cutf(x$replicate, d = "\\.")
    x$donor_group <- "healthy_bm" }
  else if (n %in% skin_only_samples) {
    x$orig.ident2 <- x$orig.ident
    x$donor_group <- "skin_only" }
  else if (n %in% bm_involvement_samples) {
    x$orig.ident2 <- x$orig.ident
    x$donor_group <- "bm_involvement" }
  return(x)
} )
names(seu_all_ls) <- c("BM", names(seu_bpdcn_ls))


# pDC proportion ----------------------------------------------------------------------------------

# Calculate proportion of pDCs over HSPC+Erythroid+Myeloid cells
celltypes.ch <- c("pDC", "cDC", "ncMono", "Mono", "ProMono", "Prog", "HSC", "EarlyE", "LateE")
seu_all_split_ls <- c(SplitObject(seu_all_ls[[1]], split = "orig.ident2"),
                      seu_all_ls[2:length(seu_all_ls)])
pdcs_num <- unlist( lapply(seu_all_split_ls, function(x)
  ncol(subset(x, CellType == "pDC")) / ncol(subset(x, CellType %in% celltypes.ch))) )
# Make tibble with results
pdc_proportion_tib <- tibble(Sample = names(pdcs_num), pDCs = pdcs_num) %>%
  mutate(Group = case_when(grepl("^BM", Sample) ~ "healthy_bm",
                           Sample %in% skin_only_samples ~ "skin_only",
                           Sample %in% bm_involvement_samples ~ "bm_involvement")) %>%
  mutate(Group = factor(Group, levels = c("healthy_bm", "skin_only", "bm_involvement")))

# Plot
pdf("1_pDC_Proportion.pdf", width = 6, height = 6)
pdc_proportion_tib %>%
  ggplot(aes(x = Group, y = pDCs, label = Sample)) +
  geom_point(size = 2) +
  geom_label_repel() +
  ylab("pDCs (% of HSPC+erythroid+myeloid)") +
  theme_bw() +
  theme(aspect.ratio = 1,
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10, color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        panel.grid = element_blank())
dev.off()

# The proportion of pDCs tends to be higher in patients with involvement
t.test(filter(pdc_proportion_tib, Group == "healthy_bm")$pDCs,
       filter(pdc_proportion_tib, Group == "skin_only")$pDCs, var.equal = T)$p.value
t.test(filter(pdc_proportion_tib, Group == "healthy_bm")$pDCs,
       filter(pdc_proportion_tib, Group == "bm_involvement")$pDCs, var.equal = T)$p.value
t.test(filter(pdc_proportion_tib, Group == "skin_only")$pDCs,
       filter(pdc_proportion_tib, Group == "bm_involvement")$pDCs, var.equal = T)$p.value


# Check for technical artifacts -------------------------------------------------------------------

# Make one object with all pDCs & check that it makes sense
seu_all <- merge(seu_all_ls[[1]], seu_all_ls[2:length(seu_all_ls)])
seu_all$donor_group <- factor(seu_all$donor_group, levels = c("healthy_bm", "skin_only", "bm_involvement"))
as_tibble(seu_all@meta.data) %>% select(orig.ident, orig.ident2, donor_group) %>%
  group_by(orig.ident, orig.ident2, donor_group) %>% count()

# Subset for pDCs
all_pdcs <- subset(seu_all, CellType == "pDC")

# Plot data complexity
VlnPlot(all_pdcs, features = c("nFeature_RNA", "nCount_RNA"), group.by = "orig.ident2")

# There are substantial differences in ribosomal protein genes in Seq-Well data (BM02, Pt9Dx) compared to 10x data
rp.ch <- rownames(seu_all)[grepl("^RPS|^RPL", rownames(seu_all))]
pdf("2_RP_genes.pdf", width = 6, height = 6)
sapply(rp.ch, function(x) print( VlnPlot(seu_all, features = x, group.by = "orig.ident2") ))
dev.off()

# This technology-driven artifact causes a number of issues later on; remove RPS and RPL genes.
all_pdcs <- subset(all_pdcs, features = setdiff(rownames(all_pdcs), rp.ch))
all_pdcs <- NormalizeData(all_pdcs)


# UMAP --------------------------------------------------------------------------------------------

# Similar to 1_Seurat_Harmony.R
all_pdcs <- FindVariableFeatures(all_pdcs)
all_pdcs <- ScaleData(all_pdcs, features = rownames(all_pdcs))
all_pdcs <- RunPCA(all_pdcs, features = VariableFeatures(object = all_pdcs))
all_pdcs <- RunUMAP(all_pdcs, reduction = "pca", dims = 1:50, seed = 42)

pdf("3_UMAPs.pdf", width = 6, height = 6)
print(
DimPlot(all_pdcs, reduction = "umap", group.by = "orig.ident", pt.size = 0.8, shuffle = T) +
  scale_color_manual(values = mycol[unique(all_pdcs$orig.ident)]) + 
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
)
print(
DimPlot(all_pdcs, reduction = "umap", group.by = "donor_group", pt.size = 0.8, shuffle = T) +
  scale_color_manual(values = c("#708090", "#81c57a", "#f08080")) + 
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
)
dev.off()


# Differential gene expression --------------------------------------------------------------------

# There are widely different numbers of cells for different samples
as_tibble(all_pdcs@meta.data) %>% select(orig.ident, orig.ident2, donor_group) %>% group_by_all() %>% count()

# Make matrices for gene expression comparison. For BPDCNs, 50 pDCs max
bm_mat <- subset(all_pdcs, donor_group == "healthy_bm") %>% GetAssayData(slot = "data")

skin_only_pdcs_ls <- subset(all_pdcs, donor_group == "skin_only") %>% SplitObject(split.by = "orig.ident")
skin_only_pdcs50 <- lapply(skin_only_pdcs_ls, function(x) x[,sample(colnames(x), size = min(50, ncol(x)))])
skin_only_pdcs_merge <- merge(skin_only_pdcs50[[1]], skin_only_pdcs50[2:length(skin_only_pdcs50)])
skin_only_mat <- GetAssayData(skin_only_pdcs_merge, slot = "data")

bm_involvement_pdcs_ls <- subset(all_pdcs, donor_group == "bm_involvement") %>% SplitObject(split.by = "orig.ident")
bm_involvement_pdcs50 <- lapply(bm_involvement_pdcs_ls, function(x) x[,sample(colnames(x), size = min(50, ncol(x)))])
bm_involvement_pdcs_merge <- merge(bm_involvement_pdcs50[[1]], bm_involvement_pdcs50[2:length(bm_involvement_pdcs50)])
bm_involvement_mat <- GetAssayData(bm_involvement_pdcs_merge, slot = "data")

# Make table with gene info
diffgenes.tib <- tibble(gene = rownames(bm_mat),
                        healthy_bm = rowMeans(bm_mat),
                        skin_only = rowMeans(skin_only_mat),
                        bm_involvement = rowMeans(bm_involvement_mat))
diffgenes.tib <- mutate(diffgenes.tib, logFC = bm_involvement - healthy_bm)

# Add P-values
diffgenes.tib$pval <- p.adjust(sapply(diffgenes.tib$gene, function(z)
  wilcox.test(x = bm_involvement_mat[z,], y = bm_mat[z,])$p.value), method = "fdr")
FeaturePlot(all_pdcs, features = "TCL1A") + theme(aspect.ratio = 1)

# Order and save
diffgenes.tib <- arrange(diffgenes.tib, desc(logFC))
write.table(diffgenes.tib, file = "DiffGenes.txt", quote = F, sep = "\t", row.names = F)
write.table(select(diffgenes.tib, gene, logFC), file = "DiffGenes.rnk", quote = F, sep = "\t", row.names = F, col.names = F)


# Run GSEA (outside of R) -------------------------------------------------------------------------

# Use the gmt file "201204_CombinedSigns.gmt" and the rnk file DiffGenes.rnk
# Pre-Ranked GSEA was run using default settings
# This GSEA shows significant enrichment KEGG_OXIDATIVE_PHOSPHORYLATION, VILLANI_BPDCN_UP, and TAKAO_RESPONSE_TO_UVB_RADIATION_UP


# Project signatures ------------------------------------------------------------------------------

signs.tib <- read_excel("11_Signatures.xlsx")[-1,]
signs.ls <- lapply(signs.tib, intersect, rownames(GetAssayData(all_pdcs)))

# Add module scores
pdf("4_ModuleScores.pdf", width = 6, height = 6)
for (n in names(signs.ls)) {
  all_pdcs <- AddModuleScore(object = all_pdcs, features = signs.ls[n], name = n)
  colnames(all_pdcs@meta.data) <- gsub(str_c(n, "1$"), str_c(n, "_Score"), colnames(all_pdcs@meta.data))
  print(
    FeaturePlot(all_pdcs, features = paste0(n, "_Score"), pt.size = 0.8, order = T) +
      ggtitle(n) +
      theme(aspect.ratio = 1)
  )
}
dev.off()

as_tibble(all_pdcs@meta.data, rownames = "cell") %>%
  mutate(UMAP1 = all_pdcs@reductions$umap@cell.embeddings[,1],
         UMAP2 = all_pdcs@reductions$umap@cell.embeddings[,2]) %>%
  ggplot(aes_string(y = paste0(names(signs.ls)[1], "_Score"), x = "orig.ident2", color = "donor_group")) +
  geom_sina() +
  theme_bw() +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# NEXT:

# 1. Look at CD123 expression
#grep("IL3RA", rownames(all_pdcs))

# 2. Quantify the pDC signature in Progenitors between patient groups and in Progenitors with or without TET2 





