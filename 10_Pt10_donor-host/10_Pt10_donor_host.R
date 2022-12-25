# Peter van Galen, 220930
# Plot Patient 10 donor and host cell calls

library(tidyverse)
library(Seurat)
library(readxl)
library(data.table)
library(ggforce)
library(cowplot)

rm(list=ls())

# Set working directory
setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/10_Pt10_donor-host")

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")
cell_colors <- popcol.tib$hex[1:21]
names(cell_colors) <- popcol.tib$pop[1:21]
donor_colors <- popcol.tib$hex[24:41]
names(donor_colors) <- popcol.tib$pop[24:41]
mut_colors <- popcol.tib$hex[44:46]
names(mut_colors) <- popcol.tib$pop[44:46]

# Load Seurat objects
seurat_files <- list.files("../04_XV-seq", pattern = "*.rds", full.names = T)
seurat_files <- seurat_files[c(1,3)]
seu_ls <- lapply(seurat_files, function(x) readRDS(x))
names(seu_ls) <- cutf(basename(seurat_files), d = "_")
seu <- merge(seu_ls[[1]], seu_ls[2:length(seu_ls)])
seu$CellType <- factor(seu$CellType, levels = levels(seu_ls[[1]]$CellType))


# Prepare data ------------------------------------------------------------------------------------

# Load Volker's donor / recipient calls (from his email 220927)
genotype.dt <- fread("712-genotype.20220927.call.txt", col.names = c("CB","NA","donor","host","call"))
genotype.dt$CB <- paste0(cutf(genotype.dt$CB, d = "/|-", f = 3), "-",
                         gsub("BPDCN712", "Pt10Dx", gsub("BPDCN712R", "Pt10Rel", cutf(genotype.dt$CB, d = "-", f = 1))), ".",
                         cutf(genotype.dt$CB, d = "/|-", f = 2))

# Make factor to color new UMAP
orig.host.donor.ch <- seu$orig.ident
orig.host.donor.ch[names(orig.host.donor.ch) %in% genotype.dt[call == "host"]$CB] <- "Pt10Rel.Host"
orig.host.donor.ch[names(orig.host.donor.ch) %in% genotype.dt[call == "donor"]$CB] <- "Pt10Rel.Donor"
orig.host.donor.ch[names(orig.host.donor.ch) %in% genotype.dt[call == "doublet"]$CB] <- "Pt10Rel.Doublet"
orig.host.donor.ch[names(orig.host.donor.ch) %in% genotype.dt[is.na(genotype.dt$call)]$CB] <- "Pt10Rel.Inconclusive"

# Add meta.data
seu$orig.host.donor <- factor(x = orig.host.donor.ch, levels = c("BM", "Pt10Rel.Host", "Pt10Rel.Donor", "Pt10Rel.Doublet", "Pt10Rel.Inconclusive"))

# For the paper: 
table(seu$orig.ident)
table(seu$orig.host.donor)
# So 7,498 total cells, of which 148 doublets (1.97%) and 234 inconclusive (3.12%)
# We were able to determine the origin of (4453+2663)/7498 = 94.9% of cells

# Exclude doublets and inconclusive and Seq-Well data (BM6)
NormalRelapse <- subset(seu, subset = orig.host.donor %in% c("BM", "Pt10Rel.Host", "Pt10Rel.Donor") & tech == "TenX")

# Compare data quality (considering I did not do batch correction)
VlnPlot(NormalRelapse, features = c("nFeature_RNA", "nCount_RNA"), group.by = "orig.ident", pt.size = 0)
as_tibble(NormalRelapse@meta.data) %>%
  ggplot(aes(x = nCount_RNA, y = nFeature_RNA, color = orig.ident)) +
  geom_point(size = 0.3) +
  theme_bw() +
  theme(aspect.ratio = 1, panel.grid = element_blank())


# Cell type proportions ---------------------------------------------------------------------------

# After excluding unclear cells (above), 2663/(2663+4453)=37.4% are donor whereas 4453/(4453+2663)=62.6% are host
table(subset(NormalRelapse, subset = orig.ident == "Pt10Rel")$orig.host.donor)
# Of the 4453 Relapse Host cells, 4209 (94.5%) are pDCs
table(subset(NormalRelapse, subset = orig.host.donor == "Pt10Rel.Host")@meta.data[,"CellType"])

# Visualize
celltype_freq <- as_tibble(NormalRelapse@meta.data) %>%
  group_by(orig.host.donor, CellType) %>% count %>% group_by(orig.host.donor) %>%
  filter(orig.host.donor != "BM") %>%
  mutate(Percent = n/sum(n)*100)

pdf("10.1_CellTypeProportions.pdf", width = 8, height = 6)

celltype_freq %>% mutate(CellType = factor(CellType, levels = rev(levels(CellType)))) %>%
  ggplot(aes(x = Percent, y = orig.host.donor, fill = CellType)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cell_colors) +
  theme_bw() +
  theme(aspect.ratio = 0.3, panel.grid = element_blank(), axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"), axis.title.y = element_blank(),
        axis.text.y = element_text(size = 11))

dev.off()


# UMAP (not used for paper) -----------------------------------------------------------------------

# Use standard Seurat functions for dimensionality reduction
NormalRelapse <- NormalizeData(NormalRelapse)
NormalRelapse <- FindVariableFeatures(NormalRelapse)
NormalRelapse <- ScaleData(NormalRelapse)
NormalRelapse <- RunPCA(NormalRelapse, features = VariableFeatures(object = NormalRelapse))

pdf("10.2_UMAP_all_dimensions.pdf", width = 20, height = 4)
#for (ndim in c(10, 13, 20, 27, 30)) {
#print(ndim)
ndim <- 13
NormalRelapse <- RunUMAP(NormalRelapse, reduction = "pca", dims = 1:ndim)

NormalRelapse$UMAP_1 <- NormalRelapse@reductions$umap@cell.embeddings[,1]
NormalRelapse$UMAP_2 <- NormalRelapse@reductions$umap@cell.embeddings[,2]

metadata_tib <- as_tibble(NormalRelapse@meta.data, rownames = "cell")

p1 <- metadata_tib %>% mutate(rowsample = sample(nrow(metadata_tib))) %>% arrange(rowsample) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = CellType)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = cell_colors) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  ggtitle(paste0(ndim, " dimensions")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank())

p2 <- metadata_tib %>% mutate(rowsample = sample(nrow(metadata_tib))) %>% arrange(rowsample) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = orig.host.donor)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = c(BM = "#c0c0c0", Pt10Rel.Host = "#fa8072", Pt10Rel.Donor = "#da70d6")) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  ggtitle(paste0(ndim, " dimensions")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank())

p3 <- metadata_tib %>%
  mutate(orig.ident2 = ifelse(orig.ident == "BM", yes = cutf(replicate, d = "\\."), no = as.character(orig.host.donor))) %>%
  mutate(rowsample = sample(nrow(metadata_tib))) %>% arrange(rowsample) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = orig.ident2)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = c(donor_colors[1:5], Pt10Rel.Donor = "#da70d6", Pt10Rel.Host = "#fa8072")) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  ggtitle(paste0(ndim, " dimensions")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank())

print(plot_grid(p1, p2, p3, ncol = 3))
#}
dev.off()


# CDKN2A and RAB9A mutations in host/donor pDCs ---------------------------------------------------

# Determine the proportion of cells with CDKN2A deletions and RAB9a mutations in host vs. donor cells
metadata_tib <- as_tibble(NormalRelapse@meta.data, rownames = "cell")
metadata_tib %>% filter(orig.ident == "Pt10Rel") %>% .$orig.host.donor %>% table(useNA = "always")

# Visualize
pdf("10.3_Donor-host_mutations.pdf", width = 6, height = 4)

p1 <- metadata_tib %>% filter(orig.ident == "Pt10Rel") %>%
  group_by(orig.host.donor, MTAP.rearr) %>% summarize(n = n()) %>% ungroup %>%
  group_by(orig.host.donor) %>% mutate(Percent = n/sum(n)*100) %>%
  mutate(MTAP.rearr = factor(MTAP.rearr, levels = c("no call" , "mutant", "wildtype"))) %>%
  ggplot(aes(x = orig.host.donor, y = Percent, fill = MTAP.rearr)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = mut_colors) +
  ylab("Percentage of cells") +
  theme_bw() +
  theme(aspect.ratio = 3, panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_line(color = "black"))
  
p2 <- metadata_tib %>% filter(orig.ident == "Pt10Rel") %>%
  group_by(orig.host.donor, RAB9A.3pUTR) %>% summarize(n = n()) %>% ungroup %>%
  group_by(orig.host.donor) %>% mutate(Percent = n/sum(n)*100) %>%
  mutate(RAB9A.3pUTR = factor(RAB9A.3pUTR, levels = c("no call" , "mutant", "wildtype"))) %>%
  ggplot(aes(x = orig.host.donor, y = Percent, fill = RAB9A.3pUTR)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = mut_colors) +
  ylab("Percentage of cells") +
  theme_bw() +
  theme(aspect.ratio = 3, panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_line(color = "black"))

print(plot_grid(p1, p2, ncol = 2))

dev.off()
