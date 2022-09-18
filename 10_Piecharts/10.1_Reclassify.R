# Peter van Galen, 220718
# Exclude cells that are suspect of misclassification

library(tidyverse)
library(Seurat)
library(readxl)
library(data.table)
library(randomForest)
library(gplots)
library(ggforce)
library(cowplot)

rm(list=ls())

# Set working directory
setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/10_Hierarchy")

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")
cell_colors <- popcol.tib$hex[1:21]
names(cell_colors) <- popcol.tib$pop[1:21]
donor_colors <- popcol.tib$hex[23:40]
names(donor_colors) <- popcol.tib$pop[23:40]
group_colors <- popcol.tib$hex[41:43]
names(group_colors) <- popcol.tib$pop[41:43]
mut_colors <- popcol.tib$hex[44:46]
names(mut_colors) <- popcol.tib$pop[44:46]

# Load Seurat objects
seurat_files <- list.files("../04_XV-seq", pattern = "*.rds", full.names = T)
seu_ls <- lapply(seurat_files, function(x) readRDS(x))
names(seu_ls) <- cutf(basename(seurat_files), d = "_")
seu <- merge(seu_ls[[1]], seu_ls[2:length(seu_ls)])
seu$orig.ident2 <- ifelse(grepl("BM", seu$orig.ident), yes = cutf(seu$replicate, d = "\\."), no = seu$orig.ident)
seu$CellType <- factor(seu$CellType, levels = levels(seu_ls[[1]]$CellType))

# Load genotyping information
genotyping_tables.tib <- read_excel("../04_XV-seq/XV-seq_overview.xlsx")
# Replace different MTAP entries with one, just as in 4.1_Add_GoT-XV_to_Seurat.R
genotyping_tables.tib$Mutation <- gsub("MTAP.rearr.*", "MTAP.rearr", genotyping_tables.tib$Mutation)

# Add malignant BPDCN cell module score, extract metadata
bpdcn_sign <- read.table("../09_pDC_expr/bpdcn_sign.txt")[,1]
seu <- AddModuleScore(seu, features = list(bpdcn_sign), name = "bpdcn_score")
metadata_tib <- as_tibble(seu@meta.data, rownames = "cell") %>%
  dplyr::select(cell, orig.ident, project.umap.x, project.umap.y, genotyping_tables.tib$Mutation,
                CellType, bm_involvement, bpdcn_score1) %>%
  rename(bpdcn_score = "bpdcn_score1")


# Add genotyping information ----------------------------------------------------------------------

# Add which cells have founder mutations
f_muts <- filter(genotyping_tables.tib, `Founder or progression mutation` == "Founder") %>% .$Mutation
metadata_tib$f_call <- ifelse(apply(dplyr::select(metadata_tib, all_of(f_muts)), 1, 
                                    function(x) sum(grepl("wildtype", x) > 0)), yes = "wildtype", no = "no call")
metadata_tib$f_call <- ifelse(apply(dplyr::select(metadata_tib, all_of(f_muts)), 1,
                                    function(x) sum(grepl("mutant", x) > 0)), yes = "mutant", no = metadata_tib$f_call)
metadata_tib$f_call <- factor(metadata_tib$f_call, levels = c("no call", "wildtype", "mutant"))
# Add which cells have progession mutations
p_muts <- filter(genotyping_tables.tib, `Founder or progression mutation` == "Progression") %>% .$Mutation
metadata_tib$p_call <- ifelse(apply(dplyr::select(metadata_tib, all_of(p_muts)), 1, 
                                    function(x) sum(grepl("wildtype", x) > 0)), yes = "wildtype", no = "no call")
metadata_tib$p_call <- ifelse(apply(dplyr::select(metadata_tib, all_of(p_muts)), 1,
                                    function(x) sum(grepl("mutant", x) > 0)), yes = "mutant", no = metadata_tib$p_call)
# Add UV-associated TC>TT mutations
prog_tctt <- filter(genotyping_tables.tib, `Founder or progression mutation` == "Progression", `TC>TT (UV-associated)` == "Yes") %>%
  .$Mutation %>% unique
metadata_tib$p_call <- ifelse(apply(dplyr::select(metadata_tib, all_of(prog_tctt)), 1,
                                    function(x) sum(grepl("mutant", x) > 0)), yes = "TC>TT", no = metadata_tib$p_call)
metadata_tib$p_call <- factor(metadata_tib$p_call, levels = c("no call", "wildtype", "mutant", "TC>TT"))


# Exclude non-pDCs with a high BPDCN score (likely misclassified BPDCN cells) ---------------------

# Visualize BPDCN score, split by patient bone marrow involvement
plot_tib <- metadata_tib[sample(nrow(metadata_tib)),]

pdf(file = "10.1.1_Scores_by_BM_involvement.pdf", width = 6, height = 5)
plot_tib %>%
  ggplot(aes(x = bm_involvement, y = bpdcn_score)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  geom_sina(aes(color = CellType, group = bm_involvement), size = 0.5) +
  scale_color_manual(values = cell_colors) +
  theme_bw() +
  theme(aspect.ratio = 1, panel.grid = element_blank(), axis.text = element_text(color = "black"))
dev.off()

# Exclude cells that are likely transformed BPDCN cells misclassified as something else
metadata_tib$exclude <- ifelse(metadata_tib$CellType != "pDC" & metadata_tib$bpdcn_score > 0.5, yes = T, no = F)

# Visualize mutations and cells to exclude
pdf(file = "10.1.2_Exclude.pdf", width = 12, height = 5)

p1 <- metadata_tib %>% arrange(p_call) %>%
  ggplot(aes(x = CellType, y = bpdcn_score)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  geom_sina(aes(color = p_call, group = CellType), size = 0.3) +
  scale_color_manual(values = c(mut_colors, `TC>TT` = "#daa520")) +
  theme_bw() +
  ggtitle("Colored by progression mutation status") +
  theme(aspect.ratio = 0.5, panel.grid = element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank(), axis.text = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5))

p2 <- metadata_tib %>%
  ggplot(aes(x = CellType, y = bpdcn_score, color = exclude)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  geom_sina(aes(color = exclude, group = CellType), size = 0.3) +
  scale_color_manual(values = c("#bdb76b", "#b22222")) +
  theme_bw() +
  ggtitle("Cells to exclude (suspect misclassification)") +
  theme(aspect.ratio = 0.5, panel.grid = element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank(), axis.text = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5))

plot_grid(p1, p2)

dev.off()

# Check accuracy of exclusion by evaluating expression of canonical marker genes ------------------

# Add gene epxression to metadata
all(metadata_tib$cell == colnames(seu))
genes <- c("CD19", "IL3RA", "CD4", "SDC1")
for (g in genes) {
  metadata_tib[,g] <- GetAssayData(seu, slot = "data")[g,]
}

pdf(file = "10.1.3_MarkerGenes.pdf", width = 5, height = 4)

p1 <- metadata_tib %>% filter(CellType == "ProB") %>%
  ggplot(aes(x = bpdcn_score > 0.5, y = CD19)) +
  geom_violin(scale = "width", fill = "#eee8aa") +
  geom_jitter(alpha = 0.4, size = 0.7) +
  ggtitle("ProB cells") +
  theme_bw() +
  theme(aspect.ratio = 2, panel.grid = element_blank(), axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"), plot.title = element_text(hjust = 0.5))

p2 <- metadata_tib %>% filter(CellType == "Plasma") %>%
  ggplot(aes(x = bpdcn_score > 0.5, y = SDC1)) +
  geom_violin(scale = "width", fill = "#eee8aa") +
  geom_jitter(alpha = 0.4, size = 0.7) +
  ggtitle("Plasma cells") +
  theme_bw() +
  theme(aspect.ratio = 2, panel.grid = element_blank(), axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"), plot.title = element_text(hjust = 0.5))

plot_grid(p1, p2)

dev.off()

# Save
metadata_subset_tib <- filter(metadata_tib, exclude == F)
metadata_subset_tib$exclude <- NULL
metadata_subset_tib$CD19 <- NULL
metadata_subset_tib$IL3RA <- NULL
metadata_subset_tib$CD4 <- NULL
metadata_subset_tib$SDC1 <- NULL
write_tsv(metadata_subset_tib, file = "subset_metadata.txt")





