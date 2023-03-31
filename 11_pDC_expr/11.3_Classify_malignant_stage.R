# Peter van Galen, 221130
# Classify malignant BPDCN cells

library(tidyverse)
library(janitor)
library(Seurat)
library(readxl)
library(ggforce)
library(cowplot)
library(viridis)

rm(list=ls())

# Set working directory
setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/11_pDC_expr")

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")
cell_colors <- popcol.tib$hex[1:21]
names(cell_colors) <- popcol.tib$pop[1:21]
mut_colors <- popcol.tib$hex[44:46]
names(mut_colors) <- popcol.tib$pop[44:46]
group_colors <- popcol.tib$hex[41:43]
names(group_colors) <- popcol.tib$pop[41:43]
malignant_stage <- popcol.tib$hex[47:50]
names(malignant_stage) <- popcol.tib$pop[47:50]

# Load Seurat objects
seurat_files <- list.files("../04_XV-seq", pattern = "*.rds", full.names = T)
seu_ls <- lapply(seurat_files, function(x) readRDS(x))
names(seu_ls) <- cutf(basename(seurat_files), d = "_")
seu <- merge(seu_ls[[1]], seu_ls[2:length(seu_ls)])
seu$CellType <- factor(seu$CellType, levels = levels(seu_ls[[1]]$CellType))

# Load genotyping information
genotyping_tables.tib <- read_excel("../04_XV-seq/XV-seq_overview.xlsx")
# Replace different MTAP entries with one, just as in 4.1_Add_GoT-XV_to_Seurat.R
genotyping_tables.tib$Mutation <- gsub("MTAP.rearr.*", "MTAP.rearr", genotyping_tables.tib$Mutation)


# Calculate BPDCN Signature Score -----------------------------------------------------------------

# Score BPDCN signature that was made in the previous script
bpdcn_sign <- read.table("bpdcn_sign.txt")[,1]
seu <- AddModuleScore(seu, features = list(bpdcn_sign), name = "bpdcn_sign_score")
colnames(seu@meta.data) <- gsub("score1$", "score", colnames(seu@meta.data))


# pDC UMAP ----------------------------------------------------------------------------------------

# Load pDC gene expression and metadata (including UMAP coordinates) from 11.1_pDC_object.R
pdcs <- readRDS("pdcs.rds")
# Add signature score to pdcs
pdcs$bpdcn_sign_score <- seu$bpdcn_sign_score[colnames(pdcs)]

# Visualize
pdf("11.3.1_ModuleScores.pdf", width = 6, height = 6)
FeaturePlot(pdcs, features = "bpdcn_sign_score", pt.size = 0.6, order = T) +
  scale_color_viridis() +
  theme_void() +
  theme(aspect.ratio = 1, panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(hjust = 0.5))
dev.off()


# Relate BPDCN Signature Score to bone marrow involvement -----------------------------------------

pdf(file = "11.3.2_Scores_by_BM_involvement.pdf", width = 6, height = 5)
seu@meta.data %>% mutate(my_order = sample(ncol(seu))) %>% arrange(my_order) %>%
  ggplot(aes(x = bm_involvement, y = bpdcn_sign_score)) +
  geom_sina(aes(color = CellType, group = bm_involvement), size = 0.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  scale_color_manual(values = cell_colors) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  theme_bw() +
  theme(aspect.ratio = 1, panel.grid = element_blank(), axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"))
dev.off()


# Add extra MALAT1 genoytping information ----------------------------------------------------
# See also 11.4_pDC_Heatmaps.R

# I'm updating MALAT1 genotyping information with cells from transcripts that were supported by only 1 read. This is because
# of the outsized importance of the MALAT1 mutation, and the low sequencing depth we achieved. 
malat1_ReadThreshold1 <- read_tsv("../04_XV-seq/Pt12Dx/MALAT1.5155/Patient12_MALAT1.5155_summTable_ReadThreshold1.txt")
wt_calls <- colnames(subset(pdcs, orig.ident == "Pt12Dx"))[cutf(colnames(subset(pdcs, orig.ident == "Pt12Dx")),
                                                                d = "-") %in% filter(malat1_ReadThreshold1, mutUMIs == 0)$BC]
mut_calls <- colnames(subset(pdcs, orig.ident == "Pt12Dx"))[cutf(colnames(subset(pdcs, orig.ident == "Pt12Dx")),
                                                                 d = "-") %in% filter(malat1_ReadThreshold1, mutUMIs > 0)$BC]
pdcs$progression <- ifelse(colnames(pdcs) %in% wt_calls, yes = "wildtype", no = pdcs$progression)
pdcs$progression <- ifelse(colnames(pdcs) %in% mut_calls, yes = "mutant", no = pdcs$progression)

# Calculate P-value (tidy to base R)
cont_tib <- as_tibble(pdcs@meta.data) %>% filter(bm_involvement == "No") %>%
  mutate(progression = gsub("no call|wildtype", "not.mutant", progression)) %>%
  group_by(bpdcn_sign_score > 0.5, progression) %>% count %>% ungroup() %>% pivot_wider(names_from = "progression", values_from = "n")
cont_df <- data.frame(cont_tib[,2:3], row.names = c("low_score", "high_score"))
stats::chisq.test(cont_df)
# Calculate P-value (tidy)
as_tibble(pdcs@meta.data) %>% filter(bm_involvement == "No") %>%
  mutate(progression = gsub("no call|wildtype", "not.mutant", progression), score = ifelse(bpdcn_sign_score > 0.5, yes = "high", no = "low")) %>%
  tabyl(progression, score) %>%
  chisq.test()


# Relate BPDCN Signature Score to progression mutations in samples without involvement ------------

# The following plots justify why pDCs from patients without marrow involvement with a BPDCN Signature
# Score of >0.5 should be classified as malignant

as_tibble(pdcs@meta.data) %>% tabyl(orig.ident, bm_involvement)

pdf("11.3.3_Non-involved_pDCs.pdf", width = 5, height = 5)

# Founder mutations
as_tibble(pdcs@meta.data) %>% filter(bm_involvement != "Yes") %>%
  mutate(founder = factor(founder, levels = c("no call", "mutant", "wildtype"))) %>%
  arrange(founder) %>%
  ggplot(aes(x = orig.ident, y = bpdcn_sign_score, color = founder)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  geom_sina(aes(group = orig.ident)) +
  scale_color_manual(values = mut_colors) +
  coord_cartesian(ylim = c(min(seu$bpdcn_sign_score), max(seu$bpdcn_sign_score))) +
  theme_bw() +
  theme(aspect.ratio = 1, panel.grid = element_blank(), axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_line(color = "black"))

# Progression mutations
as_tibble(pdcs@meta.data) %>% filter(bm_involvement != "Yes") %>%
  mutate(progression = factor(progression, levels = c("no call", "mutant", "wildtype"))) %>%
  arrange(progression) %>%
  ggplot(aes(x = orig.ident, y = bpdcn_sign_score, color = progression)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  geom_sina(aes(group = orig.ident)) +
  scale_color_manual(values = mut_colors) +
  coord_cartesian(ylim = c(min(seu$bpdcn_sign_score), max(seu$bpdcn_sign_score))) +
  theme_bw() +
  theme(aspect.ratio = 1, panel.grid = element_blank(), axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_line(color = "black"))

dev.off()

# Main text: "While the clinical significance of this finding requires further investigation, XV-seq revealed progression mutations in 19/23 (82.3%) of these cells but only 1/495 other pDCs, verifying their identity as malignant BPDCN cells (P<2.2E-16, Extended Data Fig. 8c-e)."
as_tibble(pdcs@meta.data) %>% filter(bm_involvement == "No") %>% mutate(bpdcn_exceed = bpdcn_sign_score >= 0.5) %>%
  tabyl(bpdcn_exceed, progression) %>% adorn_totals(where = "col")
contingency_table <- data.frame(mutant = c(1, 19), not_mutant = c(494, 4), row.names = c("low_score", "high_score"))
chisq.test(contingency_table)
print(chisq.test(contingency_table))$p.value


# Canonical cell type markers in non-pDCs with high BPDCN Signature Score  ------------------------

# This justifies why cells not classified as pDCs in patients with marrow involvement but with a
# BPDCN Signature Score of >0.5 should be classified as malignant.

metadata_tib <- as_tibble(seu@meta.data, rownames = "cell")

all(metadata_tib$cell == colnames(seu))
genes <- c("CD19", "IL3RA", "CD4", "SDC1")
for (g in genes) {
  metadata_tib[,g] <- GetAssayData(seu, slot = "data")[g,]
}

pdf(file = "11.3.4_MarkerGenes.pdf", width = 5, height = 4)

p1 <- metadata_tib %>% filter(CellType == "ProB") %>%
  ggplot(aes(x = bpdcn_sign_score > 0.5, y = CD19)) +
  geom_violin(scale = "width", fill = "#eee8aa") +
  geom_jitter(alpha = 0.4, size = 0.7) +
  ggtitle("ProB cells") +
  theme_bw() +
  theme(aspect.ratio = 2, panel.grid = element_blank(), axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"), plot.title = element_text(hjust = 0.5))

p2 <- metadata_tib %>% filter(CellType == "Plasma") %>%
  ggplot(aes(x = bpdcn_sign_score > 0.5, y = SDC1)) +
  geom_violin(scale = "width", fill = "#eee8aa") +
  geom_jitter(alpha = 0.4, size = 0.7) +
  ggtitle("Plasma cells") +
  theme_bw() +
  theme(aspect.ratio = 2, panel.grid = element_blank(), axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"), plot.title = element_text(hjust = 0.5))

plot_grid(p1, p2)

dev.off()


# Classify healthy / premalignant / malignant pDCs ------------------------------------------------

# Based on the plots above and downstream analysis, classify cells as healthy/premalignant/malignant cells:
# Healthy donor samples, pDCs --> Healthy pDCs
# Samples without BM involvement, pDCs with BPDCN signature score < 0.5 --> Premalignant pDCs
# Samples without BM involvement, pDCs with BPDCN signature score > 0.5 --> Malignant BPDCN cells
# Samples with BM involvement, pDCs --> Malignant BPDCN cells
# Samples with BM involvement, any cell with BPDCN signature score > 0.5 --> Malignant BPDCN cells
seu$is_malignant <- ifelse(seu$bm_involvement == "HD" & seu$CellType == "pDC", yes = "Healthy", no = "Other")
seu$is_malignant <- ifelse(seu$bm_involvement == "No" & seu$CellType == "pDC" & seu$bpdcn_sign_score < 0.5,
                           yes = "Premalignant", no = seu$is_malignant)
seu$is_malignant <- ifelse(seu$bm_involvement == "Yes" & seu$CellType == "pDC", yes = "Malignant", no = seu$is_malignant)
seu$is_malignant <- ifelse(seu$bpdcn_sign_score > 0.5, yes = "Malignant", no = seu$is_malignant)
seu$is_malignant <- factor(seu$is_malignant, levels = c("Healthy", "Premalignant", "Malignant", "Other"))

# Check / numbers for EDF7
as_tibble(seu@meta.data) %>% filter(CellType == "pDC") %>% count(orig.ident, CellType, bm_involvement, is_malignant)
as_tibble(seu@meta.data) %>% .$is_malignant %>% table(useNA = "always")
as_tibble(seu@meta.data) %>% filter(CellType == "pDC", bm_involvement == "No") %>% .$is_malignant %>% table
as_tibble(seu@meta.data) %>% filter(CellType == "pDC", bm_involvement == "Yes") %>% .$is_malignant %>% table
as_tibble(seu@meta.data) %>% filter(CellType != "pDC", bm_involvement == "Yes") %>% .$is_malignant %>% table

pdf(file = "11.3.5_Classify_Premalignant_Malignant.pdf", width = 12, height = 5)

as_tibble(seu@meta.data) %>% 
  arrange(is_malignant) %>%
  mutate(my_order = c(sample((sum(grepl("Other", seu@meta.data$is_malignant))+1):ncol(seu)),
         sample(1:sum(grepl("Other", seu@meta.data$is_malignant))))) %>%
  arrange(my_order) %>%
  ggplot(aes(x = CellType, y = bpdcn_sign_score, color = is_malignant)) +
  geom_sina(aes(group = CellType), size = 0.5) +
  scale_color_manual(values = malignant_stage) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  ggtitle(paste("Classification of all", ncol(seu), "cells")) +
  theme_bw() +
  theme(aspect.ratio = 0.5, panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(), axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5))

# Check by coloring cells according to the presence of progression mutations
p_mut <- filter(genotyping_tables.tib, `Founder or progression mutation` == "Progression") %>% .$Mutation %>% unique
metadata_tib$progression <- ifelse(apply(metadata_tib[,p_mut], 1, function(x) sum(grepl("wildtype", x) > 0)), yes = "wildtype", no = "no call")
metadata_tib$progression <- ifelse(apply(metadata_tib[,p_mut], 1, function(x) sum(grepl("mutant", x) > 0)), yes = "mutant", no = metadata_tib$progression)

metadata_tib %>% filter(progression %in% c("wildtype", "mutant")) %>%
  ggplot(aes(x = CellType, y = bpdcn_sign_score, color = progression)) +
  geom_sina(aes(group = CellType), size = 0.5) +
  scale_color_manual(values = mut_colors[-3]) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  ggtitle(paste("Classification of all", ncol(seu), "cells")) +
  theme_bw() +
  theme(aspect.ratio = 0.5, panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(), axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5))

dev.off()


# Some stats and save -----------------------------------------------------------------------------

seu@meta.data %>% filter(CellType == "ProB") %>% mutate(bpdcn_exceed = bpdcn_sign_score >= 0.5) %>%
  group_by(CellType, is_malignant, bpdcn_exceed) %>%
  summarize(n = n(), min_s = min(bpdcn_sign_score), max_s = max(bpdcn_sign_score))
seu@meta.data %>% filter(CellType == "Plasma") %>% mutate(bpdcn_exceed = bpdcn_sign_score >= 0.5) %>%
  group_by(CellType, is_malignant, bpdcn_exceed) %>%
  summarize(n = n(), min_s = min(bpdcn_sign_score), max_s = max(bpdcn_sign_score))

# Retrieve random forest prediction scores to plot on x-axis in later analyses
prediction_scores_ls <- readRDS("../03_RandomForest/Prediction_scores.rds")
pDC_scores_tib <- as_tibble(do.call(rbind, lapply(prediction_scores_ls, function(x) data.frame(x[,"pDC"]))), rownames = "cell")
pDC_scores_tib <- pDC_scores_tib %>% rename(RF_pDC_score = "x....pDC..") %>% mutate(cell = gsub("^.*?\\.", "", cell))

# Merge & save
metadata2_tib <- as_tibble(seu@meta.data, rownames = "cell") %>% left_join(pDC_scores_tib) %>%
  dplyr::select(cell, orig.ident, CellType, bpdcn_sign_score, RF_pDC_score, is_malignant)
# Quickly visualize (don't include normal cells; they don't have an RF_pDC_score)
metadata2_tib %>% filter(orig.ident != "BM") %>%
  ggplot(aes(x = RF_pDC_score, y = bpdcn_sign_score, color = is_malignant)) +
  geom_point() +
  scale_color_manual(values = malignant_stage) +
  theme_bw() +
  theme(aspect.ratio = 1, panel.grid = element_blank(),
        axis.ticks = element_line(color = "black"), axis.text = element_text(color="black"))
# It looks like a few healthy/premalignant pDCs from patients with bone marrow involvement are erroneously classified as malignant. This is a limitation of the study that does not affect overall conclusions.

write_tsv(metadata2_tib, file = "11.3_MalignantCalls_Final.txt")

