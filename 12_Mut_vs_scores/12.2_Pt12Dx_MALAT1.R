# Peter van Galen, 221130
# Visualize the occurrence of progression mutations in Patient 12 diagnosis cells
# This is similar to 12.1_Pt12Dx_pDCs but simplified and focused

library(tidyverse)
library(Seurat)
library(readxl)
library(ggforce)

rm(list=ls())

# Set working directory
setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/12_Mut_vs_scores")

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")
mut_colors <- popcol.tib$hex[44:46]
names(mut_colors) <- popcol.tib$pop[44:46]

# Load all data
seurat_files <- list.files("../04_XV-seq", pattern = "*.rds", full.names = T)
seu_ls <- lapply(seurat_files, function(x) readRDS(x))
names(seu_ls) <- cutf(basename(seurat_files), d = "_")
seu <- merge(seu_ls[[1]], seu_ls[2:length(seu_ls)])

# Add metadata from 11_pDC_expr
MalignantCalls_df <- read.table("../11_pDC_expr/11.3_MalignantCalls_Final.txt", header = T, row.names = "cell")
all(rownames(seu@meta.data) == rownames(MalignantCalls_df))
seu$bpdcn_sign_score <- MalignantCalls_df$bpdcn_sign_score
seu$RF_pDC_score <- MalignantCalls_df$RF_pDC_score
seu$is_malignant <- MalignantCalls_df$is_malignant


# Include additional mutations for MALAT1 ---------------------------------------------------------
# Similar to 11.3_Classify_malignant_stage.R and 11.4_pDC_Heatmaps.R

Pt12Dx_metadata_tib <- as_tibble(subset(seu, orig.ident == "Pt12Dx")@meta.data, rownames = "cell")

malat1_ReadThreshold1 <- read_tsv("../04_XV-seq/Pt12Dx/MALAT1.5155/Patient12_MALAT1.5155_summTable_ReadThreshold1.txt")

wt_calls <- Pt12Dx_metadata_tib$cell[cutf(Pt12Dx_metadata_tib$cell, d = "-") %in% filter(malat1_ReadThreshold1, mutUMIs == 0)$BC]
mut_calls <- Pt12Dx_metadata_tib$cell[cutf(Pt12Dx_metadata_tib$cell, d = "-") %in% filter(malat1_ReadThreshold1, mutUMIs > 0)$BC]

Pt12Dx_metadata_tib$`MALAT1.chr11:65270399:G/A` <- ifelse(Pt12Dx_metadata_tib$cell %in% wt_calls, yes = "wildtype", no = Pt12Dx_metadata_tib$`MALAT1.chr11:65270399:G/A`)
Pt12Dx_metadata_tib$`MALAT1.chr11:65270399:G/A` <- ifelse(Pt12Dx_metadata_tib$cell %in% mut_calls, yes = "mutant", no = Pt12Dx_metadata_tib$`MALAT1.chr11:65270399:G/A`)


# Plot progression mutations vs. pDC scores & BPDCN scores ----------------------------------------
# pDC scores are from the random forest classifier (03_RandomForest)
# BPDCN scores are from differential gene expression analysis (11_pDC_expr)
# Mutations are from XV-seq

# Wrangle more (mainly reorderdering)
n_call <- sum(grepl("wildtype|mutant", Pt12Dx_metadata_tib$`MALAT1.chr11:65270399:G/A`))
Pt12Dx_order_tib <- Pt12Dx_metadata_tib %>%
  mutate(`MALAT1.chr11:65270399:G/A` = factor(`MALAT1.chr11:65270399:G/A`, levels = c("wildtype", "mutant", "no call"))) %>%
  arrange(`MALAT1.chr11:65270399:G/A`) %>%
  mutate(my_order = c(sample(n_call), (n_call+1):nrow(Pt12Dx_metadata_tib))) %>%
  arrange(desc(my_order))

pdf("12.2.1_Pt12Dx_Scores-vs-Muts.pdf", width = 5, height = 5)

Pt12Dx_order_tib %>%
  ggplot(aes(x = RF_pDC_score, y = bpdcn_sign_score, color = `MALAT1.chr11:65270399:G/A`)) +
  geom_point(size = 0.8) +
  scale_color_manual(values = mut_colors) +
  coord_cartesian(ylim = c(-0.5, 1.3)) +
  theme_bw() +
  theme(aspect.ratio = 1, panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"))

dev.off()

