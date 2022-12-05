# Peter van Galen, 221119
# Evaluate the overlap between allelic structures and mutation data of single cells in the BPDCN628 / Patient 9.

# Here's how I explained it to Andy on 221119:
# "We know haplotype B is lost at some point (Volker determined this from the skin tumor exome data). However, we 
# find co-occurrence of haplotype B and the indicated mutations. Therefore, we know the mutations occurred prior
# to full transformation (which is demarcated by loss of haplotype B). Thus, the mutations were present in the
# pre-malignant clone."

library(tidyverse)
library(readxl)
library(Seurat)

# Set working directory
setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/09_Pt9_haplotype")

rm(list=ls())

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")

# Load Seurat objects
seurat_files <- list.files("../04_XV-seq", pattern = "*.rds", full.names = T)
seu_ls <- lapply(seurat_files, function(x) readRDS(x))
names(seu_ls) <- cutf(basename(seurat_files), d = "_")

# Load genotyping information
genotyping_tables.tib <- read_excel("../04_XV-seq/XV-seq_overview.xlsx")
# Replace different MTAP entries with one, just as in 4.1_Add_GoT-XV_to_Seurat.R
genotyping_tables.tib$Mutation <- gsub("MTAP.rearr.*", "MTAP.rearr", genotyping_tables.tib$Mutation)

# From the exome data, Volker determined that haplotype B is lost in the skin tumor. Building on this,
# he quantified the detection of haplotype A and B in the single-cell data. Load his assessment:
chr7q <- read_tsv("628_7q_detected.txt")
plot(rev(sort(rowSums(chr7q[,2:3]))), ylim = c(0,12), ylab = "Cumulative coverage")

# Merge XV-seq and 7q data
merge_tib <- as_tibble(seu_ls[["Pt9Dx"]]@meta.data, rownames = "cell") %>%
  mutate(cell = cutf(cell, d = "-")) %>%
  select(cell, project.umap.x, project.umap.y, unique(filter(genotyping_tables.tib, Sample == "Pt9Dx"))$Mutation) %>%
  left_join(chr7q) %>%
  mutate(covA = replace_na(covA, 0), covB = replace_na(covB, 0)) # Add 0 coverage to cells that were not in Volker's table

# Plot coverage of both haplotype A and haplotype B in each single cell. Unfortunately, the coverage is insufficient to do haplotype calling of the single cells.
merge_tib %>%
  arrange(desc(covA), covB) %>% mutate(key = row_number()) %>%
  pivot_longer(col = c(covA, covB), names_to = "allele", values_to = "cov") %>%
  ggplot(aes(x = key, y = cov, color = allele)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio = 1)

# Table of allele coverage and mutation detection. This underlies the co-occurrence plot in Extended Data Figure 6.
merge_tib %>%
  mutate(hapA_detected = ifelse(covA > 0, yes = "A-covered", no = "A-not_detected"),
         hapB_detected = ifelse(covB > 0, yes = "B-covered", no = "B-not_detected")) %>%
  group_by(hapA_detected, hapB_detected) %>%
  summarize(n_cell = n(),
            `CUX1.L911fs*` = sum(`CUX1.L911fs*` == "mutant"),
            `TET2.E1437fs*` = sum(`TET2.E1437fs*` == "mutant"),
            `TET2.Q1547*` = sum(`TET2.Q1547*` == "mutant"))

# Because the mutations co-occur with haplotype B, we know they occured prior to full transformation (which is demarcated
# by loss of haplotype B)
