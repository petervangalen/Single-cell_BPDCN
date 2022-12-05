# Peter van Galen, 221119
# Test the co-occurrence of mutations to generate Extended Data Figure 6 panel e and f

library(tidyverse)
library(readxl)
library(Seurat)

# Set working directory
setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/13_Misc")

rm(list=ls())

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")

# Load Seurat objects
seurat_files <- list.files("../04_XV-seq", pattern = "*.rds", full.names = T)
seu_ls <- lapply(seurat_files, function(x) readRDS(x))
names(seu_ls) <- cutf(basename(seurat_files), d = "_")
seu <- merge(seu_ls[[1]], seu_ls[2:length(seu_ls)])
seu$orig.ident2 <- ifelse(grepl("BM", seu$orig.ident), yes = cutf(seu$replicate, d = "\\."), no = seu$orig.ident)

# Load genotyping information
genotyping_tables.tib <- read_excel("../04_XV-seq/XV-seq_overview.xlsx")
# Replace different MTAP entries with one, just as in 4.1_Add_GoT-XV_to_Seurat.R
genotyping_tables.tib$Mutation <- gsub("MTAP.rearr.*", "MTAP.rearr", genotyping_tables.tib$Mutation)

# View tables
#pt <- "Pt1Dx"
#pt <- "Pt1Rem"
pt <- "Pt9Dx"
pt <- "Pt10Dx"
pt <- "Pt10Rel"
#pt <- "Pt12Dx"
#pt <- "Pt12Rel"
as_tibble(seu_ls[[pt]]@meta.data) %>%
  select(unique(filter(genotyping_tables.tib, Sample == pt)$Mutation)) %>%
  view()

# Calculate numbers for Pt10
mut_tib <- as_tibble(seu@meta.data) %>% filter(orig.ident %in% c("Pt10Dx", "Pt10Rel")) %>%
  select(unique(filter(genotyping_tables.tib, Sample %in% c("Pt10Dx", "Pt10Rel"))$Mutation))

mut_tib %>% filter(`TET2.S792*` == "mutant") %>% .$`TET2.Q1034*` %>% table
mut_tib %>% filter(`TET2.S792*` == "mutant") %>% .$`TET2.R1216*` %>% table
mut_tib %>% filter(`TET2.S792*` == "mutant") %>% .$`TET2.H1380Y` %>% table

mut_tib %>% filter(`TET2.H1380Y` == "mutant") %>% .$`TET2.S792*` %>% table
mut_tib %>% filter(`TET2.H1380Y` == "mutant") %>% .$`TET2.Q1034*` %>% table
mut_tib %>% filter(`TET2.H1380Y` == "mutant") %>% .$`TET2.R1216*` %>% table

mut_tib %>% filter(`ASXL1.G642fs` == "mutant") %>% .$`TET2.S792*` %>% table
mut_tib %>% filter(`ASXL1.G642fs` == "mutant") %>% .$`TET2.Q1034*` %>% table
mut_tib %>% filter(`ASXL1.G642fs` == "mutant") %>% .$`TET2.R1216*` %>% table
mut_tib %>% filter(`ASXL1.G642fs` == "mutant") %>% .$`TET2.H1380Y` %>% table

# Calculate numbers for Pt9
mut_tib <- as_tibble(seu@meta.data) %>% filter(orig.ident == "Pt9Dx") %>%
  select(unique(filter(genotyping_tables.tib, Sample == "Pt9Dx"))$Mutation)

mut_tib %>% filter(`TET2.E1437fs*` == "mutant") %>% .$`TET2.Q1547*` %>% table
mut_tib %>% filter(`TET2.E1437fs*` == "mutant") %>% .$`CUX1.L911fs*` %>% table
mut_tib %>% filter(`TET2.Q1547*` == "mutant") %>% .$`CUX1.L911fs*` %>% table






