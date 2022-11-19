# Peter van Galen, 220703
# Number of cells, UMIs, etc

library(tidyverse)
library(Seurat)

setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/05_Stats")

rm(list=ls())

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")

# Load Seurat objects
seurat_files <- list.files("../04_XV-seq", pattern = "*.rds", full.names = T)
seu_ls <- lapply(seurat_files, function(x) readRDS(x))
names(seu_ls) <- cutf(basename(seurat_files), d = "_")
seu <- merge(seu_ls[[1]], seu_ls[2:length(seu_ls)])

# Load genotyping information
genotyping_tables.tib <- read_excel("../04_XV-seq/XV-seq_overview.xlsx")
# Replace different MTAP entries with one, just as in 4.1_Add_GoT-XV_to_Seurat.R
genotyping_tables.tib$Mutation <- gsub("MTAP.rearr.*", "MTAP.rearr", genotyping_tables.tib$Mutation)


# Stats for Supplemental Table 2 ------------------------------------------------------------------

# Split bm, merge all
seu_ls$BM$orig.ident2 <- cutf(seu_ls$BM$replicate, d = "\\.")
seu_all_ls <- c(SplitObject(seu_ls$BM, split.by = "orig.ident2"), seu_ls[-1])

# Create stats table
stats_df <- data.frame(Sample_ID = names(seu_all_ls),
                       Cell_number = unlist(lapply(seu_all_ls, ncol)),
                       UMIs_mean = unlist(lapply(seu_all_ls, function(x) mean(x$nCount_RNA))),
                       UMIs_min = unlist(lapply(seu_all_ls, function(x) min(x$nCount_RNA))),
                       UMIs_max = unlist(lapply(seu_all_ls, function(x) max(x$nCount_RNA))),
                       Genes_mean = unlist(lapply(seu_all_ls, function(x) mean(x$nFeature_RNA))),
                       Genes_min = unlist(lapply(seu_all_ls, function(x) min(x$nFeature_RNA))),
                       Genes_max = unlist(lapply(seu_all_ls, function(x) max(x$nFeature_RNA))))
stats_order.df <- stats_df[c("BM1", "BM2", "BM3", "BM4", "BM5", "BM6", "Pt1Dx", "Pt1Rem", "Pt5Dx", "Pt9Dx",
                             "Pt10Dx", "Pt10Rel", "Pt12Dx", "Pt12Rel", "Pt14Dx", "Pt15Dx", "Pt16Dx"),]
#write.table(stats_order.df, file = "5.2_Stats.txt", row.names = F, sep = "\t", quote = F)


# Stats for Response to Reviewers -----------------------------------------------------------------

# "Of these data, 27,994 cells were genotyped at one or more of 40 individual DNA sites by XV-seq."
as_tibble(seu@meta.data) %>% select(all_of(unique(genotyping_tables.tib$Mutation))) %>% 
  apply(., 2, unique) %>% unlist %>% unique
as_tibble(seu@meta.data) %>% select(all_of(unique(genotyping_tables.tib$Mutation))) %>%
  apply(., 1, function(x) sum(grepl("wildtype|mutant", x)) > 0) %>% sum
length(unique(genotyping_tables.tib$Mutation))
genotyping_tables.tib %>% filter(`Initial submission or revision` == "Submission") %>% .$Mutation %>% unique %>% length

# "Our mean genotyping efficiency by XV-seq was 10.2% (range 0.1-99.0%)"
mean( na.omit(genotyping_tables.tib$`Genotyping efficiency (%)` ))
#median( na.omit(genotyping_tables.tib$`Genotyping efficiency (%)` ))
min( na.omit(genotyping_tables.tib$`Genotyping efficiency (%)` ))
max( na.omit(genotyping_tables.tib$`Genotyping efficiency (%)` ))











