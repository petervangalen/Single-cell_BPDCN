# Peter van Galen, 220707
# Compare genotyping efficiency of different approaches

library(tidyverse)
library(Seurat)
library(readxl)

setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/5_Stats")

rm(list=ls())

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")

# Load genotyping information
genotyping_tables.tib <- read_excel("../4_XV-seq/XV-seq_overview.xlsx")
# Replace different MTAP entries with one, just as in 4.1_Add_GoT-XV_to_Seurat.R
genotyping_tables.tib$Mutation <- gsub("MTAP.rearr.*", "MTAP.rearr", genotyping_tables.tib$Mutation)

# Load and subset genotyping efficiency data
subset_tib <- genotyping_tables.tib %>% dplyr::select(Mutation, `Genotyping efficiency (%)`, `Mutation identification`) %>%
  na.omit() %>% mutate(`Mutation identification` = gsub("WES|WGS", "WES/WGS", `Mutation identification`))

# What's the mean?
subset_tib %>% group_by(`Mutation identification`) %>% summarize(mean_efficiency = mean(`Genotyping efficiency (%)`))
subset_tib %>% group_by(`Mutation identification`) %>% summarize(mean_efficiency = median(`Genotyping efficiency (%)`))

pdf("5.4_Genotyping_efficiency.pdf")
subset_tib %>% 
  ggplot(aes(x = `Mutation identification`, y = `Genotyping efficiency (%)`, label = Mutation)) +
  geom_boxplot() +
  geom_point() +
  geom_label(alpha = 0.8, size = 3) +
  scale_y_continuous(trans = "log10") +
  theme(aspect.ratio = 1)
dev.off()

