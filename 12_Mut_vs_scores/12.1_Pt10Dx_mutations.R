# Peter van Galen, 221130
# Visualize the occurrence of progression mutations in pDCs from Patient 10 vs. BDPCN signature score and pDC prediction score

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
donor_colors <- popcol.tib$hex[23:40]
names(donor_colors) <- popcol.tib$pop[23:40]
group_colors <- popcol.tib$hex[41:43]
names(group_colors) <- popcol.tib$pop[41:43]
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

# Load genotyping information
genotyping_tables.tib <- read_excel("../04_XV-seq/XV-seq_overview.xlsx")
# Replace different MTAP entries with one, just as in 4.1_Add_GoT-XV_to_Seurat.R
genotyping_tables.tib$Mutation <- gsub("MTAP.rearr.*", "MTAP.rearr", genotyping_tables.tib$Mutation)


# Generate metadata tibble with mutation information ----------------------------------------------

# Define mutations of interest for this script
Pt10_muts <- filter(genotyping_tables.tib, Sample == "Pt10Dx") %>% arrange(`Founder or progression mutation`) %>%
  .$Mutation %>% unique
Pt10_F_mut <- filter(genotyping_tables.tib, `Founder or progression mutation` == "Founder", Sample == "Pt10Dx") %>%
  .$Mutation %>% unique
Pt10_P_mut <- filter(genotyping_tables.tib, `Founder or progression mutation` == "Progression", Sample == "Pt10Dx") %>%
  .$Mutation %>% unique # In the text and legend, we say n = 6, b/c DOLPP1 is uninformative

# Extract metadata
pt10dx_metadata_tib <- as_tibble(seu@meta.data, rownames = "cell") %>%
  filter(orig.ident == "Pt10Dx") %>%
  dplyr::select(CellType, bpdcn_sign_score, RF_pDC_score, all_of(Pt10_muts))

# Add summary columns for founder & progression mutations
pt10dx_metadata_tib <- pt10dx_metadata_tib %>% mutate(any_mut = ifelse(apply(pt10dx_metadata_tib[,Pt10_muts], 1, function(x)
  sum(grepl("wildtype", x) > 0)), yes = "wildtype", no = "no call"))
pt10dx_metadata_tib <- pt10dx_metadata_tib %>% mutate(any_mut = ifelse(apply(pt10dx_metadata_tib[,Pt10_muts], 1, function(x)
  sum(grepl("mutant", x) > 0)), yes = "mutant", no = .$any_mut))
pt10dx_metadata_tib <- pt10dx_metadata_tib %>% mutate(Founder = ifelse(apply(pt10dx_metadata_tib[,Pt10_F_mut], 1, function(x)
  sum(grepl("wildtype", x) > 0)), yes = "wildtype", no = "no call"))
pt10dx_metadata_tib <- pt10dx_metadata_tib %>% mutate(Founder = ifelse(apply(pt10dx_metadata_tib[,Pt10_F_mut], 1, function(x)
  sum(grepl("mutant", x) > 0)), yes = "mutant", no = .$Founder))
pt10dx_metadata_tib <- pt10dx_metadata_tib %>% mutate(Progression = ifelse(apply(pt10dx_metadata_tib[,Pt10_P_mut], 1, function(x)
  sum(grepl("wildtype", x) > 0)), yes = "wildtype", no = "no call"))
pt10dx_metadata_tib <- pt10dx_metadata_tib %>% mutate(Progression = ifelse(apply(pt10dx_metadata_tib[,Pt10_P_mut], 1, function(x)
  sum(grepl("mutant", x) > 0)), yes = "mutant", no = .$Progression))


# Violin plots of pDCs  with BPDCN Signature Scores for all mutations -----------------------------

# Wrangle
plot_tib <- pt10dx_metadata_tib %>% filter(CellType == "pDC") %>%
  pivot_longer(cols = all_of(Pt10_muts), names_to = "Mutation", values_to = "XVseq") %>%
  mutate(XVseq = factor(XVseq, levels = c("no call", "wildtype", "mutant")),
         Mutation = factor(Mutation, levels = unique(Pt10_muts)))
plot_tib <- plot_tib %>% arrange(Mutation, XVseq) %>%
  mutate(XVseq = gsub("wildtype", "#32cd32", gsub("mutant", "#dc143c", gsub("no call", "#dcdcdc", as.character(XVseq)))))

# Visualize Pt10Dx sina plots with mut/wt calls for every mutation
pdf("12.1.1_Pt10Dx_pDC_scores_XVseq.pdf", width = 10, height = 6)
plot_tib %>%
  ggplot(aes(x = Mutation, y = bpdcn_sign_score)) +
  geom_sina(color = plot_tib$XVseq) +
  geom_violin(fill = NA) +
  coord_cartesian(ylim = c(min(seu$bpdcn_sign_score), max(seu$bpdcn_sign_score))) +
  annotate(geom = "text", x = c(4, 12), y = c(1.3, 1.3), label = c("Founder", "Progression")) +
  annotate(geom = "segment", x = c(1, 9), xend = c(7, 15), y = c(1.2, 1.2), yend = c(1.2, 1.2)) +
  theme_bw() +
  theme(aspect.ratio = 0.6, panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"), axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# Wrangle to combine collapse violins from the previous visualization
plot2_tib <- pt10dx_metadata_tib %>% filter(CellType == "pDC") %>%
  pivot_longer(cols = c("Founder", "Progression"), names_to = "Mutation", values_to = "XVseq") %>%
  mutate(XVseq = factor(XVseq, levels = c("no call", "wildtype", "mutant")))
plot2_tib <- plot2_tib %>% arrange(Mutation, XVseq) %>%
  mutate(XVseq = gsub("wildtype", "#32cd32", gsub("mutant", "#dc143c", gsub("no call", "#dcdcdc", as.character(XVseq)))))

# Visualize
pdf("12.1.2_Pt10Dx_pDC_scores_combined_XVseq.pdf", width = 3, height = 4)
plot2_tib %>%
  ggplot(aes(x = Mutation, y = bpdcn_sign_score)) +
  geom_sina(color = plot2_tib$XVseq) +
  geom_violin(fill = NA) +
  coord_cartesian(ylim = c(min(seu$bpdcn_sign_score), max(seu$bpdcn_sign_score))) +
  theme_bw() +
  theme(aspect.ratio = 2, panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"), axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


# Plot progression mutations vs. pDC scores & BPDCN scores ----------------------------------------
# pDC scores are from the random forest classifier (03_RandomForest)
# BPDCN scores are from differential gene expression analysis (11_pDC_expr)
# Mutations are from XV-seq

# Wrangle (pivot)
pt10dx_mut_tib <- pt10dx_metadata_tib %>%
  pivot_longer(cols = all_of(Pt10_muts), names_to = "mut", values_to = "mut_call") #%>%
  #filter(mut != "DOLPP1.R227S") # this one's uninformative (no mutant calls at Dx)

# Wrangle more (mainly reorderdering)
n_call <- sum(grepl("wildtype|mutant", pt10dx_mut_tib$mut_call))
pt10dx_order_tib <- pt10dx_mut_tib %>%
  mutate(mut_call = factor(mut_call, levels = c("wildtype", "mutant", "no call"))) %>%
  arrange(mut_call) %>%
  mutate(my_order = c(sample(n_call), (n_call+1):nrow(pt10dx_mut_tib))) %>%
  arrange(desc(my_order)) %>%
  mutate(mut = factor(mut, levels = Pt10_muts))

pdf("12.1.3_Pt10Dx_Scores-vs-Muts.pdf", width = 8, height = 8)

pt10dx_order_tib %>%
  ggplot(aes(x = RF_pDC_score, y = bpdcn_sign_score, color = mut_call)) +
  geom_point(size = 0.8) +
  scale_color_manual(values = mut_colors) +
  coord_cartesian(ylim = c(-0.5, 1.3)) +
  theme_bw() +
  theme(aspect.ratio = 1, panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  facet_wrap(~ mut)

dev.off()

