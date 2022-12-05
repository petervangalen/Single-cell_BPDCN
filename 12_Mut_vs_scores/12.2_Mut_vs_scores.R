# Peter van Galen, 221112
# Relate pDCs gene expression scores and classification to their mutational status

library(tidyverse)
library(Seurat)
library(readxl)
library(ggrastr)

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
seu <- AddMetaData(seu, MalignantCalls_df[, c("bpdcn_sign_score", "RF_pDC_score", "is_malignant")])
seu$is_malignant <- factor(seu$is_malignant, levels = c("Healthy", "Premalignant", "Malignant"))

# Load genotyping information
genotyping_tables.tib <- read_excel("../04_XV-seq/XV-seq_overview.xlsx")
# Replace different MTAP entries with one, just as in 4.1_Add_GoT-XV_to_Seurat.R
genotyping_tables.tib$Mutation <- gsub("MTAP.rearr.*", "MTAP.rearr", genotyping_tables.tib$Mutation)

# Prepare
metadata_tib <- as_tibble(seu@meta.data, rownames = "cell")
muts <- unique(genotyping_tables.tib$Mutation)


# Visualize by mutation ---------------------------------------------------------------------------

# Create plot titles for facet_wrap
labels_tib <- genotyping_tables.tib %>% group_by(Mutation) %>%
  mutate(Samples = paste0(unique(Sample), collapse = ", "),
         CC.CT = ifelse(`CC>CT (UV-associated)` == "Yes", yes = ", CC>CT", no = ""),
         TC.TT = ifelse(`TC>TT (UV-associated)` == "Yes", yes = ", TC>TT", no = ""),
         CC.TT = ifelse(`CC>TT (UV-specific)` == "Yes", yes = ", CC>TT", no = "")) %>% #view
  summarize(Samples = Samples, Type = `Founder or progression mutation`[1], CC.CT = CC.CT[1], TC.TT = TC.TT[1], CC.TT = CC.TT[1]) %>%
  unique %>%
  mutate(title = paste0(Mutation, "\n", Samples, " (", Type, CC.CT, TC.TT, CC.TT, ")")) %>%
  select(Mutation, title) %>% rename(mut = "Mutation")

pdf("12.2.1_Scores-vs-Muts.pdf", width = 20, height = 20)
metadata_tib %>%
  select(RF_pDC_score, bpdcn_sign_score, all_of(muts)) %>%
  pivot_longer(cols = all_of(muts), names_to = "mut", values_to = "mut_call") %>%
  left_join(labels_tib) %>%
  mutate(mut_call = factor(mut_call, levels = c("no call", "wildtype", "mutant")),
         title = factor(title, levels = unique(title))) %>%
  arrange(mut_call) %>% na.omit %>%
  ggplot(aes(x = RF_pDC_score, y = bpdcn_sign_score, color = mut_call)) +
  #geom_point_rast() +
  geom_point() +
  scale_color_manual(values = mut_colors) +
  theme_bw() +
  theme(aspect.ratio = 1, panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  facet_wrap(~ title)
dev.off()


# By mutation type --------------------------------------------------------------------------------

# Define mutation types
f_mut <- filter(genotyping_tables.tib, `Founder or progression mutation` == "Founder") %>% .$Mutation %>% unique
p_mut <- filter(genotyping_tables.tib, `Founder or progression mutation` == "Progression") %>% .$Mutation %>% unique
cc.ct <- filter(genotyping_tables.tib, `Founder or progression mutation` != "Founder", `CC>CT (UV-associated)` == "Yes") %>%
  .$Mutation %>% unique
tc.tt <- filter(genotyping_tables.tib, `Founder or progression mutation` != "Founder", `TC>TT (UV-associated)` == "Yes") %>%
  .$Mutation %>% unique
cc.tt <- filter(genotyping_tables.tib, `Founder or progression mutation` != "Founder", `CC>TT (UV-specific)` == "Yes") %>%
  .$Mutation %>% unique

# Add metadata of mutation types to Seurat object
seu$f_call <- ifelse(apply(dplyr::select(seu@meta.data, all_of(f_mut)), 1,
                           function(x) sum(grepl("wildtype", x) > 0)), yes = "wildtype", no = "no call")
seu$f_call <- ifelse(apply(dplyr::select(seu@meta.data, all_of(f_mut)), 1,
                           function(x) sum(grepl("mutant", x) > 0)), yes = "mutant", no = seu$f_call)
seu$p_call <- ifelse(apply(dplyr::select(seu@meta.data, all_of(p_mut)), 1,
                           function(x) sum(grepl("wildtype", x) > 0)), yes = "wildtype", no = "no call")
seu$p_call <- ifelse(apply(dplyr::select(seu@meta.data, all_of(p_mut)), 1,
                           function(x) sum(grepl("mutant", x) > 0)), yes = "mutant", no = seu$p_call)
seu$cc.ct <- ifelse(apply(dplyr::select(seu@meta.data, all_of(cc.ct)), 1,
                          function(x) sum(grepl("wildtype", x) > 0)), yes = "wildtype", no = "no call")
seu$cc.ct <-  ifelse(apply(dplyr::select(seu@meta.data, all_of(cc.ct)), 1,
                           function(x) sum(grepl("mutant", x) > 0)), yes = "mutant", no = seu$cc.ct)
seu$tc.tt <- ifelse(apply(dplyr::select(seu@meta.data, all_of(tc.tt)), 1,
                          function(x) sum(grepl("wildtype", x) > 0)), yes = "wildtype", no = "no call")
seu$tc.tt <-  ifelse(apply(dplyr::select(seu@meta.data, all_of(tc.tt)), 1,
                           function(x) sum(grepl("mutant", x) > 0)), yes = "mutant", no = seu$tc.tt)
seu$cc.tt <- ifelse(apply(dplyr::select(seu@meta.data, all_of(cc.tt)), 1,
                          function(x) sum(grepl("wildtype", x) > 0)), yes = "wildtype", no = "no call")
seu$cc.tt <-  ifelse(apply(dplyr::select(seu@meta.data, all_of(cc.tt)), 1,
                           function(x) sum(grepl("mutant", x) > 0)), yes = "mutant", no = seu$cc.tt)

# Extract metadata
metadata_tib <- as_tibble(seu@meta.data, rownames = "cell")

# Wrangle. Note the inelegant use of the mut_type variable, which should be corrected.
mutation_types <- c("f_call", "p_call", "cc.ct", "tc.tt", "cc.tt")
plot_tib <- metadata_tib %>% filter(orig.ident != "BM") %>%
  select(orig.ident, cell, RF_pDC_score, bpdcn_sign_score, all_of(mutation_types)) %>%
  pivot_longer(cols = all_of(mutation_types), names_to = "mut_type", values_to = "mut_call") %>%
  mutate(mut_type = factor(mut_type, levels = c("f_call", "p_call", "cc.ct", "tc.tt", "cc.tt")))
n_call <- sum(grepl("mutant|wildtype", plot_tib$mut_call))

# Visualize
pdf("12.2.2_Scores-vs-Mut_type.pdf", width = 16, height = 4)

plot_tib %>%
  mutate(mut_call = factor(mut_call, levels = c("wildtype", "mutant", "no call"))) %>%
  arrange(mut_call) %>%
  mutate(my_order = c(sample(n_call), (n_call+1):nrow(plot_tib))) %>%
  arrange(desc(my_order)) %>%
  ggplot(aes(x = RF_pDC_score, y = bpdcn_sign_score, color = mut_call)) +
  #geom_point_rast(size = 0.5) +
  geom_point(size = 0.5) +
  scale_color_manual(values = mut_colors) +
  theme_bw() +
  theme(aspect.ratio = 1, panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  facet_wrap(~ mut_type, nrow = 1)

dev.off()


# ETV6 only ---------------------------------------------------------------------------------------

pdf("12.2.3_Pt14Dx_ETV6.R369W.pdf", width = 5, height = 5)

metadata_tib %>% filter(orig.ident == "Pt14Dx") %>%
  select(RF_pDC_score, bpdcn_sign_score, ETV6.R369W) %>%
  mutate(ETV6.R369W = factor(ETV6.R369W, levels = c("no call", "wildtype", "mutant"))) %>%
  arrange(ETV6.R369W) %>%
  ggplot(aes(x = RF_pDC_score, y = bpdcn_sign_score, color = ETV6.R369W)) +
  geom_point(size = 1.5) +
  scale_color_manual(values = mut_colors) +
  theme_bw() +
  theme(aspect.ratio = 1, panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"))

dev.off()


