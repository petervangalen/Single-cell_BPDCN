# Peter van Galen, 221204
# Visualize UV-associated mutations in relation to cell types or BPDCN scores

library(tidyverse)
library(Seurat)
library(readxl)
library(ggforce)
library(cowplot)

rm(list=ls())

# Set working directory
setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/12_Mut_vs_scores")

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")
cell_colors <- popcol.tib$hex[1:21]
names(cell_colors) <- popcol.tib$pop[1:21]
mut_colors <- popcol.tib$hex[44:46]
names(mut_colors) <- popcol.tib$pop[44:46]
group_colors <- popcol.tib$hex[41:43]
names(group_colors) <- popcol.tib$pop[41:43]

# Load Seurat objects
seurat_files <- list.files("../04_XV-seq", pattern = "*.rds", full.names = T)
seu_ls <- lapply(seurat_files, function(x) readRDS(x))
names(seu_ls) <- cutf(basename(seurat_files), d = "_")
seu <- merge(seu_ls[[1]], seu_ls[2:length(seu_ls)])
seu$CellType <- factor(seu$CellType, levels = levels(seu_ls[[1]]$CellType))

# Add metadata from 11_pDC_expr
MalignantCalls_df <- read.table("../11_pDC_expr/11.3_MalignantCalls_Final.txt", header = T, row.names = "cell")
all(rownames(seu@meta.data) == rownames(MalignantCalls_df))
seu <- AddMetaData(seu, MalignantCalls_df[, c("bpdcn_sign_score", "RF_pDC_score", "is_malignant")])
# For this analysis, healthy pDCs and healthy other cells can be considered together:
seu$is_malignant[seu$is_malignant == "Other"] <- "Healthy"
seu$is_malignant <- factor(seu$is_malignant, levels = c("Healthy", "Premalignant", "Malignant"))

# Load genotyping information
genotyping_tables.tib <- read_excel("../04_XV-seq/XV-seq_overview.xlsx")
# Replace different MTAP entries with one, just as in 4.1_Add_GoT-XV_to_Seurat.R
genotyping_tables.tib$Mutation <- gsub("MTAP.rearr.*", "MTAP.rearr", genotyping_tables.tib$Mutation)


# Wrangle metadata --------------------------------------------------------------------------------

metadata_tib <- as_tibble(seu@meta.data, rownames = "cell")

# Text: "Out of 4,263 bone marrow cells with UV progression mutations, 4,100 (96.2%) were classified as malignant BPDCN cells"
all_uv_muts <- genotyping_tables.tib %>%
  filter(`CC>CT (UV-associated)`	== "Yes" | `TC>TT (UV-associated)` == "Yes" | `CC>TT (UV-specific)` == "Yes") %>%
  .$Mutation %>% unique
metadata_tib$any_uv_mut <- apply(dplyr::select(metadata_tib, all_of(all_uv_muts)), 1, function(x) sum(x %in% "mutant") > 0)
metadata_tib %>% filter(any_uv_mut == T) %>% nrow
metadata_tib %>% filter(any_uv_mut == T, is_malignant == "Malignant") %>% nrow

# Select UV-associated mutations to show. See also the email thread "Bar plot" around 220918
selected_uv_mut <- c("SMARCC1.chr3:47627735:G/A", # Patient 1 CC>CT
                     "FAM98C.chr19:38899400:G/A", # Patient 1 CC>CT
                     "SETX.chr9:135136947:G/A",   # Patient 1 CC>CT
                     "HNRNPUL1.F559F",            # Patient 10 TC>TT
                     "MAP4K5.P667S",              # Patient 10 TC>TT
                     "ACAP2.L97M",                # Patient 10 CC>CT
                     "IDH2.R140Q",                # Patient 14 CC>CT
                     "ETV6.R369W")                # Patient 14 CC>TT
# Omit - reason:
#Pt10 "TET2.H1380Y" - founder, not specific to BPDCN cells
#Pt10 "PRKDC.chr8:48713410:C/T" - subclonal expansion that disappears
#Pt10 "MALAT1.chr11:65270399:G/A" - very few enriched transcripts at diagnosis
#Pt15 "U2AF1.S34F" - looks good but UMAP suggests false positives by soup RNA or doublets, possibly due to very high expression
#Pt14 "EZH2.P132L" - only detected in non-pDCs with a low BPDCN score

# Add genotyping summary
metadata_tib$UV_mutations <- ifelse(apply(dplyr::select(metadata_tib, all_of(selected_uv_mut)), 1,
                                     function(x) sum(grepl("wildtype", x) > 0)), yes = "wildtype", no = "no call")
metadata_tib$UV_mutations <- ifelse(apply(dplyr::select(metadata_tib, all_of(selected_uv_mut)), 1,
                                     function(x) sum(grepl("mutant", x) > 0)), yes = "mutant", no = metadata_tib$UV_mutations)

# Factorize & subset relevant information
metadata_tib$UV_mutations <- factor(metadata_tib$UV_mutations, levels = c("no call", "wildtype", "mutant"))
metadata_tib <- metadata_tib %>% dplyr::select(cell, orig.ident, CellType, bpdcn_sign_score, is_malignant, RF_pDC_score, all_of(selected_uv_mut), UV_mutations)


# Visualize barplots of each mutation in healthy, premalignant, and malignant cells ---------------

# Also start a list for the next section
split_ls <- vector(mode = "list")

# Wrangle
for (m in selected_uv_mut) {
#m <- selected_uv_mut[1]

  # Subset metadata for cells that were genotyped for the current mutation
  metadata_subset_tib <- metadata_tib %>% filter(get(m) %in% c("wildtype", "mutant", "no call")) %>%
    mutate(call = factor(.[, m, drop=T], levels = c("no call", "wildtype", "mutant"))) %>%
    arrange(call)
  
  for (p in unique(metadata_subset_tib$orig.ident)) {
    #p <- "Pt1Dx"
    
    # For the current mutation/patient, count the total number of cells, then determine the proportion,
    # then add labels for plotting
    summary_tib <- metadata_subset_tib %>% filter(orig.ident == p) %>%
      group_by(is_malignant, call, .drop = F) %>% count() %>% ungroup %>% group_by(is_malignant) %>%
      mutate(`% of cells` = n/sum(n)*100) %>%
      mutate(axis_label = paste0(is_malignant, " (n = ", formatC(sum(n), big.mark=","), ")")) %>%
      mutate(axis_label = factor(axis_label, levels = unique(.$axis_label))) %>%
      mutate(plot_title = paste0(cutf(m, d = "\\."), "\nin ", p))

    # Save
    split_ls <- c(split_ls, list(summary_tib))
    names(split_ls) <- c(names(split_ls)[-length(split_ls)], paste0(m, "\nin ", p))
  }
}


# Visualize bar plots

pdf("12.4.1_Split_combined.pdf", width = 15, height = 6)

do.call(rbind, split_ls) %>% na.omit() %>% filter(`% of cells` != 0) %>%
  mutate(plot_title = factor(plot_title, levels = unique(plot_title))) %>%
  ggplot(aes(x = axis_label, y = `% of cells`, fill = call)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = mut_colors) +
  facet_wrap(~ plot_title, scales = "free_x", nrow = 1) +
  theme_bw() +
  annotate("segment", x = 0, xend = 0, y = -Inf, yend = Inf) +
  theme(aspect.ratio = 4, panel.grid = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.ticks = element_line(color = "black", size = 0.3),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        panel.border = element_blank(), axis.line = element_line(size = 0.3),
        strip.background = element_blank(),
        strip.text.x = element_text(angle = 0, hjust = 0.5, face = "bold"),
        panel.spacing = unit(1, "lines"))

dev.off()




