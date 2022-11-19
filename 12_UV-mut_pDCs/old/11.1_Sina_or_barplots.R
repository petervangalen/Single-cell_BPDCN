# Peter van Galen, 220914
# Visualize UV-associated mutations in relation to cell types or BPDCN scores

library(tidyverse)
library(Seurat)
library(readxl)
library(ggforce)
library(cowplot)

rm(list=ls())

# Set working directory
setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/11_UV-mut_pDCs")

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

# Add malignant BPDCN cell module score
bpdcn_sign <- read.table("../09_pDC_expr/bpdcn_sign.txt")[,1]
seu <- AddModuleScore(seu, features = list(bpdcn_sign), name = "bpdcn_score")
colnames(seu@meta.data) <- gsub("bpdcn_score1", "bpdcn_score", colnames(seu@meta.data))

# Load genotyping information
genotyping_tables.tib <- read_excel("../04_XV-seq/XV-seq_overview.xlsx")
# Replace different MTAP entries with one, just as in 4.1_Add_GoT-XV_to_Seurat.R
genotyping_tables.tib$Mutation <- gsub("MTAP.rearr.*", "MTAP.rearr", genotyping_tables.tib$Mutation)


# Classify malignant/normal pDCs ------------------------------------------------------------------

# Sina plots of bpdcn_score across cell types colored by bone marrow involvement
seu@meta.data %>% ggplot(aes(x = CellType, y = bpdcn_score, color = bm_involvement)) +
  geom_sina(aes(group = CellType)) +
  scale_color_manual(values = group_colors) +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black"), axis.ticks = element_line(color = "black"))

# Based on this and looking at the data for weeks, it is fair to say malignant BPDCN cells should be classified as:
#    (1) All cells classified as pDCs in patients with bone marrow involvement
#    (2) All cells in any BPDCN patient  (with or without bone marrow involvement) with a BPDCN signature score exceeding 0.5
seu$Malignant_cell <- ifelse(seu$bm_involvement == "Yes" & seu$CellType == "pDC", yes = "Malignant_cell", no = "Normal_cell")
seu$Malignant_cell <- ifelse(seu$bm_involvement != "HD" & seu$bpdcn_score >= 0.5, yes = "Malignant_cell", no = seu$Malignant_cell)
seu$Malignant_cell <- factor(seu$Malignant_cell, levels = c("Normal_cell", "Malignant_cell"))

# Wrangle metadata --------------------------------------------------------------------------------

metadata_tib <- as_tibble(seu@meta.data, rownames = "cell")

# Select mutations to show. This is based on a variety of considerations. Outlined in the email thread "Bar plot" around 220918
# The first line is TC>TT (and CC>TT for ETV6) and the second line is CC>CT
current_mut <- c("HNRNPUL1.F559F", "MAP4K5.P667S", "MALAT1.chr11:65270399:G/A", "ETV6.R369W", "U2AF1.S34F",
                 "SMARCC1.chr3:47627735:G/A", "FAM98C.chr19:38899400:G/A", "SETX.chr9:135136947:G/A", "ACAP2.L97M", "IDH2.R140Q")

# Add genotyping summary
metadata_tib$UV_mutations <- ifelse(apply(dplyr::select(metadata_tib, all_of(current_mut)), 1,
                                     function(x) sum(grepl("wildtype", x) > 0)), yes = "wildtype", no = "no call")
metadata_tib$UV_mutations <- ifelse(apply(dplyr::select(metadata_tib, all_of(current_mut)), 1,
                                     function(x) sum(grepl("mutant", x) > 0)), yes = "mutant", no = metadata_tib$UV_mutations)

# Factorize & subset relevant information
metadata_tib$UV_mutations <- factor(metadata_tib$UV_mutations, levels = c("no call", "wildtype", "mutant"))
metadata_tib <- metadata_tib %>% dplyr::select(orig.ident, CellType, bpdcn_score, Malignant_cell, all_of(current_mut), UV_mutations)

# Visualize BPDCN scores vs. mutations ------------------------------------------------------------

pdf("11.1.1_SinaPlot.pdf", width = 3.5, height = 3.5)
metadata_tib %>%
  mutate(orig.order = sample(nrow(metadata_tib))) %>% arrange(orig.order) %>%
  filter(UV_mutations != "no call") %>%
  ggplot(aes(x = paste0("All genotyped cells (n = ", sum(metadata_tib$UV_mutations != "no call"), ")"),
             y = bpdcn_score, color = UV_mutations)) +
  geom_sina(aes(group = 1), size = 0.1) +
  scale_color_manual(values = c("#b0c4de", "#ffa500")) +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  theme_bw() +
  theme(aspect.ratio = 2, panel.grid = element_blank(),
        axis.title.x = element_blank(), axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"))
dev.off()


# Visualize UV-associated mutations in pDCs and non-pDCs ------------------------------------------
plot_tib <- metadata_tib %>%
  mutate(n_na = apply(metadata_tib, 1, function(x) sum(is.na(x)))) %>% filter(n_na < length(current_mut)) %>%
  group_by(Malignant_cell, UV_mutations) %>% count %>% ungroup %>% group_by(Malignant_cell) %>%
  mutate(`% of cells` = n/sum(n)*100) %>% 
  mutate(Malignant_cell = paste0(Malignant_cell, "\n(n = ", formatC(sum(n), big.mark=","), ")")) %>%
  mutate(Malignant_cell = factor(Malignant_cell, levels = unique(.$Malignant_cell)))
  
pdf("11.1.2_UV_in_Malignant_vs_Normal.pdf", width = 5, height = 6)
plot_tib %>% filter(UV_mutations != "no call") %>%
  ggplot(aes(x = Malignant_cell, y = `% of cells`, fill = UV_mutations)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#b0c4de", "#ffa500")) +
  theme_bw() +
  theme(aspect.ratio = 3, panel.grid = element_blank(), axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"), axis.title.x = element_blank())
dev.off()


# Same but split by mutation ----------------------------------------------------------------------

# Also start a list for the next section
split_ls <- vector(mode = "list")

# Visualize plots
pdf("11.1.3_Split_by_mutation.pdf", width = 10, height = 6)

for (m in current_mut) {
#m <- current_mut[1]
#m <- "FAM98C.chr19:38899400:G/A"

  # Wrangle and plot
  metadata_subset_tib <- metadata_tib %>% filter(get(m) %in% c("wildtype", "mutant", "no call")) %>%
    mutate(call = factor(.[,m,drop=T], levels = c("no call", "wildtype", "mutant"))) %>%
    arrange(call)
  
  for (p in unique(metadata_subset_tib$orig.ident)) {
    #p <- "Pt10Dx"
    #p <- "Pt1Dx"
    
    # Wrangle and plot
    metadata_subset_tib2 <- metadata_subset_tib %>% filter(orig.ident == p)
    
    p1 <- metadata_subset_tib2 %>% filter(get(m) %in% c("wildtype", "mutant")) %>%
      ggplot(aes(x = paste0("All genotyped cells (n = ", nrow(metadata_subset_tib2), ")"), y = bpdcn_score, color = call)) +
      geom_sina(aes(group = 1), size = 1) +
      scale_color_manual(values = c("#b0c4de", "#ffa500")) +
      guides(colour = guide_legend(override.aes = list(size=3))) +
      ggtitle(label = paste0(m, "\nin ", p)) +
      theme_bw() +
      theme(aspect.ratio = 2, panel.grid = element_blank(),
            axis.title.x = element_blank(), axis.ticks = element_line(color = "black"),
            axis.text = element_text(color = "black"),
            plot.title = element_text(hjust = 0.5))
    
    # Wrangle
    summary_tib <- metadata_subset_tib2 %>%
      group_by(Malignant_cell, call, .drop = F) %>% count() %>% ungroup %>% group_by(Malignant_cell) %>%
      mutate(`% of cells` = n/sum(n)*100) %>%
      mutate(Malignant_cell = paste0(Malignant_cell, "\n(n = ", formatC(sum(n), big.mark=","), ")")) %>%
      mutate(Malignant_cell = factor(Malignant_cell, levels = unique(.$Malignant_cell))) %>%
      mutate(plot_title = paste0(m, "\nin ", p))

    # Plot
    p2 <- summary_tib %>% filter(call != "no call") %>%
      ggplot(aes(x = Malignant_cell, y = `% of cells`, fill = call)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = c("#b0c4de", "#ffa500")) +
      ggtitle(label = paste0(m, "\nin ", p)) +
      theme_bw() +
      theme(aspect.ratio = 3, panel.grid = element_blank(), axis.ticks = element_line(color = "black"),
            axis.text = element_text(color = "black"), axis.title.x = element_blank(),
            plot.title = element_text(hjust = 0.5))

    # Visualize
    print(plot_grid(p1, p2))

    # Save
    split_ls <- c(split_ls, list(summary_tib))
    names(split_ls) <- c(names(split_ls)[-length(split_ls)], paste0(m, "\nin ", p))
  }
}
dev.off()


# Barplot split by mutation but shown together ----------------------------------------------------

pdf("11.1.4_Split_combined.pdf", width = 8, height = 16)
do.call(rbind, split_ls) %>% filter(call != "no call") %>%
  mutate(plot_title = factor(plot_title, levels = unique(plot_title))) %>%
  ggplot(aes(x = `% of cells`, y = Malignant_cell, fill = call)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#b0c4de", "#ffa500")) +
  facet_wrap(~ plot_title, scales = "free", ncol = 1, strip.position = "right") +
  theme_bw() +
  theme(aspect.ratio = 0.3, panel.grid = element_blank(),
        axis.ticks = element_line(color = "black", size = 0.3),
        axis.text = element_text(color = "black"),
        axis.title.y = element_blank(),
        panel.border = element_blank(), axis.line = element_line(size = 0.3),
        strip.background = element_blank(),
        strip.text.y = element_text(angle = 0, hjust = 0, face = "bold"))
dev.off()




