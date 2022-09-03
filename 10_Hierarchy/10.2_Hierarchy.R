# Peter van Galen, 220718
# Assess proportions of different cell types that harbor founder and progression mutations

library(tidyverse)
library(Seurat)
library(readxl)
library(data.table)
library(ggrepel)
library(ggforce)
library(cowplot)
#library(viridis)

rm(list=ls())

# Set working directory
setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/10_Hierarchy")

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")
cell_colors <- popcol.tib$hex[1:21]
names(cell_colors) <- popcol.tib$pop[1:21]
donor_colors <- popcol.tib$hex[23:40]
names(donor_colors) <- popcol.tib$pop[23:40]
group_colors <- popcol.tib$hex[41:43]
names(group_colors) <- popcol.tib$pop[41:43]
mut_colors <- popcol.tib$hex[44:46]
names(mut_colors) <- popcol.tib$pop[44:46]

# Load Seurat objects
seurat_files <- list.files("../04_XV-seq", pattern = "*.rds", full.names = T)
seu_ls <- lapply(seurat_files, function(x) readRDS(x))
names(seu_ls) <- cutf(basename(seurat_files), d = "_")
seu <- merge(seu_ls[[1]], seu_ls[2:length(seu_ls)])
seu$orig.ident2 <- ifelse(grepl("BM", seu$orig.ident), yes = cutf(seu$replicate, d = "\\."), no = seu$orig.ident)
seu$CellType <- factor(seu$CellType, levels = levels(seu_ls[[1]]$CellType))

# Define sample groups
#metadata_tib <- as_tibble(seu@meta.data, rownames = "cell")
#skin_only_samples <- unique(filter(metadata_tib, bm_involvement == "No")$orig.ident) %>% .[c(3,4,5,1,2)]
#bm_involvement_samples <- unique(filter(metadata_tib, bm_involvement == "Yes")$orig.ident) %>% .[c(6,1,2,3,4,5)]

# Load genotyping information
genotyping_tables.tib <- read_excel("../04_XV-seq/XV-seq_overview.xlsx")
# Replace different MTAP entries with one, just as in 4.1_Add_GoT-XV_to_Seurat.R
genotyping_tables.tib$Mutation <- gsub("MTAP.rearr.*", "MTAP.rearr", genotyping_tables.tib$Mutation)

# Add malignant BPDCN cell module score, extract metadata
bpdcn_sign <- read.table("../09_pDC_expr/bpdcn_sign.txt")[,1]
seu <- AddModuleScore(seu, features = list(bpdcn_sign), name = "bpdcn_score")
metadata_tib <- as_tibble(seu@meta.data, rownames = "cell") %>%
  dplyr::select(cell, orig.ident, project.umap.x, project.umap.y, genotyping_tables.tib$Mutation,
                CellType, bm_involvement, bpdcn_score1) %>%
  rename(bpdcn_score = "bpdcn_score1")


# Add genotyping information ----------------------------------------------------------------------

# Add which cells have founder mutations
f_muts <- filter(genotyping_tables.tib, `Founder or progression mutation` == "Founder") %>% .$Mutation
metadata_tib$f_call <- ifelse(apply(dplyr::select(metadata_tib, all_of(f_muts)), 1, 
                                    function(x) sum(grepl("wildtype", x) > 0)), yes = "wildtype", no = "no call")
metadata_tib$f_call <- ifelse(apply(dplyr::select(metadata_tib, all_of(f_muts)), 1,
                                    function(x) sum(grepl("mutant", x) > 0)), yes = "mutant", no = metadata_tib$f_call)
metadata_tib$f_call <- factor(metadata_tib$f_call, levels = c("no call", "wildtype", "mutant"))
# Add which cells have progession mutations
p_muts <- filter(genotyping_tables.tib, `Founder or progression mutation` == "Progression") %>% .$Mutation
metadata_tib$p_call <- ifelse(apply(dplyr::select(metadata_tib, all_of(p_muts)), 1, 
                                    function(x) sum(grepl("wildtype", x) > 0)), yes = "wildtype", no = "no call")
metadata_tib$p_call <- ifelse(apply(dplyr::select(metadata_tib, all_of(p_muts)), 1,
                                    function(x) sum(grepl("mutant", x) > 0)), yes = "mutant", no = metadata_tib$p_call)
metadata_tib$p_call <- factor(metadata_tib$p_call, levels = c("no call", "wildtype", "mutant"))


# Add a new cell type (BPDCN) and reclassify cells accordingly ------------------------------------

# Add new cell type: BPDCN. The validity of these procedures is described in 10.1_Reclassify.R
metadata_tib <- metadata_tib %>% mutate(CellTypeRefined = case_when(
  CellType == "pDC" & bm_involvement == "Yes" ~ "BPDCN",
  bpdcn_score > 0.5 ~ "BPDCN",
  TRUE ~ CellType)) %>%
  mutate(CellTypeRefined = factor(CellTypeRefined, levels = c(levels(metadata_tib$CellType), "BPDCN")))
metadata_tib$CellType %>% table
metadata_tib$CellTypeRefined %>% table


# Hierarchy of pie charts -------------------------------------------------------------------------


# Summarize to get fraction mutated cells ... Two Options
  # Option 1: Founder
  pdf_name <- "Founder"
  summary_tib <- metadata_tib %>% filter(f_call != "no call") %>%
    group_by(CellTypeRefined) %>% summarize(n = n(), wildtype = sum(f_call == "wildtype"), mutant = sum(f_call == "mutant")) %>%
    mutate(mutated_cells = mutant/(wildtype+mutant))

  # Option 2: Progression
  pdf_name <- "Progression"
  summary_tib <- metadata_tib %>% filter(p_call != "no call") %>%
    group_by(CellTypeRefined) %>% summarize(n = n(), wildtype = sum(p_call == "wildtype"), mutant = sum(p_call == "mutant")) %>%
    mutate(mutated_cells = mutant/(wildtype+mutant))

# Stacked bar plot
pdf(paste0("10.2.1_Barplot_", pdf_name, ".pdf"))
summary_tib %>% pivot_longer(cols = c(wildtype, mutant), names_to = "call", values_to = "count") %>%
  ggplot(aes(x = CellTypeRefined, y = count, fill = call)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = mut_colors) +
  ggtitle(paste("Number of cells with", tolower(pdf_name), "mutations")) +
  theme_bw() +
  theme(aspect.ratio = 2/3, axis.title.x = element_blank(), axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))
dev.off()

# Use facet wrap to plot pie charts
pdf(paste0("10.2.2_PieArray_", pdf_name, ".pdf"))
summary_tib %>% mutate(wildtype_cells = 1-mutated_cells) %>%
  pivot_longer(cols = c(mutated_cells, wildtype_cells), names_to = "call", values_to = "proportion") %>%
  ggplot(aes(x = "", y = proportion, fill = call)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(wildtype_cells = "#32cd32", mutated_cells = "#dc143c")) +
  coord_polar("y") +
  ggtitle(paste("Proportion of cells with", tolower(pdf_name), "mutations")) +
  facet_wrap(~ CellTypeRefined) +
  theme_bw() +
  theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
        panel.grid = element_blank())
dev.off()

# To customize size, make a list of ggplots
pies_ls <- lapply(summary_tib$CellTypeRefined, function(x) {
  summary_tib %>% mutate(wildtype_cells = 1-mutated_cells) %>%
    filter(CellTypeRefined == x) %>%
    pivot_longer(cols = c(mutated_cells, wildtype_cells), names_to = "call", values_to = "proportion") %>%
    ggplot(aes(x = "", y = proportion, fill = call)) + 
    geom_bar(stat = "identity", show.legend = F) +
    scale_fill_manual(values = c(wildtype_cells = "#32cd32", mutated_cells = "#dc143c")) +
    coord_polar("y") +
    ggtitle(x) +
    theme_bw() +
    theme(axis.title = element_blank(), axis.ticks = element_blank(),
          axis.text = element_blank(), panel.border = element_blank(),
          plot.title = element_text(hjust = 0.5),
          panel.grid = element_blank())
    })
names(pies_ls) <- summary_tib$CellTypeRefined

# Save PDF
pdf(file = paste0("10.2.3_PieSizes_", pdf_name, ".pdf"), width = 20, height = 3)
plot_grid(plotlist = pies_ls, nrow = 1, rel_widths = log(summary_tib$n))
dev.off()

  

