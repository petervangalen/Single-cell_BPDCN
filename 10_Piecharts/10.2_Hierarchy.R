# Peter van Galen, 220718
# Assess proportions of different cell types that harbor founder and progression mutations

library(tidyverse)
library(Seurat)
library(readxl)
library(data.table)
library(ggrepel)
library(ggforce)
library(cowplot)
library(ggrepel)
library(scatterpie)

rm(list=ls())

# Set working directory
setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/10_Piecharts")

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")
cell_colors <- popcol.tib$hex[1:21]
names(cell_colors) <- popcol.tib$pop[1:21]
donor_colors <- popcol.tib$hex[23:40]
names(donor_colors) <- popcol.tib$pop[23:40]
group_colors <- popcol.tib$hex[41:43]
names(group_colors) <- popcol.tib$pop[41:43]
mut_colors <- c(popcol.tib$hex[44:46], "#daa520")
names(mut_colors) <- c(popcol.tib$pop[44:46], "TC>TT")

# Load data
metadata_tib <- read_tsv("subset_metadata.txt")
metadata_tib$CellType <- factor(metadata_tib$CellType, levels = names(cell_colors))
metadata_tib$orig.ident <- factor(metadata_tib$orig.ident, levels = names(donor_colors)[c(1,8:18)])

# Load genotyping information
genotyping_tables.tib <- read_excel("../04_XV-seq/XV-seq_overview.xlsx")
# Replace different MTAP entries with one, just as in 4.1_Add_GoT-XV_to_Seurat.R
genotyping_tables.tib$Mutation <- gsub("MTAP.rearr.*", "MTAP.rearr", genotyping_tables.tib$Mutation)
f_muts <- filter(genotyping_tables.tib, `Founder or progression mutation` == "Founder") %>% .$Mutation
p_muts <- filter(genotyping_tables.tib, `Founder or progression mutation` == "Progression") %>% .$Mutation

# Coordinates for scatter pie charts
celltype_coordinates <- tribble(~CellType, ~x, ~y,
                                "HSPC", 4, 5,
                                "EarlyEry", 3, 4.5,
                                "LateEry", 2, 4,
                                "GMP", 3.5, 4,
                                "ProMono", 2.75, 3,
                                "Mono", 2, 2,
                                "ncMono", 3, 2,
                                "cDC", 4, 2,
                                "pDC", 5, 2,
                                "ProB", 5, 4.25,
                                "PreB", 5+1/3, 3.5,
                                "B", 5+2/3, 2.75,
                                "Plasma", 6, 2,
                                "CD4Naive", 7, 4,
                                "CD4Memory", 8, 4,
                                "CD8Naive", 7, 3,
                                "CD8Memory", 8, 3,
                                "CD8TermExh", 9, 3,
                                "GammaDeltaLike", 9, 4,
                                "NKT", 7, 2,
                                "NK", 8, 2)

# Hierarchy of pie charts -------------------------------------------------------------------------

# Summarize to get fraction mutated cells ... Two Options
  # Option 1: Founder
  pdf_name <- "Founder"
  mut_types <- c("wildtype", "mutant")
  summary_tib <- metadata_tib %>% filter(f_call != "no call") %>% group_by(CellType) %>%
    summarize(n = n(), wildtype = sum(f_call == "wildtype"), mutant = sum(f_call == "mutant"))

  # Option 2: Progression
  pdf_name <- "Progression"
  mut_types <- c("wildtype", "mutant", "TC>TT")
  summary_tib <- metadata_tib %>% filter(p_call != "no call") %>% group_by(CellType) %>%
    summarize(n = n(), wildtype = sum(p_call == "wildtype"), mutant = sum(p_call == "mutant"), `TC>TT` = sum(p_call == "TC>TT"))
  # For the text: "97.6% of the cells with progression mutations were pDCs"
  metadata_tib %>% filter(p_call %in% c("mutant", "TC>TT")) %>% 
    mutate(is_pdc = ifelse(.$CellType == "pDC", yes = TRUE, no = FALSE)) %>% .$is_pdc %>% mean
    mean(.$is_pdc)
  
# Visualize
pdf(paste0("10.2.1_", pdf_name, "_Mutations_in_CellTypes.pdf"))

# Merge data frames, plot scatterpie
plot_tib <- summary_tib %>% left_join(celltype_coordinates) %>%
  mutate(radius = log(n)/log(max(n))*0.5) %>% filter(n > 5)
ggplot() +
  geom_scatterpie(aes(x = x, y = y, r = radius), data = plot_tib, cols = mut_types, color = NA) +
  geom_scatterpie_legend(plot_tib$radius, x = 9, y = 1.5, n = 3) +
  scale_fill_manual(values = mut_colors[mut_types]) +
  geom_label(aes(x = x, y = y+radius+0.15, label = CellType), data = plot_tib, size = 3, label.padding = unit(0.1, "lines")) +
  coord_cartesian(xlim = c(1,10), ylim = c(1,6)) +
  ggtitle(label = paste0(pdf_name, " mutations")) +
  theme_bw() +
  theme(aspect.ratio = 5/9, panel.grid = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank())

# Stacked bar plot
summary_tib %>% pivot_longer(cols = all_of(mut_types), names_to = "call", values_to = "count") %>%
  ggplot(aes(x = CellType, y = count, fill = call)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = mut_colors[mut_types]) +
  ggtitle(paste("Number of cells with", tolower(pdf_name), "mutations")) +
  theme_bw() +
  theme(aspect.ratio = 2/3, axis.title.x = element_blank(), axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))

dev.off()

# For the legend (for the progression mutations, resize the smallest circle to match GMP size)
min(plot_tib$n); max(plot_tib$n)


# Split pie scatterplots by patients --------------------------------------------------------------

# This really only makes sense for patients with FOUNDER and PROGRESSION mutations
patient_list <- list(Pt1 = c("Pt1Dx", "Pt1Rem"), Pt10 = c("Pt10Dx", "Pt10Rel"), Pt12 = c("Pt12Dx", "Pt12Rel"))

for (p in names(patient_list)) {
  #p <- "Pt10"
  metadata_subset_tib <- metadata_tib %>% filter(orig.ident %in% patient_list[[p]])
  
  for (pdf_name in c("Founder", "Progression")) {
    #pdf_name <- "Progression"

if (pdf_name == "Founder") {
  mut_types <- c("wildtype", "mutant")
  summary_tib <- metadata_subset_tib %>% filter(f_call != "no call") %>% group_by(CellType) %>%
    summarize(n = n(), wildtype = sum(f_call == "wildtype"), mutant = sum(f_call == "mutant"))
} else if (pdf_name == "Progression") {
  mut_types <- c("wildtype", "mutant", "TC>TT")
  summary_tib <- metadata_subset_tib %>% filter(p_call != "no call") %>% group_by(CellType) %>%
    summarize(n = n(), wildtype = sum(p_call == "wildtype"), mutant = sum(p_call == "mutant"), `TC>TT` = sum(p_call == "TC>TT"))
}
    
# Visualize
pdf(paste0("split_by_patient/10.2.1_", p, "_", pdf_name, ".pdf"))

# Merge data frames, plot scatterpie
plot_tib <- summary_tib %>% left_join(celltype_coordinates) %>%
  mutate(radius = log(n)/log(max(n))*0.5)# %>% filter(n > 5)

print( 
  ggplot() +
    geom_scatterpie(aes(x = x, y = y, r = radius), data = plot_tib, cols = mut_types, color = NA) +
    geom_scatterpie_legend(plot_tib$radius, x = 9, y = 1.5, n = 3) +
    scale_fill_manual(values = mut_colors[mut_types]) +
    geom_label(aes(x = x, y = y+radius+0.15, label = CellType), data = plot_tib, size = 3, label.padding = unit(0.1, "lines")) +
    coord_cartesian(xlim = c(1,10), ylim = c(1,6)) +
    ggtitle(label = paste0(pdf_name, " mutations")) +
    theme_bw() +
    theme(aspect.ratio = 5/9, panel.grid = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), axis.title = element_blank())
)
dev.off()

  }
}
