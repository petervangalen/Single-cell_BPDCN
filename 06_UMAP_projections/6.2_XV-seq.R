# Peter van Galen, 220707
# Compare genotyping efficiency of different approaches

library(tidyverse)
library(Seurat)
library(readxl)
library(data.table)
library(gridExtra)
library(KernSmooth)
library(MASS)

setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/06_UMAP_projections")

rm(list=ls())

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")
cell_colors <- popcol.tib$hex[1:22]
names(cell_colors) <- popcol.tib$pop[1:22]
mut_colors <- popcol.tib$hex[45:47]
names(mut_colors) <- popcol.tib$pop[45:47]

# Load Seurat objects
seurat_files <- list.files("../04_XV-seq", pattern = "*.rds", full.names = T)
seu.ls <- lapply(seurat_files, function(x) readRDS(x))
names(seu.ls) <- cutf(basename(seurat_files), d = "_")

# Merge BPDCN sample Seurat objects (omit BM)
seu_bpdcn_merge <- merge(seu.ls[[2]], seu.ls[3:length(seu.ls)])

# Load genotyping information
genotyping_tables.tib <- read_excel("../04_XV-seq/XV-seq_overview.xlsx")
# Replace different MTAP entries with one, just as in 4.1_Add_GoT-XV_to_Seurat.R
genotyping_tables.tib$Mutation <- gsub("MTAP.rearr.*", "MTAP.rearr", genotyping_tables.tib$Mutation)

# Calculate healthy donors contour density
n <- 100
bins <- bkde2D(seu.ls$BM@reductions$umap@cell.embeddings, bandwidth = c(0.5, 0.5), gridsize = c(n, n))
meltbins <- data.table( reshape2::melt(bins$fhat) )
meltbins$Var1 <- rep(bins$x1, n)
meltbins$Var2 <- rep(bins$x2, each = n)

# CHOOSE ONE
muts <- genotyping_tables.tib$Mutation
pdf_name <- "all_mutations"
# OR
muts <- filter(genotyping_tables.tib, `Founder or progression mutation` == "Progression") %>% .$Mutation
pdf_name <- "progression_mutations"


# Project UMAPs with wildtype and mutated cells  --------------------------------------------------

# Calculate which cells were mutated
metadata_tib <- as_tibble(seu_bpdcn_merge@meta.data) %>%
  dplyr::select(orig.ident, project.umap.x, project.umap.y, genotyping_tables.tib$Mutation, CellType, bm_involvement)
# Add metadata column indicating the mutational status of each cell
metadata_tib$call <- ifelse(apply(dplyr::select(metadata_tib, all_of(muts)), 1, 
                                  function(x) sum(grepl("wildtype", x) > 0)), yes = "wildtype", no = "no call")
metadata_tib$call <- ifelse(apply(dplyr::select(metadata_tib, all_of(muts)), 1,
                                  function(x) sum(grepl("mutant", x) > 0)), yes = "mutant", no = metadata_tib$call)

# Reorder for plotting
no_call <- metadata_tib %>% filter(call == "no call")
called <- metadata_tib %>% filter(call != "no call")
called <- called %>% mutate(my_order = sample(1:nrow(called))) %>% arrange(my_order) %>% mutate(my_order = NULL)
metadata_order_tib <- rbind(no_call, called)

pdf(paste0("6.2_XV-seq_", pdf_name, ".pdf"), width = 12, height = 5)

# Plot all cells
p1 <- ggplot() +
  geom_contour(data = meltbins, mapping = aes(x = Var1, y = Var2, z = value, color = after_stat(level)), bins = 10, size = 1) +
  scale_color_gradient2(low = "grey", high = "black", guide="none") +
  geom_point(data = metadata_order_tib, mapping = aes(x = project.umap.x, y = project.umap.y),
             color = mut_colors[metadata_order_tib$call], size = 0.8) +
  ggtitle("All") +
  theme(aspect.ratio = 1, axis.line=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(),
        axis.ticks=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(),
        panel.border=element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14))

# Plot skin-only patients
metadata_order_skin_only <- metadata_order_tib %>% filter(bm_involvement == "No")
p2 <- ggplot() +
  geom_contour(data = meltbins, mapping = aes(x = Var1, y = Var2, z = value, color = after_stat(level)), bins = 10, size = 1) +
  scale_color_gradient2(low = "grey", high = "black", guide="none") +
  geom_point(data = metadata_order_skin_only, mapping = aes(x = project.umap.x, y = project.umap.y),
             color = mut_colors[metadata_order_skin_only$call], size = 0.8) +
  ggtitle("Skin-only") +
  theme(aspect.ratio = 1, axis.line=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(),
        axis.ticks=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(),
        panel.border=element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14))

# Plot BM involvement patients
metadata_order_bm_involvement <- metadata_order_tib %>% filter(bm_involvement == "Yes")
p3 <- ggplot() +
  geom_contour(data = meltbins, mapping = aes(x = Var1, y = Var2, z = value, color = after_stat(level)), bins = 10, size = 1) +
  scale_color_gradient2(low = "grey", high = "black", guide="none") +
  geom_point(data = metadata_order_bm_involvement, mapping = aes(x = project.umap.x, y = project.umap.y),
             color = mut_colors[metadata_order_bm_involvement$call], size = 0.8) +
  ggtitle("BM involvement") +
  theme(aspect.ratio = 1, axis.line=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(),
        axis.ticks=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(),
        panel.border=element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14))

grid.arrange(p1, p2, p3, ncol = 3)

dev.off()

# For the text: how many mutated and wild-type cells in skin-only patients? (use for all mutations)
metadata_tib %>% filter(bm_involvement == "No") %>% .$call %>% table
# How many mutated and wild-type cells in bm involvement patients? (use for progression mutations)
metadata_tib %>% filter(bm_involvement == "Yes") %>% .$call %>% table


# Overall heatmap of mutation ratio ---------------------------------------------------------------

# See related scripts in 8.1_All_Mutation_ratios.R

# Add a column specifying more broad cell types
metadata_tib$lineage[metadata_tib$CellType %in% c("HSPC", "GMP")] <- "HSPC"
metadata_tib$lineage[metadata_tib$CellType %in% c("EarlyEry", "LateEry")] <- "Ery"
metadata_tib$lineage[metadata_tib$CellType %in% c("ProMono", "Mono", "ncMono", "cDC", "pDC")] <- "Myeloid"
metadata_tib$lineage[metadata_tib$CellType %in% c("ProB", "PreB", "B", "Plasma")] <- "B.cells"
metadata_tib$lineage[metadata_tib$CellType %in% c("CD4Naive", "CD4Memory", "CD8Naive", "CD8Memory", "CD8TermExh", "GammaDeltaLike")] <- "T.cells"
metadata_tib$lineage[metadata_tib$CellType %in% c("NKT", "NK")] <- "NK"
metadata_tib <- metadata_tib %>% mutate(lineage = factor(lineage, levels = c("HSPC", "Ery", "Myeloid", "NK", "B.cells", "T.cells")))

# Visualize
pdf(paste0("6.2_Mutation_ratio_heatmap_", pdf_name, ".pdf"), width = 10, height = 2)

# All: Make table, then visualize
count_tib <- metadata_tib %>% group_by(lineage, call) %>% count()
count2_tib <- count_tib %>% pivot_wider(id_cols = lineage, names_from = call, values_from = n) %>%
  mutate(total = sum(mutant, wildtype))
p1 <- count2_tib %>%
  ggplot(aes(x = lineage, y = 1, fill = mutant/total, label = paste0(mutant, "\n", total))) +
  geom_tile() +
  geom_text(size = 3) + 
  ggtitle("All") +
  theme_classic() +
  theme(axis.line = element_blank(), aspect.ratio = 1/6, axis.title = element_blank(), axis.text.y = element_blank(),
        axis.text.x = element_text(color = "black"), axis.ticks = element_blank(),
        panel.background = element_rect(colour = "black", size=1)) +
  scale_fill_gradient(low = "white", high = "#4B0092", limits = c(0, 1))

# Skin-only: Make table, then visualize
count_tib <- metadata_tib %>% filter(bm_involvement == "No") %>% group_by(lineage, call) %>% count()
count2_tib <- count_tib %>% pivot_wider(id_cols = lineage, names_from = call, values_from = n) %>%
  mutate(total = sum(mutant, wildtype))
p2 <- count2_tib %>%
  ggplot(aes(x = lineage, y = 1, fill = mutant/total, label = paste0(mutant, "\n", total))) +
  geom_tile() +
  geom_text(size = 3) +
  ggtitle("Skin-only") +
  theme_classic() +
  theme(axis.line = element_blank(), aspect.ratio = 1/6, axis.title = element_blank(), axis.text.y = element_blank(),
        axis.text.x = element_text(color = "black"), axis.ticks = element_blank(),
        panel.background = element_rect(colour = "black", size=1)) +
  scale_fill_gradient(low = "white", high = "#4B0092", limits = c(0, 1))

# BM involvement: Make table, then visualize
count_tib <- metadata_tib %>% filter(bm_involvement == "Yes")  %>% group_by(lineage, call) %>% count()
count2_tib <- count_tib %>% pivot_wider(id_cols = lineage, names_from = call, values_from = n) %>%
  mutate(total = sum(mutant, wildtype))
p3 <- count2_tib %>%
  ggplot(aes(x = lineage, y = 1, fill = mutant/total, label = paste0(mutant, "\n", total))) +
  geom_tile() +
  geom_text(size = 3) + 
  ggtitle("BM involvement") +
  theme_classic() +
  theme(axis.line = element_blank(), aspect.ratio = 1/6, axis.title = element_blank(), axis.text.y = element_blank(),
        axis.text.x = element_text(color = "black"), axis.ticks = element_blank(),
        panel.background = element_rect(colour = "black", size=1)) +
  scale_fill_gradient(low = "white", high = "#4B0092", limits = c(0, 1))

grid.arrange(p1, p2, p3, ncol = 3)

dev.off()

# Chi square test
count_mat <- count2_tib %>% mutate(lineage_broad = ifelse(lineage %in% c("HSPC", "Ery", "Myeloid"), yes = "myeloerythroid", no = "lymphoid")) %>%
  group_by(lineage_broad) %>% summarize(wt = sum(wildtype), mut = sum(mutant)) %>% data.frame(row.names = 1) %>% t()

chisq.test(count_mat)


# Heatmap of mutation ratio ---------------------------------------------------------------

stopifnot(pdf_name == "progression_mutations")

# Add a column specifying pDCs vs. all others
metadata_tib <- metadata_tib %>% filter(bm_involvement == "Yes") %>%
  mutate(pDC = ifelse(CellType == "pDC", yes = "pDC", no = "Other")) %>%
  mutate(pDC = factor(pDC, levels = c("pDC", "Other")))

count_tib <- metadata_tib %>% group_by(pDC, call) %>% count()
count2_tib <- count_tib %>% pivot_wider(id_cols = pDC, names_from = call, values_from = n) %>%
  mutate(total = sum(mutant, wildtype))

# Make visualize
pdf(paste0("6.2_Mutation_ratio_heatmap_progression_mutations_pDCs.pdf"), width = 3, height = 2)

# All: Make table, then visualize
count2_tib %>%
  ggplot(aes(x = pDC, y = 1, fill = mutant/total, label = paste0(mutant, "\n", total))) +
  geom_tile() +
  geom_text(size = 3) + 
  ggtitle("Progression mutations \nin pDCs vs. all others") +
  theme_classic() +
  theme(axis.line = element_blank(), aspect.ratio = 1/4, axis.title = element_blank(), axis.text.y = element_blank(),
        axis.text.x = element_text(color = "black"), axis.ticks = element_blank(),
        panel.background = element_rect(colour = "black", size=1)) +
  scale_fill_gradient(low = "white", high = "#4B0092", limits = c(0, 1))

dev.off()
