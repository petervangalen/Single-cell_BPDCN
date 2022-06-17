# Peter van Galen, 220520
# Make heatmaps of TET2 mutations in non-involved marrows

library(tidyverse)
library(Seurat)
library(readxl)
library(data.table)

rm(list=ls())

# Set working directory
setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/Github/9_Mutation_ratios")

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")

# Get results matrix for each non-involved sample
Samples <- c("Pt1Mrd", "Pt9Dx", "Pt10Dx", "Pt12Dx")

results.ls <- vector(mode = "list", length = length(Samples))
names(results.ls) <- Samples

for (s in Samples) {
  # Load data from previous script (9.1_All_Mutation_ratios.R)
  counts.df <- read.table(paste0("transcript_counts/", s, "_counts.txt"))
  wt_counts.df <- read.table(paste0("transcript_counts/", s, "_wt_counts.txt"))
  mut_counts.df <- read.table(paste0("transcript_counts/", s, "_mut_counts.txt"))
  
  # Reshape into a rectangular data format
  summary.dt <- data.table( reshape2::melt(t(counts.df)) )
  colnames(summary.dt) <- c("lineage", "mutation", "count")
  summary.dt$wt_counts <- t(wt_counts.df)
  summary.dt$mut_counts <- t(mut_counts.df)
  summary.dt$fraction <- summary.dt$mut_counts / summary.dt$count
  summary.dt$fraction[summary.dt$count <= 3] <- NA
  summary.dt$Sample <- s
  
  results.ls[[s]] <- summary.dt
}

# Prepare data table for plotting
tet2_results.ls <- lapply(results.ls, function(x) x[grepl("TET2", x$mutation),])
tet2_results.dt <- do.call(rbind, tet2_results.ls)
tet2_results.dt$mutation <- paste0(tet2_results.dt$Sample, "-", tet2_results.dt$mutation)

# Plot
pdf(file = paste0("TET2_uninvolved.pdf"), width = 6, height = 6)
print(
ggplot(data = tet2_results.dt, aes(x = lineage, y = mutation, fill = fraction, label = paste0(mut_counts, "/", count))) +
  geom_tile() +
  geom_text(size = 3, color = "#4D4D4D") + 
  scale_y_discrete(limits = rev(unique(tet2_results.dt$mutation))) +
  theme_classic() +
  theme(axis.line = element_blank(), aspect.ratio = 8/6,
        axis.text = element_text(color="black"),
        plot.title = element_text(hjust = 0.5),
        axis.ticks=element_blank(), panel.background = element_rect(colour = "black", size=1)) +
  ggtitle("TET2 mutations in uninvolved bone marrows") +
  scale_fill_gradient(low = "white", high = "#4B0092", na.value = "#DCDCDC", limits = c(0, 1))
)
dev.off()


