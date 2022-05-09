# Peter van Galen, 220508
# Make heatmaps of the ratio of wild-type to mutated transcripts in each population

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

# Load Seurat objects
seurat_files <- list.files("../7_XV-seq", pattern = "*.rds", full.names = T)
seu_bpdcn.ls <- lapply(seurat_files, function(x) readRDS(x))
names(seu_bpdcn.ls) <- gsub("_Seurat_Anno.rds", "", cutf(seurat_files, d = "/", f = 3))

# Generate data frame for heatmap annotations
genotyping_tables.tib <- read_excel("../7_XV-seq/FilteredCells_files.xlsx")
# Replace different MTAP primers with one, just as in 7_Add_GoT-XV_to_Seurat.R
genotyping_tables.tib$Mut <- gsub("MTAP.rearr.*", "MTAP.rearr", genotyping_tables.tib$Mut)
genotyping_tables.tib <- genotyping_tables.tib %>% select(Sample, Mut) %>% unique

# -------------------------------------------------------------------------------------------------

for (n in 1:nrow(genotyping_tables.tib)) {

sample <- genotyping_tables.tib$Sample[n]
seu <- seu_bpdcn.ls[[sample]]
muts <- filter(genotyping_tables.tib, Sample == sample)$Mut

# Extract the necessary columns from Seurat metadata
info.df <- seu@meta.data[,c("CellType", muts)]

# Add a column specifying more broad cell types: info.df$lineage
info.df$lineage[info.df$CellType %in% c("HSC", "Prog")] <- "HSPC"
info.df$lineage[info.df$CellType %in% c("EarlyE", "LateE")] <- "Ery"
info.df$lineage[info.df$CellType %in% c("ProMono", "Mono", "ncMono", "cDC", "pDC")] <- "Myeloid"
info.df$lineage[info.df$CellType %in% c("ProB", "PreB", "B", "Plasma")] <- "B.cells"
info.df$lineage[info.df$CellType %in% c("T", "CTL")] <- "T.cells"
info.df$lineage[info.df$CellType == "NK"] <- "NK"

# 3. Split data frame into a list of data frames.
split.ls <- split(info.df, f = info.df$lineage)
# Each data frame in the list has a subset of cells as rows, and the same columns as split.df
# Remove celltype and lineage columns which we don't need
split.ls <- lapply(split.ls, function(x) x[,-(grep("CellType|lineage", colnames(x)))])

# Make a data frame with the number of genotyped cells in each lineage
counts.df <- as.data.frame(do.call(cbind, lapply(split.ls, function(x) apply(x, 2, function(y) sum(y != "no call")))))
# Reorder columns as HSPC-Ery-Myeloid-NK-B-T
counts.df <- counts.df[,c("HSPC", "Ery", "Myeloid", "NK", "B.cells", "T.cells")]

# Make a data frame with the number of wild-type cells for each mutation in each lineage
wt_counts.df <- as.data.frame(do.call(cbind, lapply(split.ls, function(x) apply(x, 2, function(y) sum(y == "wildtype")))))
wt_counts.df <- wt_counts.df[,c("HSPC", "Ery", "Myeloid", "NK", "B.cells", "T.cells")]

# Make a data frame with the number of mutated cells for each mutation in each lineage
mut_counts.df <- as.data.frame(do.call(cbind, lapply(split.ls, function(x) apply(x, 2, function(y) sum(y == "mutant")))))
mut_counts.df <- mut_counts.df[,c("HSPC", "Ery", "Myeloid", "NK", "B.cells", "T.cells")]

# Check that it makes sense
mut_counts.df + wt_counts.df == counts.df

# Reshape into a rectangular data format
summary.dt <- data.table( reshape2::melt(t(counts.df)) )
colnames(summary.dt) <- c("lineage", "mutation", "count")
summary.dt$wt_counts <- t(wt_counts.df)
summary.dt$mut_counts <- t(mut_counts.df)
summary.dt$fraction <- summary.dt$mut_counts / summary.dt$count
summary.dt$fraction[summary.dt$count <= 3] <- NA


pdf(file = paste0(sample, "_Heatmap.pdf"), width = 10, height = 10)
print(
ggplot(data = summary.dt, aes(x = lineage, y = mutation, fill = fraction, label = paste0(mut_counts, "/", count))) +
  geom_tile() +
  geom_text(size = 3) + 
  scale_y_discrete(limits = rev(muts)) +
  theme_classic() +
  theme(axis.line=element_blank(), aspect.ratio = nrow(counts.df)/ncol(counts.df),
        plot.title = element_text(hjust = 0.5),
        axis.ticks=element_blank(), panel.background = element_rect(colour = "black", size=1)) +
  ggtitle(sample) +
  scale_fill_gradient(low = "white", high = "#4B0092", limits = c(0, 1))
)
dev.off()

}


