# Peter van Galen, 220207
# Make dotplots and heatmaps of cell type marker genes

library(tidyverse)
library(Seurat)
library(readxl)

setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/Github/5_Dotplots")

rm(list=ls())

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
#popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")
#cell_colors <- popcol.tib$hex[1:17]
#names(cell_colors) <- popcol.tib$pop[1:17]
heatcol <- read_excel("../Single-cell_BPDCN_colors.xlsx", sheet = 2, col_names = F) %>% pull("...1")

# Load Seurat objects, merge
bm <- readRDS("../2_Annotate/BM_Seurat_CellTypes.rds")
seurat_files <- list.files("../3_RandomForest", pattern = "*.rds", full.names = T)
bpdcn_seu.ls <- lapply(seurat_files, function(x) readRDS(x))
merge.seu <- merge(x = bm, y = bpdcn_seu.ls)
merge.seu$CellType <- factor(merge.seu$CellType, levels = levels(bm$CellType))
merge.seu@active.ident <- merge.seu$CellType

# Dot Plots ---------------------------------------------------------------------------------------

# Population signature genes
markerGenes.tib <- read_excel("../2_Annotate/markerGenes_select_two.xlsx") %>%
  select(! c(ProB, PreB, Plasma))
showgenes <- do.call(c, markerGenes.tib[1:2,])

# Matrix for gene expression heatmap
plot.mat <- as.matrix( GetAssayData(merge.seu, slot = "data") )[rev(showgenes),]
plot.mat <- plot.mat - rowMeans(plot.mat)

# Plot that have >= 3 cells in all samples
ncells.tib <- as_tibble(merge.seu@meta.data, rownames = "cell") %>% group_by(orig.ident, CellType) %>%
  summarize(n = n()) %>% pivot_wider(id_cols = orig.ident, names_from = CellType, values_from = n)
ncells.df <- data.frame(ncells.tib, row.names = 1)
cell_types <- na.omit( colnames(ncells.df)[apply(ncells.df, 2, function(x) all(x >= 2))] )

for (s in unique(merge.seu@meta.data$orig.ident)) {
#s <- "BM"
print(s)

# Subset Seurat object by cells from current sample
current.seu <- subset(merge.seu, subset = orig.ident == s)

# Subset gene expression data and split by cell type (note that rowMeans were subtracted on the merged data)
plot.s.mat <- plot.mat[,colnames(current.seu)]
plot.s.mat.ls <- lapply(split(colnames(plot.s.mat), f = current.seu@active.ident), function(x) plot.s.mat[,x])
# Exclude cell types with very few cells
plot.s.mat.ls <- plot.s.mat.ls[cell_types]

# Make data frame with expression level and remove extremes
expr.df <- do.call(rbind, lapply(plot.s.mat.ls, rowMeans))
colnames(expr.df) <- paste0("e_", colnames(expr.df))
z.lim <- c(-1, 4)
expr.df[expr.df < z.lim[1]] <- z.lim[1]
expr.df[expr.df > z.lim[2]] <- z.lim[2]

# Make data frame with fraction positive cells & combine data table for ggplot
frac.df <- do.call(rbind, lapply(plot.s.mat.ls, function(x) apply(x, 1, function(y) return(sum(y > 0) / length(y)))))
colnames(frac.df) <- paste0("f_", colnames(frac.df))
all.dt <- data.table(cbind(data.frame(CellType = factor(names(plot.s.mat.ls), levels = names(plot.s.mat.ls))), expr.df, frac.df))

# Melt / prepare for ggplot
melt.dt <- melt.data.table(all.dt, id.vars = "CellType", measure = patterns("^e_", "^f_"))
colnames(melt.dt) <- c("CellType", "Gene", "Expression", "Fraction")
levels(melt.dt$Gene) <- rev(showgenes)

pdf(paste0(s, "_dotplot.pdf"), width = 6, height = 6)
print(
ggplot(melt.dt, aes(x = CellType, y = Gene, color = Expression, size = Fraction)) +
    geom_point(alpha = 0.8) +
    scale_colour_gradientn(limits = z.lim, colors = c("#000000", heatcol[c(2,8,9,10,11)])) +
    scale_radius(range = c(0.1, 5), trans = "exp") +
    theme(aspect.ratio = 2, panel.border=element_rect(colour = "black", fill=NA, size=0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), plot.title = element_text(hjust = 0.5)) +
    RotatedAxis() +
    ggtitle(label = s)
)
dev.off()

}



