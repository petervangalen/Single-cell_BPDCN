# Peter van Galen, 220706
# Bar plot and t-test of cell type proportions

library(tidyverse)
library(Seurat)
library(readxl)

setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/05_Stats")

rm(list=ls())

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")
cell_colors <- popcol.tib$hex[1:21]
names(cell_colors) <- popcol.tib$pop[1:21]
group_colors <- popcol.tib$hex[41:43]
names(group_colors) <- popcol.tib$pop[41:43]

# Load Seurat objects
seurat_files <- list.files("../04_XV-seq", pattern = "*.rds", full.names = T)
seu_ls <- lapply(seurat_files, function(x) readRDS(x))
names(seu_ls) <- cutf(basename(seurat_files), d = "_")


# Plot cell type proportion barplot ---------------------------------------------------------------

# Similar procedures are used in 3_RandomForest.R
as_tibble(unlist(lapply(seu_ls, function(x) unique(x@meta.data$bm_involvement))), rownames = "sample") %>%
  mutate(value = factor(value, levels = c("HD", "No", "Yes"))) %>% arrange(value)
bm_samples <- c("BM1", "BM2", "BM3", "BM4", "BM5", "BM6")
skin_only_samples <- c("Pt1Rem", "Pt5Dx", "Pt9Dx", "Pt10Dx", "Pt12Dx")
bm_involvement_samples <- c("Pt1Dx", "Pt10Rel", "Pt12Rel", "Pt14Dx", "Pt15Dx", "Pt16Dx")
healthy_donor_freq <- do.call(cbind, lapply(split(seu_ls$BM@meta.data, f = cutf(seu_ls$BM@meta.data$replicate, d = "\\.")), function(x) table(x$CellType)))
skin_only_freq <- do.call(cbind, lapply(seu_ls[skin_only_samples], function(x) table(x$CellType)))
bm_involvement_freq <- do.call(cbind, lapply(seu_ls[bm_involvement_samples], function(x) table(x$CellType)))

# Normalize to 100
all_freq <- cbind(healthy_donor_freq, skin_only_freq, bm_involvement_freq)
freq_norm <- sweep(all_freq, 2, colSums(all_freq), "/")*100

pdf("5.3_CellTypeProportions.pdf", width = 8, height = 6)
par(mar = c(8,4,8,12), xpd = T)

barplot(freq_norm[nrow(freq_norm):1,], col = rev(cell_colors[rownames(freq_norm)]),
        xaxt = "n", ylab = "Population frequency (%)", border = NA)
axis(side = 1, at = seq(1,ncol(freq_norm))*1.2-0.5, labels = colnames(freq_norm), las = 2, srt = 45)
legend(x = ncol(freq_norm)*1.2+0.5, y = 100, legend = rownames(freq_norm), fill = cell_colors[rownames(freq_norm)], bty = "n", border = NA)

dev.off()


# Plot fold change and P values -------------------------------------------------------------------

# Test if there are significant differences in cell type frequencies between the groups
skin_only_sign_tib <- tibble(CellType = rownames(freq_norm),
  FC = apply(freq_norm[,skin_only_samples], 1, mean) / apply(freq_norm[,bm_samples], 1, mean),
  P = apply(freq_norm, 1, function(x) t.test(x[bm_samples], x[skin_only_samples])$p.value))
bm_involvement_sign_tib <- tibble(CellType = rownames(freq_norm),
  FC = apply(freq_norm[,bm_involvement_samples], 1, mean) / apply(freq_norm[,bm_samples], 1, mean),
  P = apply(freq_norm, 1, function(x) t.test(x[bm_samples],x[bm_involvement_samples])$p.value))

pdf("5.3_Cell_type_proportion_fold_change.pdf")

rbind(mutate(skin_only_sign_tib, donor_group = "skin_only"), mutate(bm_involvement_sign_tib, donor_group = "bm_involvement")) %>%
  mutate(CellType = factor(CellType, levels = levels(seu_ls$BM$CellType))) %>%
  ggplot(aes(x = CellType, y = FC, fill = donor_group, color = P < 0.05)) +
  geom_hline(yintercept = 1, linetype="dashed", color = "grey") +
  scale_fill_manual(values = group_colors) +
  geom_point(pch = 21, alpha = 0.5, size = 5, stroke = 1) +
  scale_color_manual(values = c("white", "black")) +
  ylab("Fold change relative to healthy controls") +
  theme_bw() +
  theme(aspect.ratio = 0.4, axis.line = element_line(color = "black"), axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())

dev.off()


# Test for Y chromosome expression ----------------------------------------------------------------

# Make a new list of Seurat objects; split healthy donors
seu <- merge(seu_ls[[1]], seu_ls[2:length(seu_ls)])
seu$orig.ident2 <- ifelse(grepl("BM", seu$orig.ident), yes = cutf(seu$replicate, d = "\\."), no = seu$orig.ident)
seu_ls2 <- SplitObject(seu, split.by = "orig.ident2")

# Count Y genes (Esther Rheinbay 220908)
ygenes <- c("RPS4Y1", "KDM5D", "DDX3Y", "UTY", "USP9Y", "EIF1AY", "ZFY")
mean_counts <- do.call(cbind, lapply(seu_ls2, function(x) rowMeans(GetAssayData(x, slot = "counts")[ygenes,])))
#write_tsv(as_tibble(mean_counts, rownames = "gene"), file = "ycounts.txt")

# Bar plot
mean_tib <- as_tibble(t(mean_counts), rownames = "Sample") 
mean_tib %>% mutate(y_count = RPS4Y1+KDM5D+DDX3Y+UTY+USP9Y+EIF1AY+ZFY) %>%
  ggplot(aes(x = Sample, y = y_count)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(aspect.ratio = 0.5, axis.text.x = element_text(angle = 45, hjust = 1))

# Split by cell type
counts_ls <- lapply(seu_ls2, function(x) t(as.matrix(GetAssayData(x, slot = "counts")[ygenes,])))
counts_ls[[1]][1:3,]
counts_df <- do.call(rbind, counts_ls)
all(colnames(seu) == rownames(counts_df))

all_metadata_tib <- as_tibble(seu@meta.data, rownames = "cell") %>% left_join(as_tibble(counts_df, rownames = "cell")) %>%
  dplyr::select(orig.ident, orig.ident2, CellType, all_of(ygenes)) %>%
  mutate(y_count = RPS4Y1+KDM5D+DDX3Y+UTY+USP9Y+EIF1AY+ZFY) %>%
  mutate(CellType = factor(CellType, levels = levels(seu_ls[[1]]$CellType)))

#pdf("ycounts_by_celltype.pdf", width = 20, height = 12)
all_metadata_tib %>% group_by(orig.ident2, CellType) %>% summarize(y_count_mean = mean(y_count)) %>%
  ggplot(aes(x = CellType, y = y_count_mean)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  facet_wrap(~ orig.ident2, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
#dev.off()




