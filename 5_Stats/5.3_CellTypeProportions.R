# Peter van Galen, 220706
# Bar plot and t-test of cell type proportions

library(tidyverse)
library(Seurat)

setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/Github/5_Stats")

rm(list=ls())

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")
cell_colors <- popcol.tib$hex[1:23]
names(cell_colors) <- popcol.tib$pop[1:23]
group_colors <- popcol.tib$hex[42:44]
names(group_colors) <- popcol.tib$pop[42:44]

# Load Seurat objects
seurat_files <- list.files("../4_XV-seq", pattern = "*.rds", full.names = T)
seu.ls <- lapply(seurat_files, function(x) readRDS(x))
names(seu.ls) <- cutf(basename(seurat_files), d = "_")


# Plot cell type proportion barplot ---------------------------------------------------------------

# Similar procedures are used in 3_RandomForest.R
as_tibble(unlist(lapply(seu.ls, function(x) unique(x@meta.data$donor_group))), rownames = "sample") %>%
  mutate(value = factor(value, levels = c("healthy_bm", "skin_only", "bm_involvement"))) %>% arrange(value)
bm_samples <- c("BM1", "BM2", "BM3", "BM4", "BM5", "BM6")
skin_only_samples <- c("Pt1Mrd", "Pt5Dx", "Pt9Dx", "Pt10Dx", "Pt12Dx")
bm_involvement_samples <- c("Pt1Dx", "Pt10Rel", "Pt12Rel", "Pt14Dx", "Pt15Dx", "Pt16Dx")
healthy_donor_freq <- do.call(cbind, lapply(split(seu.ls$BM@meta.data, f = cutf(seu.ls$BM@meta.data$replicate, d = "\\.")), function(x) table(x$CellType)))
skin_only_freq <- do.call(cbind, lapply(seu.ls[skin_only_samples], function(x) table(x$CellType)))
bm_involvement_freq <- do.call(cbind, lapply(seu.ls[bm_involvement_samples], function(x) table(x$CellType)))

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
  mutate(CellType = factor(CellType, levels = levels(seu.ls$BM$CellType))) %>%
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




