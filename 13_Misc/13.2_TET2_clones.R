# Peter van Galen, 221020
# On various occasions, I tested if the TET2 mutated clones were biased towards pDC differentation, as we saw in the mouse experiments, consistent with published work (https://doi.org/10.1016/j.stemcr.2020.02.011)
# However, our data *do not* reveal such a bias.

library(tidyverse)
library(Seurat)

rm(list = ls())

# Set working directory
setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/13_Misc")

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")

# Load Seurat objects
seurat_files <- list.files("../04_XV-seq", pattern = "*.rds", full.names = T)
seu_ls <- lapply(seurat_files, function(x) readRDS(x))
names(seu_ls) <- cutf(basename(seurat_files), d = "_")
seu <- merge(seu_ls[[1]], seu_ls[2:length(seu_ls)])


# pDC proportion in TET2 mutated clones -----------------------------------------------------------

# Extract metadata and add summary column for TET2 mutation status
metadata_tib <- as_tibble(seu@meta.data, rownames = "cell")
tet_mut <- metadata_tib %>% filter(if_any(.cols = contains("TET2"), .fns = ~ grepl("mutant", .))) %>% .$cell
tet_wt <- metadata_tib %>% filter(if_any(.cols = contains("TET2"), .fns = ~ grepl("wildtype", .))) %>% .$cell
metadata_tib$TET2.summary <- ifelse(metadata_tib$cell %in% tet_wt, yes = "wt", no = "no call")
metadata_tib$TET2.summary <- ifelse(metadata_tib$cell %in% tet_mut, yes = "mut", no = metadata_tib$TET2.summary)

# Calculate pDC proportion in cells +/- TET2 mutation call
summary_tib <- metadata_tib %>%
  #filter(CellType %in% c("pDC", "cDC", "ncMono", "Mono", "ProMono", "GMP", "HSPC", "EarlyEry", "LateEry")) %>%
  group_by(orig.ident, TET2.summary) %>%
  summarize(pDC_prop = sum(CellType == "pDC")/n(), n = n(), bm_involvement = unique(bm_involvement)) %>%
  select(orig.ident, pDC_prop, TET2.summary, bm_involvement, n) %>%
  filter(bm_involvement %in% c("HD", "No"))

# Visualize
summary_tib %>% 
  ggplot(aes(x = orig.ident, y = pDC_prop, color = TET2.summary)) + geom_point(size = 3) + coord_cartesian(ylim = c(0,1))


# Quantify pDC priming in Progenitors with or without TET2 mutations ------------------------------

#prog <- subset(seu, CellType %in% c("HSPC", "GMP"))
#prog <- NormalizeData(prog)
#pdc_feats <- read.table("../2_Annotate/markerGenes.txt", header = T)[,"pDC"]
#prog <- AddModuleScore(prog, features = list(pdc_feats), name = "pDC_Score")
#prog@meta.data %>% head

#prog_metdata.tib <- as_tibble(prog@meta.data, rownames = "cell")
#apply(dplyr::select(prog_metdata.tib, contains("TET2")), 2, unique)
#tet2_mut_prog <- prog_metdata.tib %>% filter(if_any(.cols = contains("TET2"), .fns = ~ grepl("mutant", .))) %>% .$cell

#prog_metdata.tib %>% mutate(TET2.mut = cell %in% tet2_mut_prog) %>%
#  ggplot(aes_string(y = "pDC_Score1", x = "orig.ident2", color = "TET2.mut")) + # color = "bm_involvement"
#  geom_sina() +
#  theme_bw() +
#  ylab("pDC signature score") +
#  theme(aspect.ratio = 1,
#        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# This analysis does not show that TET2 mutated progenitors are primed to pDC differentiation.
