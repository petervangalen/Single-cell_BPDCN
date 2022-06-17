# Peter van Galen, 220516
# Make radar plots of the cell type frequency in "clones" defined by TET2 mutations

library(tidyverse)
library(Seurat)
library(readxl)

rm(list=ls())

# Set working directory
setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/Github/10_RadarPlot")

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")

# Function for coord_polar with straight lines (see https://github.com/petervangalen/MAESTER-2021/blob/main/4_CH_sample/4.4_LineageBias.R)
coord_radar <- function (theta = "x", start = 0, direction = 1) {
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x") "y" else "x"
  ggproto("CoordRadar", CoordPolar, theta = theta, r = r, start = start, 
          direction = sign(direction),
          is_linear = function(coord) TRUE)
}

# Load Seurat objects
seurat_files <- list.files("../7_XV-seq", pattern = "*.rds", full.names = T)
seu_bpdcn.ls <- lapply(seurat_files, function(x) readRDS(x))
names(seu_bpdcn.ls) <- gsub("_Seurat_Anno.rds", "", cutf(seurat_files, d = "/", f = 3))

# Check TET2 expression
bm <- readRDS("../2_Annotate/BM_Seurat_CellTypes.rds")
FeaturePlot(bm, "TET2")
VlnPlot(bm, "TET2")
bm@active.ident %>% head



# Loop over each sample and save genotyping count tables ------------------------------------------

for (Sample in names(seu_bpdcn.ls)) {
#Sample <- "Pt10Dx"
#Sample <- "Pt10Rel"
#Sample <- "Pt12Dx"

# Select Seurat object and TET2 mutations therein. All samples have TET2 mutations; otherwise this would crash
seu <- seu_bpdcn.ls[[Sample]]
TET2.Muts <- colnames(seu@meta.data)[grepl("TET2", colnames(seu@meta.data))]
# Exclude mutations with low detection frequency
TET2.Muts <- TET2.Muts[colSums((seu@meta.data[,TET2.Muts,drop=F] == "mutant")) > 10]

# Select necessary information
metadata_tib <- as_tibble(seu@meta.data, rownames = "cell") %>%
  select(CellType, all_of(TET2.Muts))

# Determine total cell numbers
celltype_numbers_tib <- metadata_tib %>% group_by(CellType) %>% dplyr::count(name = "all_cells") %>% ungroup

# Add mutated cell counts
for (Mut in TET2.Muts) {
  #Mut <- TET2.Muts[1]
  current_freq <- metadata_tib %>% filter(get(Mut) == "mutant") %>% group_by(CellType) %>% dplyr::count(name = Mut) %>% ungroup
  celltype_numbers_tib <- celltype_numbers_tib %>% left_join(current_freq, by = "CellType")
}

# Replace NA with 0
celltype_numbers_tib <- celltype_numbers_tib %>% replace(is.na(.), 0)

# Add frequency and normalized frequency
cell_type_freq.tib <- celltype_numbers_tib %>%
  mutate(across(-CellType, .fns = list(freq = ~ .x / sum(.x) * 100,
                                       norm = ~ (.x / sum(.x)) / (all_cells / sum(all_cells)))))
#cell_type_freq.tib %>% view


# Plot normalized frequency
radar.tib <- cell_type_freq.tib %>% dplyr::select(CellType, contains("norm")) %>%
  add_row(CellType = "empty", .before = 1) %>%
  replace(is.na(.), 0) %>%
  pivot_longer(cols = -CellType) %>%
  mutate(name = factor(gsub("_norm", "", name), levels = unique(gsub("_norm", "", name)))) %>%
  mutate(CellType = factor(CellType, levels = unique(CellType)))

# Prepare colors
#clone_colors <- mycol.ch[intersect(names(mycol.ch), gsub("_norm", "", levels(radar.tib$name)))]
#clone_colors <- c(set_names(clone_colors, str_c(names(clone_colors), "_norm")), all_cells_norm = "black")

pdf(paste0(Sample, "_RadarPlot.pdf"))
print(
radar.tib  %>% mutate(value = log2(value)) %>%
  ggplot(aes(x = as.numeric(CellType), y = value, fill = name, color = name)) +
  coord_radar() +
  geom_line(size = 2) +
  geom_polygon(alpha = 0.2) +
#  scale_color_manual(values = clone_colors) +
#  scale_fill_manual(values = clone_colors) +
  scale_x_continuous(breaks = 1:length(levels(radar.tib$CellType)),
                     labels = c(setdiff(levels(radar.tib$CellType), "empty"), "")) +
  theme_bw() +
  xlab("") +
  ylab("Normalized cell type frequency (log2)") +
  ggtitle(Sample) +
  theme(aspect.ratio = 1,
        panel.border = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "#D3D3D3"),
        axis.text = element_text(size = 14, color = "black"), axis.title.y = element_text(size = 14, color = "black"),
        axis.ticks.y = element_line(color = "black"))
)
dev.off()

}

















