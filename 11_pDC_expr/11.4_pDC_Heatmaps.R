# Peter van Galen, 220722
# Relate gene expression in pDCs and malignant BPDCN cells to their mutational status

library(tidyverse)
library(Seurat)
library(readxl)
library(ComplexHeatmap)
library(circlize)
library(ggforce)

rm(list=ls())

# Set working directory
setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/11_pDC_expr")

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")
donor_colors <- popcol.tib$hex[23:40]
names(donor_colors) <- popcol.tib$pop[23:40]
group_colors <- popcol.tib$hex[41:43]
names(group_colors) <- popcol.tib$pop[41:43]
mut_colors <- popcol.tib$hex[44:46]
names(mut_colors) <- popcol.tib$pop[44:46]

# Load all data
seurat_files <- list.files("../04_XV-seq", pattern = "*.rds", full.names = T)
seu_ls <- lapply(seurat_files, function(x) readRDS(x))
names(seu_ls) <- cutf(basename(seurat_files), d = "_")
seu <- merge(seu_ls[[1]], seu_ls[2:length(seu_ls)])

# Add metadata from 11.3_Classify_malignant_BPDCN.R
MalignantCalls_df <- read.table("11.3_MalignantCalls_Final.txt", header = T, row.names = "cell")
all(rownames(seu@meta.data) == rownames(MalignantCalls_df))
seu$bpdcn_sign_score <- MalignantCalls_df$bpdcn_sign_score
seu$RF_pDC_score <- MalignantCalls_df$RF_pDC_score
seu$is_malignant <- MalignantCalls_df$is_malignant

# Signature
bpdcn_sign <- read.table("bpdcn_sign.txt")[,1]

# Load genotyping information
genotyping_tables.tib <- read_excel("../04_XV-seq/XV-seq_overview.xlsx")
# Replace different MTAP entries with one, just as in 4.1_Add_GoT-XV_to_Seurat.R
genotyping_tables.tib$Mutation <- gsub("MTAP.rearr.*", "MTAP.rearr", genotyping_tables.tib$Mutation)

# New Seurat metadata columns: were any mutations detected?
all_mutations <- unique(genotyping_tables.tib$Mutation)
seu$any_mut <- ifelse(apply(seu@meta.data[,all_mutations], 1, function(x)
  sum(grepl("wildtype|mutant", x) > 0)), yes = "Yes", no = "No")

# Ensure samples are nicely ordered later
seu$orig.ident <- factor(seu$orig.ident, levels = c("BM", "Pt15Dx", "Pt16Dx", "Pt1Dx", "Pt1Rem", "Pt5Dx", "Pt9Dx",
                                                    "Pt10Dx", "Pt10Rel", "Pt12Dx", "Pt12Rel", "Pt14Dx"))


# Subset data to make two heatmaps ----------------------------------------------------------------

# Extract metadata
metadata_tib <- as_tibble(seu@meta.data, rownames = "cell") %>%
  arrange(orig.ident)

# First, subset for a heatmap showing pDCs from healthy donors and those with malignant infiltration
healthy_ids <- metadata_tib %>% filter(is_malignant == "Healthy") %>% .$cell
malignant_ids <- metadata_tib %>% filter(is_malignant == "Malignant") %>%
  filter(bm_involvement == "Yes", any_mut != "No") %>%
  group_by(orig.ident) %>% slice_sample(n = 30) %>% .$cell
pdcs_healthy_malignant <- subset(seu, cells = c(healthy_ids, malignant_ids))

# Then, subset data for for a heatmap showing pDCs of skin-only patients, separated into premalignant and
# presumed rare circulating tumor cells  
premalignant_ids <- metadata_tib %>% filter(is_malignant == "Premalignant") %>%
  filter(any_mut != "No") %>% # comment out this line to show more premalignant pDCs
  .$cell
ctc_ids <- metadata_tib %>% filter(bm_involvement == "No", is_malignant == "Malignant") %>% .$cell
pdcs_skin_subset <- subset(seu, cells = c(premalignant_ids, ctc_ids))

# What mutations to show (for both heatmaps)?
mutation_detection <- as_tibble(merge(pdcs_healthy_malignant, pdcs_skin_subset)@meta.data) %>%
  dplyr::select(unique(na.omit(genotyping_tables.tib$Mutation))) %>%
  apply(., 2, function(x) sum(grepl("mutant|wildtype", x)))
show_f_mut <- intersect(names(mutation_detection)[mutation_detection > 3],
                        filter(genotyping_tables.tib, `Founder or progression mutation` == "Founder")$Mutation)
show_p_mut <- intersect(names(mutation_detection)[mutation_detection > 3],
                        filter(genotyping_tables.tib, `Founder or progression mutation` == "Progression")$Mutation)
show_p_mut

# Load gene expression data, then normalize (like in 07_DEG_Heatmaps)
expr_mat <- as.matrix(GetAssayData(pdcs_healthy_malignant, slot = "data"))[bpdcn_sign,]
expr_mat <- expr_mat - rowMeans(expr_mat)
z.lim <- c(-2, 4)
expr_mat[expr_mat < z.lim[1]] <- z.lim[1]
expr_mat[expr_mat > z.lim[2]] <- z.lim[2]

# Adjust some colors
col_bw <- colorRamp2(breaks = c(min(seu$bpdcn_sign_score), max(seu$bpdcn_sign_score)), colors = c("white", "black"))

# Define annotation objects
hm_anno_df <- pdcs_healthy_malignant@meta.data[,c("orig.ident", "bpdcn_sign_score", show_f_mut, show_p_mut)]
top_anno.ha <- HeatmapAnnotation(Donor = as.character(hm_anno_df$orig.ident),
                                 Score = hm_anno_df$bpdcn_sign_score,
                                 col = list(Donor = donor_colors, Score = col_bw),
                                 annotation_name_gp = gpar(fontsize = 10),
                                 border = T, simple_anno_size = unit(3, "mm"))
bottom_anno.ha <- HeatmapAnnotation(Founder = as.matrix(hm_anno_df[,show_f_mut]),
                                    Prog = as.matrix(hm_anno_df[,show_p_mut]),
                                    col = list(Founder = mut_colors, Prog = mut_colors),
                                    annotation_name_gp = gpar(fontsize = 6),
                                    border = T, na_col = "white", simple_anno_size = unit(1.8, "mm"))

# Create Heatmap object
hm <- Heatmap(as.matrix(expr_mat),
              col = colItay(c(1:11))[3:11],
              cluster_rows = F,
              cluster_columns = F,
              row_names_gp = gpar(fontsize = 6),
              show_column_names = F,
              column_split = factor(ifelse(grepl("BM", hm_anno_df$orig.ident),
                                           yes = "BM", no = "BPDCN"), levels = c("BM", "BPDCN")),
              top_annotation = top_anno.ha,
              bottom_annotation = bottom_anno.ha,
              name = "Expr",
              column_title = "pDCs from healthy donors and\nsamples with known marrow involvement",
              column_title_gp = gpar(fontsize = 10),
              border = T,
              use_raster = T,
              raster_quality = 10)

pdf("11.4.1_Healthy_malignant_pDC_heatmap.pdf", width = 5, height = 6)
print(hm)
dev.off()


# Skin-only pDC heatmap ---------------------------------------------------------------------------

# Load gene expression data, then normalize
expr_mat <- as.matrix(GetAssayData(pdcs_skin_subset, slot = "data"))[bpdcn_sign,]
expr_mat <- expr_mat - rowMeans(expr_mat)
z.lim <- c(-2, 4)
expr_mat[expr_mat < z.lim[1]] <- z.lim[1]
expr_mat[expr_mat > z.lim[2]] <- z.lim[2]

  # Add extra MALAT1 genoytping information -----------------------------------
  # Because of the outsized importance of the MALAT1 mutation, and the low sequencing depth we achieved,
  # I am including mut/wt transcript calls were supported by 1 read (instead of the regular 3 read threshold).
  malat1_ReadThreshold1 <- read_tsv("../04_XV-seq/Pt12Dx/MALAT1.5155/Patient12_MALAT1.5155_summTable_ReadThreshold1.txt")
  pdcs_skin_subset@meta.data$`MALAT1.chr11:65270399:G/A` <- case_when(
    cutf(colnames(pdcs_skin_subset), d = "-") %in% filter(malat1_ReadThreshold1, mutUMIs == 0)$BC ~ "wildtype",
    cutf(colnames(pdcs_skin_subset), d = "-") %in% filter(malat1_ReadThreshold1, mutUMIs > 0)$BC ~ "mutant",
    .default = pdcs_skin_subset$`MALAT1.chr11:65270399:G/A`)
  # ---------------------------------------------------------------------------

# Define annotation objects
hm_anno_df <- pdcs_skin_subset@meta.data[,c("orig.ident", "bpdcn_sign_score", "is_malignant", show_f_mut, show_p_mut)]
top_anno.ha <- HeatmapAnnotation(Donor = as.character(hm_anno_df$orig.ident),
                                 Score = hm_anno_df$bpdcn_sign_score,
                                 col = list(Donor = donor_colors, Score = col_bw),
                                 annotation_name_gp = gpar(fontsize = 10),
                                 border = T, simple_anno_size = unit(3, "mm"))
bottom_anno.ha <- HeatmapAnnotation(Founder = as.matrix(hm_anno_df[,show_f_mut]),
                                    Prog = as.matrix(hm_anno_df[,show_p_mut]),
                                    col = list(Founder = mut_colors, Prog = mut_colors),
                                    annotation_name_gp = gpar(fontsize = 6),
                                    border = T, na_col = "white", simple_anno_size = unit(1.8, "mm"))

# Create Heatmap object
hm <- Heatmap(as.matrix(expr_mat),
              col = colItay(c(1:11))[3:11],
              cluster_rows = F,
              cluster_columns = F,
              row_names_gp = gpar(fontsize = 6),
              show_column_names = F,
              column_split = factor(hm_anno_df$is_malignant, levels = c("Premalignant", "Malignant")),
              top_annotation = top_anno.ha,
              bottom_annotation = bottom_anno.ha,
              name = "Expr",
              column_title = "pDCs from samples without\nknown bone marrow involvement",
              column_title_gp = gpar(fontsize = 10),
              border = T,
              use_raster = T,
              raster_quality = 10)

pdf("11.4.2_Skin_only_pDC_heatmap.pdf", width = 5, height = 6)
print(hm)
dev.off()


# Simplified skin-only pDC heatmap for Keynote ----------------------------------------------------

# Collapse mut/wt calls
hm_anno_df$founder_summary <- ifelse(apply(hm_anno_df[,show_f_mut], 1, function(x)
  sum(grepl("wildtype", x) > 0)), yes = "wildtype", no = "no call")
hm_anno_df$founder_summary <- ifelse(apply(hm_anno_df[,show_f_mut], 1, function(x)
  sum(grepl("mutant", x) > 0)), yes = "mutant", no = hm_anno_df$founder_summary)
hm_anno_df$prog_summary <- ifelse(apply(hm_anno_df[,show_p_mut], 1, function(x)
  sum(grepl("wildtype", x) > 0)), yes = "wildtype", no = "no call")
hm_anno_df$prog_summary <- ifelse(apply(hm_anno_df[,show_p_mut], 1, function(x)
  sum(grepl("mutant", x) > 0)), yes = "mutant", no = hm_anno_df$prog_summary)

# Define annotation objects
mut_colors[3] <- "white"
top_anno.ha <- HeatmapAnnotation(Donor = as.character(hm_anno_df$orig.ident),
                                 Score = hm_anno_df$bpdcn_sign_score,
                                 col = list(Donor = donor_colors, Score = col_bw),
                                 annotation_name_gp = gpar(fontsize = 10),
                                 border = T, simple_anno_size = unit(3, "mm"))
bottom_anno.ha <- HeatmapAnnotation(Founder = as.matrix(hm_anno_df$founder_summary),
                                    Prog = as.matrix(hm_anno_df$prog_summary),
                                    col = list(Founder = mut_colors, Prog = mut_colors),
                                    annotation_name_gp = gpar(fontsize = 6),
                                    border = T, na_col = "white", simple_anno_size = unit(3, "mm"))

# Create Heatmap object
hm <- Heatmap(as.matrix(expr_mat),
              col = colItay(c(1:11))[3:11],
              cluster_rows = F,
              cluster_columns = F,
              row_names_gp = gpar(fontsize = 6),
              show_column_names = F,
              column_split = factor(hm_anno_df$is_malignant, levels = c("Premalignant", "Malignant")),
              top_annotation = top_anno.ha,
              bottom_annotation = bottom_anno.ha,
              name = "Expr",
              column_title = "pDCs from samples without\nknown bone marrow involvement",
              column_title_gp = gpar(fontsize = 10),
              border = T,
              use_raster = T,
              raster_quality = 10)

pdf("11.4.3_Skin_only_pDC_heatmap_simple.pdf", width = 5, height = 6)
print(hm)
dev.off()



