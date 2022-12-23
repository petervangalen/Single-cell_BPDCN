# Peter van Galen, 220707
# Compare genotyping efficiency of different approaches

library(tidyverse)
library(Seurat)
library(readxl)
library(janitor)
library(ggrepel)
library(cowplot)

setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/05_Stats")

rm(list=ls())

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")
donor_colors <- popcol.tib$hex[23:40]
names(donor_colors) <- popcol.tib$pop[23:40]
cell_colors <- popcol.tib$hex[1:21]
names(cell_colors) <- popcol.tib$pop[1:21]

# Load Seurat objects
seurat_files <- list.files("../04_XV-seq", pattern = "*.rds", full.names = T)
seu_ls <- lapply(seurat_files, function(x) readRDS(x))
names(seu_ls) <- cutf(basename(seurat_files), d = "_")
seu <- merge(seu_ls[[1]], seu_ls[2:length(seu_ls)])

# Load genotyping information
genotyping_tables.tib <- read_excel("../04_XV-seq/XV-seq_overview.xlsx")
# Replace different MTAP entries with one, just as in 4.1_Add_GoT-XV_to_Seurat.R
genotyping_tables.tib$Mutation <- gsub("MTAP.rearr.*", "MTAP.rearr", genotyping_tables.tib$Mutation)

# Text: "In total, we enriched 40 mutations from the single-cell libraries"
genotyping_tables.tib$Mutation %>% unique()

# Text: "RAB9A, CDKN2A, and RPS24 in 37%, 20%, and 99% of cells, respectively"
as_tibble(subset(seu, orig.ident %in% c("Pt10Dx", "Pt10Rel"))@meta.data) %>% .$RAB9A.3pUTR %>% tabyl
(1843+4591)/(1843+11170+4591) # 0.3655
as_tibble(subset(seu, orig.ident %in% c("Pt10Dx", "Pt10Rel"))@meta.data) %>% .$MTAP.rearr %>% tabyl
(1724+1750)/(1724+14130+1750) # 0.197
as_tibble(subset(seu, orig.ident %in% c("Pt10Dx", "Pt10Rel"))@meta.data) %>% .$`RPS24.chr10:79795273:T/C` %>% tabyl
(678+9326)/(678+102+9326) # 0.9899

# Text: "Investigating 16 founder mutations in uninvolved bone marrow from five BPDCN patients, we detected a total of 10,245 wild-type and 1,204 mutated cells"
f_mut <- filter(genotyping_tables.tib, `Founder or progression mutation` == "Founder") %>% .$Mutation %>% unique
uninvolved_tib <- as_tibble(seu@meta.data) %>% filter(bm_involvement == "No")
uninvolved_tib$orig.ident %>% unique
f_summary <- tibble(n_mut = dplyr::select(uninvolved_tib, f_mut) %>% apply(., 1, function(x) sum(grepl("mutant", x))),
                    n_wt = dplyr::select(uninvolved_tib, f_mut) %>% apply(., 1, function(x) sum(grepl("wildtype", x))))
f_summary %>% filter(n_mut == 0, n_wt > 0)
f_summary %>% filter(n_mut > 0)

# Cover letter: "The revision includes genotyping for 40 DNA mutations across n=27,994 cells (up from 15 mutations in n=12,613 cells )."
all_muts <- unique(genotyping_tables.tib$Mutation)
all_muts %in% colnames(seu@meta.data)
calls_per_cell <- apply(select(seu@meta.data, all_muts), 1, function(x) sum(grepl("mutant|wildtype", x)))
plot(rev(sort(calls_per_cell)))
sum(calls_per_cell) # 45,576 total transcripts called
sum(calls_per_cell > 0) # 27,994 cells called
# Submission only
submission_muts <- filter(genotyping_tables.tib, `Initial submission or revision` == "Submission") %>% .$Mutation %>% unique
calls_per_cell <- apply(select(seu@meta.data, submission_muts), 1, function(x) sum(grepl("mutant|wildtype", x)))
sum(calls_per_cell) # 21,613 total transcripts called
sum(calls_per_cell > 0) # 12,963 cells called


# Plot genotyping efficiency for all mutations ----------------------------------------------------

pdf("5.4.1_Genotyping_efficiency.pdf", width = 10, height = 6)
genotyping_tables.tib %>% mutate(Mutation = factor(Mutation, levels = unique(Mutation))) %>%
  ggplot(aes(x = Mutation, y = `Genotyping efficiency (%)`, color = Sample)) +
  geom_point(size = 3) +
  scale_color_manual(values = donor_colors[unique(genotyping_tables.tib$Sample)]) +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(aspect.ratio = 0.5, axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(color = "black"))
dev.off()


# Expression vs. genotyping efficiency ------------------------------------------------------------

# Make a new table to assess genotyping efficiency vs. gene expression levels
genotyping_vs_expr_tib <- genotyping_tables.tib %>% dplyr::select(Sample, Mutation, `Genotyping efficiency (%)`) %>%
  mutate(gene = cutf(Mutation, d = "\\."), .before = 3) %>% na.omit %>%
  filter(gene %in% rownames(seu)) %>% mutate(Expression = NA)

# Fill in the expression 
for (n in 1:nrow(genotyping_vs_expr_tib)) {
  #n <- 1
  Sample <- genotyping_vs_expr_tib$Sample[n]
  gene <- genotyping_vs_expr_tib$gene[n]
  genotyping_vs_expr_tib$Expression[n] <- mean(GetAssayData(seu_ls[[Sample]])[gene,])
}

pdf(file = "5.4.2_Genotyping-vs-Expr.pdf", height = 4, width = 5)
genotyping_vs_expr_tib %>%
  ggplot(aes(x = Expression, y = `Genotyping efficiency (%)`, color = Sample, label = Mutation)) +
  geom_point() +
  geom_text_repel(show.legend = F) +
  scale_color_manual(values = donor_colors[unique(genotyping_vs_expr_tib$Sample)]) +
  annotate(geom = "text", x = 0.3, y = 100,
           label = paste("r = ", round(cor(genotyping_vs_expr_tib$Expression, genotyping_vs_expr_tib$`Genotyping efficiency (%)`), 2))) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"))
dev.off()


# Mutated cell fraction vs. VAF from bulk sequencing ----------------------------------------------

# Wrangle RHP targeted sequencing, WES/WGS columns, add mutant UMI fraction
cor_tib <- genotyping_tables.tib %>%
  mutate(`Bulk VAF (RHP)` = as.numeric(`Bulk VAF (RHP)`)) %>% 
  mutate(`Bulk VAF (WES/WGS)` = as.numeric(gsub("nd", "0", `Bulk VAF (WES/WGS)`))) %>%
  mutate(`Mutant UMI fraction (%)` = `Number of mutant UMIs` / (`Number of wildtype UMIs`+`Number of mutant UMIs`) * 100) %>%
  dplyr::select(Sample, Mutation, `Bulk VAF (RHP)`, `Bulk VAF (WES/WGS)`, `Mutant cell fraction (%)`, `Mutant UMI fraction (%)`)

# Generate four plots
plot1_tib <- cor_tib %>% filter(! is.na(`Bulk VAF (RHP)`))
p1 <- plot1_tib %>% 
  ggplot(aes(x = `Mutant cell fraction (%)`, y = `Bulk VAF (RHP)`, color = Sample, label = Mutation)) +
  geom_point() +
  geom_text_repel(show.legend = F, size = 2, max.overlaps = 10) +
  scale_color_manual(values = donor_colors[unique(plot1_tib$Sample)]) +
  coord_cartesian(xlim = c(0,100), ylim = c(0,100)) +
  annotate(geom = "text", x = 10, y = 100,
           label = paste("r = ", round(cor(plot1_tib$`Bulk VAF (RHP)`, plot1_tib$`Mutant cell fraction (%)`), 2))) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"))

p2 <- plot1_tib %>% 
  ggplot(aes(x = `Mutant UMI fraction (%)`, y = `Bulk VAF (RHP)`, color = Sample, label = Mutation)) +
  geom_point() +
  geom_text_repel(show.legend = F, size = 2, max.overlaps = 10) +
  scale_color_manual(values = donor_colors[unique(plot1_tib$Sample)]) +
  coord_cartesian(xlim = c(0,100), ylim = c(0,100)) +
  annotate(geom = "text", x = 10, y = 100,
           label = paste("r = ", round(cor(plot1_tib$`Bulk VAF (RHP)`, plot1_tib$`Mutant UMI fraction (%)`), 2))) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"))

plot2_tib <- cor_tib %>% filter(! is.na(`Bulk VAF (WES/WGS)`))
p3 <- plot2_tib %>%
  ggplot(aes(x = `Mutant cell fraction (%)`, y = `Bulk VAF (WES/WGS)`, color = Sample, label = Mutation)) +
  geom_point() +
  geom_text_repel(show.legend = F, size = 2, max.overlaps = 10) +
  scale_color_manual(values = donor_colors[unique(plot2_tib$Sample)]) +
  coord_cartesian(xlim = c(0,100), ylim = c(0,100)) +
  annotate(geom = "text", x = 10, y = 100,
           label = paste("r = ", round(cor(plot2_tib$`Bulk VAF (WES/WGS)`, plot2_tib$`Mutant cell fraction (%)`), 2))) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"))

p4 <- plot2_tib %>%
  ggplot(aes(x = `Mutant UMI fraction (%)`, y = `Bulk VAF (WES/WGS)`, color = Sample, label = Mutation)) +
  geom_point() +
  geom_text_repel(show.legend = F, size = 2, max.overlaps = 10) +
  scale_color_manual(values = donor_colors[unique(plot2_tib$Sample)]) +
  coord_cartesian(xlim = c(0,100), ylim = c(0,100)) +
  annotate(geom = "text", x = 10, y = 100,
           label = paste("r = ", round(cor(plot2_tib$`Bulk VAF (WES/WGS)`, plot2_tib$`Mutant UMI fraction (%)`), 2))) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"))

pdf(file = "5.4.3_Bulk-VAF_vs_XV-seq.pdf", height = 8, width = 10)
print( plot_grid(p1, p2, p3, p4) )
dev.off()


# Expression vs. genotyping efficiency; TET2 for rebuttal letter ----------------------------------

# Prepare metadata
seu$TET2_expr <- GetAssayData(seu)["TET2",]
metadata_tib <- as_tibble(seu@meta.data)

# Make a list of patients and mutations
tet2_muts_per_sample <- genotyping_tables.tib %>% filter(grepl("TET2", Mutation)) %>% dplyr::select(Sample, Mutation) %>%
  split(f = .$Sample) %>% lapply(., function(x) x$Mutation)
patients <- rep(names(tet2_muts_per_sample), lengths(tet2_muts_per_sample))

# Fill in a list with expression levels and genotyping efficiency for patients in which we genotyped TET2
my_ls <- vector(mode = "list", length = length(patients))

for (l in seq_along(patients)) {
  #l <- 1
  pt <- patients[l]
  for (m in tet2_muts_per_sample[[pt]]) {
    #m <- "TET2.S792*"
    
    # Subset metdata for current patient and mutation, maintain CellType and TET2 expression
    metadata_subset_tib <- metadata_tib %>% filter(orig.ident == pt) %>% dplyr::select(CellType, all_of(m), TET2_expr)
    
    # Complate list
    names(my_ls)[l] <- paste(m, "in", pt)
    my_ls[[l]] <- metadata_subset_tib %>% group_by(CellType) %>% filter(n() > 10) %>% rename(mut = m) %>%
      summarize(mean_TET2_expr = mean(TET2_expr), genotyping_efficiency = mean(grepl("mutant|wildtype", mut)))
  }
}

# Plot
unlist_my_ls <- do.call(rbind, my_ls)
unlist_my_ls$mut_pt <- rep(names(my_ls), lapply(my_ls, nrow))

pdf(file = "5.4.4_TET2_expr_vs_genotyping.pdf", height = 8, width = 10)
unlist_my_ls %>%
  mutate(CellType = factor(CellType, levels = levels(seu_ls[[1]]$CellType))) %>%
  mutate(pt = cutf(mut_pt, d = " ", f = 3)) %>%
  mutate(pt = factor(pt, levels = intersect(names(donor_colors), pt))) %>%
  arrange(pt) %>%
  mutate(mut_pt = factor(mut_pt, levels  = unique(mut_pt))) %>%
  ggplot(aes(x = genotyping_efficiency, y = mean_TET2_expr, color = CellType)) +
  geom_point() +
  scale_color_manual(values = cell_colors) +
  facet_wrap(~ mut_pt, scales = "free") +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"))
dev.off()








# Scripts below did not make it into the paper or letter

# Plot genotyping efficiency depending on how the mutation was detected ---------------------------

# Load and subset genotyping efficiency data
subset_tib <- genotyping_tables.tib %>% dplyr::select(Mutation, `Genotyping efficiency (%)`, `Mutation identification`) %>%
  na.omit() %>% mutate(`Mutation identification` = gsub("WES|WGS", "WES/WGS", `Mutation identification`))

# What's the mean?
subset_tib %>% group_by(`Mutation identification`) %>% summarize(mean_efficiency = mean(`Genotyping efficiency (%)`))
subset_tib %>% group_by(`Mutation identification`) %>% summarize(mean_efficiency = median(`Genotyping efficiency (%)`))

pdf("5.4_Genotyping_efficiency.pdf")
subset_tib %>% 
  ggplot(aes(x = `Mutation identification`, y = `Genotyping efficiency (%)`, label = Mutation)) +
  geom_boxplot() +
  geom_point() +
  geom_label(alpha = 0.8, size = 3) +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(aspect.ratio = 1, panel.grid = element_blank())
dev.off()

