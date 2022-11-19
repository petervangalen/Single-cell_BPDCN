# Peter van Galen, 220707
# Compare genotyping efficiency of different approaches

library(tidyverse)
library(Seurat)
library(readxl)

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
as_tibble(subset(seu, orig.ident %in% c("Pt10Dx", "Pt10Rel"))@meta.data) %>% .$RAB9A.3pUTR %>% table
(1843+4591)/(1843+11170+4591) # 0.3655
as_tibble(subset(seu, orig.ident %in% c("Pt10Dx", "Pt10Rel"))@meta.data) %>% .$MTAP.rearr %>% table
(1724+1750)/(1724+14130+1750) # 0.197
as_tibble(subset(seu, orig.ident %in% c("Pt10Dx", "Pt10Rel"))@meta.data) %>% .$`RPS24.chr10:79795273:T/C` %>% table
(678+9326)/(678+102+9326) # 0.9899

# Text: "Investigating 16 founder mutations in uninvolved bone marrow from five BPDCN patients, we detected a total of 10,245 wild-type and 1,204 mutated cells"
f_mut <- filter(genotyping_tables.tib, `Founder or progression mutation` == "Founder") %>% .$Mutation %>% unique
uninvolved_tib <- as_tibble(seu@meta.data) %>% filter(bm_involvement == "No")
uninvolved_tib$orig.ident %>% unique
f_summary <- tibble(n_mut =  select(uninvolved_tib, f_mut) %>% apply(., 1, function(x) sum(grepl("mutant", x))),
                    n_wt = select(uninvolved_tib, f_mut) %>% apply(., 1, function(x) sum(grepl("wildtype", x))))
f_summary %>% filter(n_mut == 0, n_wt > 0)
f_summary %>% filter(n_mut > 0)


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

# Prepare metadata
seu$TET2_expr <- GetAssayData(seu)["TET2",]
metadata_tib <- as_tibble(seu@meta.data)

# Make a list of patients and mutations
tet2_muts_per_sample <- genotyping_tables.tib %>% filter(grepl("TET2", Mutation)) %>% select(Sample, Mutation) %>%
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
    metadata_subset_tib <- metadata_tib %>% filter(orig.ident == pt) %>% select(CellType, all_of(m), TET2_expr)
    
    # Complate list
    names(my_ls)[l] <- paste(m, "in", pt)
    my_ls[[l]] <- metadata_subset_tib %>% group_by(CellType) %>% filter(n() > 10) %>% rename(mut = m) %>%
      summarize(mean_TET2_expr = mean(TET2_expr), genotyping_efficiency = mean(grepl("mutant|wildtype", mut)))
  }
}

# Plot
unlist_my_ls <- do.call(rbind, my_ls)
unlist_my_ls$mut_pt <- rep(names(my_ls), lapply(my_ls, nrow))

pdf(file = "5.4.2_TET2_expr_vs_genotyping.pdf", height = 8, width = 10)
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

