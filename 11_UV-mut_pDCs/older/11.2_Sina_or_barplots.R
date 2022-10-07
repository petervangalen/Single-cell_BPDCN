# Peter van Galen, 220914
# Visualize UV-associated mutations in relation to cell types or BPDCN scores

library(tidyverse)
library(Seurat)
library(readxl)
library(ggforce)
library(cowplot)

rm(list=ls())

# Set working directory
setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/11_UV-mut_pDCs")

# Functions & colors
source("../Single-cell_BPDCN_Functions.R")
popcol.tib <- read_excel("../Single-cell_BPDCN_colors.xlsx")
cell_colors <- popcol.tib$hex[1:21]
names(cell_colors) <- popcol.tib$pop[1:21]
mut_colors <- popcol.tib$hex[44:46]
names(mut_colors) <- popcol.tib$pop[44:46]

# Load Seurat objects
seurat_files <- list.files("../04_XV-seq", pattern = "*.rds", full.names = T)
seu_ls <- lapply(seurat_files, function(x) readRDS(x))
names(seu_ls) <- cutf(basename(seurat_files), d = "_")
seu <- merge(seu_ls[[1]], seu_ls[2:length(seu_ls)])
seu$CellType <- factor(seu$CellType, levels = levels(seu_ls[[1]]$CellType))

# Add malignant BPDCN cell module score
bpdcn_sign <- read.table("../09_pDC_expr/bpdcn_sign.txt")[,1]
seu <- AddModuleScore(seu, features = list(bpdcn_sign), name = "bpdcn_score")
colnames(seu@meta.data) <- gsub("bpdcn_score1", "bpdcn_score", colnames(seu@meta.data))

# Load genotyping information
genotyping_tables.tib <- read_excel("../04_XV-seq/XV-seq_overview.xlsx")
# Replace different MTAP entries with one, just as in 4.1_Add_GoT-XV_to_Seurat.R
genotyping_tables.tib$Mutation <- gsub("MTAP.rearr.*", "MTAP.rearr", genotyping_tables.tib$Mutation)


# Wrangle metadata --------------------------------------------------------------------------------

metadata_tib <- as_tibble(seu@meta.data, rownames = "cell")

# Select mutations to show. This is based on a variety of considerations. Outlined in my email "Bar plot" on 220918
tctt_mut <- c("HNRNPUL1.F559F", "MAP4K5.P667S", "MALAT1.chr11:65270399:G/A", "ETV6.R369W", "U2AF1.S34F")

# Add genotyping summary
metadata_tib$`TC>TT mutations` <- ifelse(apply(dplyr::select(metadata_tib, all_of(tctt_mut)), 1,
                                     function(x) sum(grepl("wildtype", x) > 0)), yes = "wildtype", no = "no call")
metadata_tib$`TC>TT mutations` <- ifelse(apply(dplyr::select(metadata_tib, all_of(tctt_mut)), 1,
                                     function(x) sum(grepl("mutant", x) > 0)), yes = "mutant", no = metadata_tib$`TC>TT mutations`)

# Factorize & subset relevant information
metadata_tib$`TC>TT mutations` <- factor(metadata_tib$`TC>TT mutations`, levels = c("no call", "wildtype", "mutant"))
metadata_tib <- metadata_tib %>% dplyr::select(orig.ident, CellType, bpdcn_score, all_of(tctt_mut), `TC>TT mutations`)

# Visualize BPDCN scores vs. mutations ------------------------------------------------------------

pdf("11.2.1_TC-TT_SinaPlot.pdf", width = 3.5, height = 3.5)
metadata_tib %>%
  mutate(orig.order = sample(nrow(metadata_tib))) %>% arrange(orig.order) %>%
  filter(`TC>TT mutations` != "no call") %>%
  ggplot(aes(x = paste0("All genotyped cells (n = ", sum(metadata_tib$`TC>TT mutations` != "no call"), ")"),
             y = bpdcn_score, color = `TC>TT mutations`)) +
  geom_sina(aes(group = 1), size = 0.1) +
  scale_color_manual(values = c("#b0c4de", "#ffa500")) +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  theme_bw() +
  theme(aspect.ratio = 2, panel.grid = element_blank(),
        axis.title.x = element_blank(), axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"))
dev.off()


# Visualize UV-associated mutations in pDCs and non-pDCs ------------------------------------------
plot_tib <- metadata_tib %>%
  mutate(n_na = apply(metadata_tib, 1, function(x) sum(is.na(x)))) %>% filter(n_na < 5) %>%
  mutate(pDC_or_not = ifelse(CellType == "pDC", yes = "All pDCs", no = "All other cells")) %>%
  group_by(pDC_or_not, `TC>TT mutations`) %>% count %>% ungroup %>% group_by(pDC_or_not) %>%
  mutate(`% of cells` = n/sum(n)*100) %>% 
  mutate(pDC_or_not = paste0(pDC_or_not, "\n(n = ", formatC(sum(n), big.mark=","), ")"))
  
pdf("11.2.2_TC-TT_All_pDC_vs_All_other.pdf", width = 5, height = 6)
plot_tib %>% filter(`TC>TT mutations` != "no call") %>%
  ggplot(aes(x = `TC>TT mutations`, y = `% of cells`, fill = `TC>TT mutations`)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#b0c4de", "#ffa500")) +
  facet_wrap(~ pDC_or_not) +
  theme_bw() +
  theme(aspect.ratio = 3, panel.grid = element_blank(), axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"), axis.title.x = element_blank())
dev.off()


# Same but split by mutation ----------------------------------------------------------------------

# Also start a list for the next section
split_by_mutation_tib_ls <- vector(mode = "list")

# Visualize plots
pdf("11.2.3_Split_by_mutation.pdf", width = 10, height = 6)

for (m in tctt_mut) {
#m <- tctt_mut[1]

# Wrangle and plot
metadata_subset_tib <- metadata_tib %>% filter(get(m) %in% c("wildtype", "mutant", "no call")) %>%
  mutate(call = factor(.[,m,drop=T], levels = c("no call", "wildtype", "mutant"))) %>%
  arrange(call) %>%
  mutate(plot_title = paste0(m, "\nin ", paste(unique(orig.ident), collapse = " and ")))

p1 <- metadata_subset_tib %>% filter(get(m) %in% c("wildtype", "mutant")) %>%
  ggplot(aes(x = paste0("All genotyped cells (n = ", nrow(.), ")"), y = bpdcn_score, color = call)) +
  geom_sina(aes(group = 1), size = 1) +
  scale_color_manual(values = c("#b0c4de", "#ffa500")) +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  ggtitle(label = unique(metadata_subset_tib$plot_title)) +
  theme_bw() +
  theme(aspect.ratio = 2, panel.grid = element_blank(),
        axis.title.x = element_blank(), axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5))

# Wrangle
plot_tib <- metadata_subset_tib %>%
  mutate(pDC_or_not = ifelse(CellType == "pDC", yes = "All pDCs", no = "All other cells")) %>%
  group_by(pDC_or_not, call, plot_title) %>% count %>% ungroup %>% group_by(pDC_or_not) %>%
  mutate(`% of cells` = n/sum(n)*100) %>%
  mutate(pDC_or_not = paste0(pDC_or_not, "\n(n = ", formatC(sum(n), big.mark=","), ")"))

# Plot
p2 <- plot_tib %>% filter(call != "no call") %>%
  ggplot(aes(x = call, y = `% of cells`, fill = call)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#b0c4de", "#ffa500")) +
  ggtitle(label = unique(plot_tib$plot_title)) +
  facet_wrap(~ pDC_or_not) +
  theme_bw() +
  theme(aspect.ratio = 3, panel.grid = element_blank(), axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"), axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))

# Visualize
print(plot_grid(p1, p2))

# Save
split_by_mutation_tib_ls <- c(split_by_mutation_tib_ls, list(plot_tib))
names(split_by_mutation_tib_ls) <- c(names(split_by_mutation_tib_ls)[-length(split_by_mutation_tib_ls)],
                                     unique(plot_tib$plot_title))

}
dev.off()


# Barplot split by mutation but shown together ----------------------------------------------------

pdf("11.2.4_Split_combined.pdf")
do.call(rbind,split_by_mutation_tib_ls) %>% filter(call != "no call") %>%
  ggplot(aes(x = `% of cells`, y = pDC_or_not, fill = call)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#b0c4de", "#ffa500")) +
  facet_wrap(~ plot_title, scales = "free", ncol = 1, strip.position = "right") +
  theme_bw() +
  theme(aspect.ratio = 0.3, panel.grid = element_blank(),
        axis.ticks = element_line(color = "black", size = 0.3),
        axis.text = element_text(color = "black"),
        axis.title.y = element_blank(),
        panel.border = element_blank(), axis.line = element_line(size = 0.3),
        strip.background = element_blank(),
        strip.text.y = element_text(angle = 0, hjust = 0, face = "bold"))
dev.off()







# THE FOLLOWING IS A QUICK CHECK IN RESPONSE TO VOLKER'S QUESTION ON 220919



# Also start a list for the next section
split_by_mutation_tib_ls <- vector(mode = "list")

for (m in tctt_mut) {
  #m <- tctt_mut[1]
  
  metadata_subset_tib <- metadata_tib %>% filter(get(m) %in% c("wildtype", "mutant", "no call")) %>%
    mutate(call = factor(.[,m,drop=T], levels = c("no call", "wildtype", "mutant"))) %>%
    arrange(call)
  
  for (p in unique(metadata_subset_tib$orig.ident)) {
    #p <- "Pt10Dx"
    
    # Wrangle
    metadata_subset_tib2 <- metadata_subset_tib %>% filter(orig.ident == p)
  
    # Wrangle
    plot_tib <- metadata_subset_tib2 %>%
      mutate(pDC_or_not = ifelse(CellType == "pDC", yes = "All pDCs", no = "All other cells")) %>%
      group_by(pDC_or_not, call, orig.ident) %>% count %>% ungroup %>% group_by(pDC_or_not) %>%
      mutate(`% of cells` = n/sum(n)*100) %>%
      mutate(pDC_or_not = paste0(pDC_or_not, "\n(n = ", formatC(sum(n), big.mark=","), ")")) %>%
      mutate(plot_title = paste0(m, "\nin ", p))
  
    # Save
    split_by_mutation_tib_ls <- c(split_by_mutation_tib_ls, list(plot_tib))
    names(split_by_mutation_tib_ls) <- c(names(split_by_mutation_tib_ls)[-length(split_by_mutation_tib_ls)], paste0(m, "\nin ", p))
    }
  }


split_by_mutation_tib_ls


plot2_tib <- do.call(rbind,split_by_mutation_tib_ls) %>% filter(call != "no call")%>%
  mutate(plot_title = factor(plot_title, levels = unique(plot_title)))


plot2_tib %>%
  ggplot(aes(x = `% of cells`, y = pDC_or_not, fill = call)) +
  geom_bar(stat = "identity", ) +
  scale_fill_manual(values = c("#b0c4de", "#ffa500")) +
  facet_wrap(~ plot_title, scales = "free", ncol = 1, strip.position = "right") +
  theme_bw() +
  theme(aspect.ratio = 0.3, panel.grid = element_blank(),
        axis.ticks = element_line(color = "black", size = 0.3),
        axis.text = element_text(color = "black"),
        axis.title.y = element_blank(),
        panel.border = element_blank(), axis.line = element_line(size = 0.3),
        strip.background = element_blank(),
        strip.text.y = element_text(angle = 0, hjust = 0, face = "bold"))





