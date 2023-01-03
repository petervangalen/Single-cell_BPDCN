# Peter van Galen, 220508
# Create Seurat objects for all samples with all genotyping information

library(tidyverse)
library(Seurat)
library(readxl)

rm(list=ls())

# Set working directory
setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/04_XV-seq")

# Load commonly used function
cutf <- function(x, f=1, d="/") sapply(strsplit(x, d), function(i) paste(i[f], collapse=d))

# Load BPDCN Seurat objects
bm <- readRDS("../02_Annotate/BM_Seurat_CellTypes.rds")
seurat_files <- list.files("../03_RandomForest", pattern = "Pt.*.rds", full.names = T)
seu_predict.ls <- lapply(seurat_files, function(x) readRDS(x))
names(seu_predict.ls) <- gsub("_Seurat_Predict.rds", "", cutf(seurat_files, "/", f = 3))

# Location of IronThrone results
genotyping_tables.tib <- read_excel("XV-seq_overview.xlsx")


# Function to add genotyping data from FilteredCells.txt to Seurat object -------------------------

genotyping2seu <- function(Seurat_object, FilteredCells, Mutation) {
    #Seurat_object <- seu_predict.ls[[1]]
    #FilteredCells <- "Pt10Dx/RPS24.40/Patient10_RPS24.40_summTable.txt"
    #Mutation <- "RPS24.chr10:79795273:T/C"
    
    # Read and summarize UMI detection
    UMIs.df <- read.table(FilteredCells, header = T)
    calls.tib <- tibble(CB = UMIs.df$BC, call = factor( ifelse(UMIs.df$mutUMIs > 0, yes = "mutant", no = "wildtype"), levels = c("mutant", "wildtype", "no call")))
    
    # Generate a tibble with the same CBs as the Seurat object and genotyping calls
    genotyping.tib <- tibble(CB = cutf(colnames(Seurat_object), d = "-", f = 1))
    genotyping.tib <- genotyping.tib %>% left_join(calls.tib)
    genotyping.tib <- genotyping.tib %>% mutate(call = replace_na(call, "no call"))
    
    # Any cell barcode that occurred in multiple libraries will show up here as a duplicated cell barcode and cannot be called because we can't know which well it came from
    dups.ch <- genotyping.tib$CB[duplicated(genotyping.tib$CB)]
    genotyping.tib$call[genotyping.tib$CB %in% dups.ch] <- "no call" 
    
    # Check that all the cells are still in the same order
    stopifnot( all( genotyping.tib$CB == cutf(colnames(Seurat_object), d = "-") ) )
    
    # Add metadata to Seurat object
    Seurat_object@meta.data[,Mutation] <- genotyping.tib$call
    
    return(Seurat_object)
}


# Add genotyping data to Seurat objects -----------------------------------------------------------
seu_anno.ls <- seu_predict.ls

for (s in unique(genotyping_tables.tib$Sample)) {
  #s <- "Pt1Dx"
  print(s)
  current_tables_tib <- filter(genotyping_tables.tib, Sample == s)
  for (m in current_tables_tib$Mutation) {
    #m <- current_tables_tib$Mutation[1]
    seu_anno.ls[[s]] <- genotyping2seu(Seurat_object = seu_anno.ls[[s]],
                                       FilteredCells = filter(current_tables_tib, Mutation == m)$FilteredCells,
                                       Mutation = m)
  }
}

# Compare without and with genotyping annotation
seu_predict.ls[[1]]@meta.data %>% head
seu_anno.ls[[1]]@meta.data %>% head


# Special processing for MTAP
any.mutant <- seu_anno.ls$Pt10Dx@meta.data %>% dplyr::select(contains("MTAP")) %>% apply(1, function(x) any(grepl("mutant", x)))
# Wild type is by definition the same for all of them, and here I'm combining all mutant calls:
seu_anno.ls$Pt10Dx@meta.data$MTAP.rearr <- seu_anno.ls$Pt10Dx@meta.data$MTAP.rearr.1
seu_anno.ls$Pt10Dx@meta.data$MTAP.rearr[any.mutant] <- "mutant"
# Remove separate MTAP columns
seu_anno.ls$Pt10Dx@meta.data[,which(grepl("MTAP.rearr.\\d", colnames(seu_anno.ls$Pt10Dx@meta.data)))] <- NULL
# Same for Relapse
any.mutant <- seu_anno.ls$Pt10Rel@meta.data %>% dplyr::select(contains("MTAP")) %>% apply(1, function(x) any(grepl("mutant", x)))
# Wild type is by definition the same for all of them, and here I'm combining all mutant calls:
seu_anno.ls$Pt10Rel@meta.data$MTAP.rearr <- seu_anno.ls$Pt10Rel@meta.data$MTAP.rearr.1
seu_anno.ls$Pt10Rel@meta.data$MTAP.rearr[any.mutant] <- "mutant"
# Remove separate MTAP columns
seu_anno.ls$Pt10Rel@meta.data[,which(grepl("MTAP.rearr.\\d", colnames(seu_anno.ls$Pt10Rel@meta.data)))] <- NULL


# Save BPDCN Seurat objects -----------------------------------------------------------------------

# Save all Seurat objects.
# Add metadata about donor group. Also exclude the "Doublets" level (doublets were already removed in 3.1_RandomForest.R)
lapply(names(seu_anno.ls), function(x) {
  #x <- names(seu_anno.ls)[1]
  if (x %in% c("Pt1Rem", "Pt5Dx", "Pt9Dx", "Pt10Dx", "Pt12Dx")) {
    seu_anno.ls[[x]]$bm_involvement <- "No"
  } else if (x %in% c("Pt1Dx", "Pt10Rel", "Pt12Rel", "Pt14Dx", "Pt15Dx", "Pt16Dx")) {
    seu_anno.ls[[x]]$bm_involvement <- "Yes"
  }
  seu_anno.ls[[x]]$CellType <- factor(seu_anno.ls[[x]]$CellType, levels = setdiff(levels(seu_anno.ls[[x]]$CellType), "Doublets"))
  saveRDS(seu_anno.ls[[x]], file = paste0(x, "_Seurat_Final.rds"))
})


# Save BM Seurat object ---------------------------------------------------------------------------
# Since these are the final Seurat objects I will use for subsequent analyses, I will similarly save the healthy controls
bm_subset <- subset(bm, subset = CellType != "Doublets")
bm_subset$CellType <- factor(bm_subset$CellType, levels = setdiff(levels(bm_subset$CellType), "Doublets"))
bm_subset$Doublets <- NULL
bm_subset$bm_involvement <- "HD"
bm_subset@commands <- list()
saveRDS(bm_subset, file = "BM_Seurat_Final.rds")



