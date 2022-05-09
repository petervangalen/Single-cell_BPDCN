# Peter van Galen, 220508
# Create Seurat objects for all samples with all genotyping information

library(tidyverse)
library(Seurat)
library(readxl)

rm(list=ls())

# Set working directory
setwd("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/Github/7_XV-seq")

# Load commonly used function
cutf <- function(x, f=1, d="/") sapply(strsplit(x, d), function(i) paste(i[f], collapse=d))

# Load BPDCN Seurat objects
seurat_files <- list.files("../3_RandomForest", pattern = "*.rds", full.names = T)
seu_predict.ls <- lapply(seurat_files, function(x) readRDS(x))
names(seu_predict.ls) <- gsub("_Seurat_Predict.rds", "", cutf(seurat_files, "/", f = 3))

# Location of IronThrone results
genotyping_tables.tib <- read_excel("FilteredCells_files.xlsx")


# Function to add genotyping data from FilteredCells.txt to Seurat object -------------------------

genotyping2seu <- function(Seurat_object, FilteredCells, Mut) {
    #Seurat_object <- seu_predict.ls[[1]]
    #FilteredCells <- "Pt10Dx/RPS24.40/Patient10_RPS24.40_summTable.txt"
    #Mut <- "RPS24.chr10:79795273:T/C"
    
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
    Seurat_object@meta.data[,Mut] <- genotyping.tib$call
    
    return(Seurat_object)
}


# Add genotyping data to Seurat objects & save ----------------------------------------------------
seu_anno.ls <- seu_predict.ls

for (s in unique(genotyping_tables.tib$Sample)) {
  #s <- "Pt1Dx"
  print(s)
  current_tables_tib <- filter(genotyping_tables.tib, Sample == s)
  for (m in current_tables_tib$Mut) {
    #m <- current_tables_tib$Mut[1]
    seu_anno.ls[[s]] <- genotyping2seu(Seurat_object = seu_anno.ls[[s]],
                                       FilteredCells = filter(current_tables_tib, Mut == m)$FilteredCells,
                                       Mut = m)
  }
}

# Compare without and with genotyping annotation
seu_predict.ls[[1]]@meta.data %>% head
seu_anno.ls[[1]]@meta.data %>% head


# Special processing for MTAP
any.mutant <- seu_anno.ls$Pt10Dx@meta.data %>% select(contains("MTAP")) %>% apply(1, function(x) any(grepl("mutant", x)))
# Wild type is by definition the same for all of them, and here I'm combining all mutant calls:
seu_anno.ls$Pt10Dx@meta.data$MTAP.rearr <- seu_anno.ls$Pt10Dx@meta.data$MTAP.rearr.1
seu_anno.ls$Pt10Dx@meta.data$MTAP.rearr[any.mutant] <- "mutant"
# Remove separate MTAP columns
seu_anno.ls$Pt10Dx@meta.data[,which(grepl("MTAP.rearr.\\d", colnames(seu_anno.ls$Pt10Dx@meta.data)))] <- NULL
# Same for Relapse
any.mutant <- seu_anno.ls$Pt10Rel@meta.data %>% select(contains("MTAP")) %>% apply(1, function(x) any(grepl("mutant", x)))
# Wild type is by definition the same for all of them, and here I'm combining all mutant calls:
seu_anno.ls$Pt10Rel@meta.data$MTAP.rearr <- seu_anno.ls$Pt10Rel@meta.data$MTAP.rearr.1
seu_anno.ls$Pt10Rel@meta.data$MTAP.rearr[any.mutant] <- "mutant"
# Remove separate MTAP columns
seu_anno.ls$Pt10Rel@meta.data[,which(grepl("MTAP.rearr.\\d", colnames(seu_anno.ls$Pt10Rel@meta.data)))] <- NULL


# Save all Seurat objects
lapply(names(seu_anno.ls), function(x) {
  saveRDS(seu_anno.ls[[x]], file = paste0(x, "_Seurat_Anno.rds"))
})





### THE FOLLOWING IS OBSOLETE ###

# Number of transcripts ---------------------------------------------------------------------------

#transcripts.tib <- tibble(Mut = unique(FilteredCells_files.tib$Mut), BPDCN628_wt = 0, BPDCN628_mut = 0,
#                          BPDCN712_wt = 0, BPDCN712_mut = 0, BPDCN712R_wt = 0, BPDCN712R_mut = 0)
#
#for (n in 1:nrow(FilteredCells_files.tib)) {
#    #n <- 1
#    UMIs.df <- read.table(FilteredCells_files.tib[n,]$FilteredCells, header = T)
#    
#    # Filter for unique high-quality cell barcodes
#    seu_bcs <- cutf(colnames(get(FilteredCells_files.tib$Sample[n])), d = "-")
#    seu_bcs <- setdiff(seu_bcs, seu_bcs[duplicated(seu_bcs)])
#    UMIs.df <- UMIs.df[UMIs.df$BC %in% seu_bcs,]
#    
#    mut <- FilteredCells_files.tib[n,]$Mut
#    wt_col <- str_c(FilteredCells_files.tib[n,]$Sample, "_wt")
#    mut_col <- str_c(FilteredCells_files.tib[n,]$Sample, "_mut")
#    
#    transcripts.tib[match(mut, transcripts.tib$Mut), match(wt_col, colnames(transcripts.tib))] <- sum(UMIs.df$wtUMIs)
#    transcripts.tib[match(mut, transcripts.tib$Mut), match(mut_col, colnames(transcripts.tib))] <- sum(UMIs.df$mutUMIs)
#}
#
## Save
#write_tsv(transcripts.tib, file = str_c(dir.ch, "TranscriptNumbers.txt"))

