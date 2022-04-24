# Peter van Galen, 190601
# General functions for scRNA-seq analysis in the Single-cell_BPDCN project
# Based on "~/DropboxPartners/Pipelines/scRNAseq_SeqWell/190601_FunctionsGeneral.R"


# Sync RStudio project to Github repository (added 220423)
#install.packages("gitcreds")
#library(gitcreds)
#gitcreds_set()
# Also see https://stackoverflow.com/questions/15843937/git-push-hangs-after-total-line for `git config --global http.postBuffer 157286400` and https://gist.github.com/nepsilon/156387acf9e1e72d48fa35c4fabef0b4 for `git rebase -i HEAD~X` and `git push --force`


message("cutf()")
cutf <- function(x, f=1, d="/") sapply(strsplit(x, d), function(i) paste(i[f], collapse=d))


message("scoreSignature()")
scoreSignature <- function(CM, signatures, n = 100, CM.mean = NULL, verbose = F) {
  if(verbose) {
    message("cells: ", ncol(CM))
    message("genes: ", nrow(CM))
    message("genes in signature: ", length(signatures))
    message("processing")
  }
  
  # Calculate average gene expression to define a set of control genes with similar levels as signature genes
  # This can be supplied as an argument to save time
  if(is.null(CM.mean)) { CM.mean <- rowMeans(CM) }
  
  # Loop over each gene
  s.score <- colMeans(do.call(rbind, lapply(signatures, function(gene) {
    if(verbose) message(".", appendLF = FALSE)
    gene.n <- names(sort(abs(CM.mean[gene] - CM.mean))[2:(n+1)])  # Find n control genes with most similar average gene expression value
    CM[gene, ] - colMeans(CM[gene.n, ])  # Substract average value of the n control genes
  })))
  if(verbose) message(" done")
  return(s.score)
}


message("scaleMinMax()")
scaleMinMax <- function(x, min=0, max=1, z=NULL, keepwithin=TRUE) {
  if(!is.null(z)) {   # use z-score normalization instead
    min <- mean(x) - z*sd(x)
    max <- mean(x) + z*sd(x)
  }
  x <- (x-min) / (max-min)
  if(keepwithin) {
    x[x < 0] <- 0
    x[x > 1] <- 1
  }
  x
}


message("plotTSNE()")
plotTSNE <- function(x, pch = 16, ...) {
  par(mfrow=c(1,1),mar=c(4,4,4,4))
  plot(x, xlab = NA, ylab = NA, tck = F, yaxt = "n", xaxt = "n", pch = pch, ...)
}


colCustom <- function(x, z=NULL, colors=c("#DDDDDD", "red")) {   # use grey to red as default
  if(is.null(z)) {   # just scale to min and max value
    x <- scaleMinMax(x, min(x, na.rm=TRUE), max(x, na.rm=TRUE))
  } else if(length(z) == 1) {   # zscore
    x <- scaleMinMax(x, z = z, keepwithin = TRUE)
  } else {  # scale to min max as provided
    x <- scaleMinMax(x, min = z[1], max = z[2], keepwithin = TRUE)
  }
  
  m <- is.na(x)
  x[m] <- 0.5
  r <- colorRamp(colors)(x)
  y <- apply(r, 1, function(rr) rgb(rr[1], rr[2], rr[3], maxColorValue = 255))
  y[m] <- NA
  y
}


message("colItay()")
colItay <- function(x, z=NULL) {
  colCustom(x, z=z, colors=rev(c("#660220", "#b01b2f", "#d46151", "#f2a585", "#fcdbc8", "#f7f7f7", "#d2e5ef", "#94c5dd", "#4794c1", "#2668aa", "#083160")))  # color pick from publication
}

