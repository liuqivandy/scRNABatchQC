##' prepareSCRNAData
##'
##' The function prepare statistics information from single cell RNAseq data table.
##'
##' @param count the count table with first column as gene name
##' @return a named list containing scRNA data, 
##'                                 hvg:high variance genes, 
##'                                 pc1genes: genes in first principle component,
##'                                 varTrend: variance trend fit result
##' @importFrom scran trendVar decomposeVar
##' @importFrom limma lmFit eBayes topTable
##' @importFrom Matrix Matrix
##' @export prepareSCRNAData
##' @examples 
##' #count1 <- as.matrix(read.csv("sample1.csv", header = F, row.names = 1))
##' #sce1 <- prepareSCRNAData(count1)
prepareSCRNAData <- function(counts, organism) {
  if(!is.matrix(counts)){
    counts <- as.matrix(counts)
  }
  stopifnot(is.matrix(counts))
  
  counts <- as(counts, "RsparseMatrix")
  
  scdata <- list()
  
  scdata$rawdata <- counts
  
  scdata$total_counts <- Matrix::colSums(counts)
  scdata$total_features <- Matrix::colSums(counts != 0)
  
  scdata$is.mito <- grepl("^mt-|^MT-", rownames(counts))
  
  scdata$total_counts_Mt <- Matrix::colSums(counts[scdata$is.mito, ])
  scdata$log10_total_counts_Mt <- log10(scdata$total_counts_Mt)
  scdata$pct_counts_Mt <- 100 * Matrix::colSums(counts[scdata$is.mito, ]) / Matrix::colSums(counts)
  
  scdata$libsize.drop <- .findOutlier(scdata$total_counts, nmads = 3, type = "lower", log = TRUE)
  scdata$feature.drop <- .findOutlier(scdata$total_features, nmads = 3, type = "lower", log = TRUE)
  scdata$mito.drop <- .findOutlier(scdata$pct_counts_Mt, nmads = 3, type = "higher")
  
  scdata$counts <- scdata$rawdata[, !(scdata$libsize.drop | scdata$feature.drop | scdata$mito.drop)]
  
  scdata$total_counts <- Matrix::colSums(scdata$counts)
  scdata$log10_total_counts <- log10(scdata$total_counts)
  scdata$total_features <- Matrix::colSums(scdata$counts != 0)
  scdata$log10_total_features <- log10(scdata$total_features)
  
  scdata$is.mito <- grepl("^mt-|^MT-", rownames(scdata$counts))
  
  scdata$total_counts_Mt <- Matrix::colSums(scdata$counts[scdata$is.mito, ])
  scdata$log10_total_counts_Mt <- log10(scdata$total_counts_Mt)

  scdata$lib_size <- scdata$total_counts/mean(scdata$total_counts)
  
  counts_norm_lib_size <- t(apply(scdata$counts, 1, function(x) x/scdata$lib_size ))
  
  scdata$ave.counts <- apply(counts_norm_lib_size, 1, mean)
  
  scdata$num.cells <- Matrix::rowSums(scdata$counts != 0)
  to.keep <- scdata$num.cells > 0
  
  scdata$data <- log2(counts_norm_lib_size[to.keep, ] + 1)
  
  scdata$mean <- apply(scdata$data, 1, mean)
  scdata$var <- apply(scdata$data, 1, var)
  
  scdata$varTrend <- trendVar(scdata$data, parametric = TRUE)
  scdata$trends <- scdata$varTrend$trend(scdata$mean)
  
  scdata$hvg <- decomposeVar(scdata$data, scdata$varTrend)
  
  feature_set <- head(order(scdata$var, decreasing = T), n = 500)
  scdata$pca <- prcomp(t(scdata$data[feature_set, , drop = FALSE]), rank. = 10)
  
  design <- model.matrix( ~ scdata$pca$x[, 1])
  fit <- lmFit(scdata$data, design)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  
  scdata$pc1genes <- topTable(fit, coef = 2, n = dim(scdata$data)[1], sort.by = "none")

  if(!missing(organism)){
    scdata$hvgPathway <- .getIndividualPathway(scdata$hvg, "FDR", organism)
    scdata$pc1Pathway <- .getIndividualPathway(scdata$pc1genes, "adj.P.Val", organism)
  }

  return(scdata)
}
