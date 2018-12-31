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
  scdata$log10_total_counts<-log10(scdata$total_counts)
  scdata$log10_total_features <- log10(scdata$total_features)
  
  
  scdata$is.mito <- grepl("^mt-|^MT-", rownames(counts))
  
  scdata$total_counts_Mt <- Matrix::colSums(counts[scdata$is.mito, ])
  scdata$log10_total_counts_Mt <- log10(scdata$total_counts_Mt)
  scdata$pct_counts_Mt <- 100 * scdata$total_counts_Mt/scdata$total_counts
  
  scdata$is.rRNA<-grepl("^Rp[sl][[:digit:]]|^RP[SL][[:digit:]]",rownames(counts))
  scdata$total_counts_rRNA <- Matrix::colSums(counts[scdata$is.rRNA, ])
  scdata$pct_counts_rRNA<-100*scdata$total_counts_rRNA/scdata$total_counts
  
  scdata$libsize.drop <- .findOutlier(scdata$log10_total_counts,  type = "lower")
  ##filter cells with less than 200 genes or 3 mad lower
  scdata$feature.drop <- .findOutlier(scdata$log10_total_features, type = "lower", lower.limit=2)
  ##filter cells with larger than 20% mtRNA or 3mad higher
  scdata$mito.drop <- .findOutlier(scdata$pct_counts_Mt, type = "higher",upper.limit=20)
  
  is.drop<- (scdata$libsize.drop | scdata$feature.drop | scdata$mito.drop)
  
  scdata$counts <- scdata$rawdata[, !is.drop]
  

##normalize to 10000
  scdata$lib_size <- scdata$total_counts/10000
  
  counts_norm_lib_size <- t(apply(scdata$counts, 1, function(x) x/scdata$lib_size ))
  
  scdata$ave.counts <- apply(counts_norm_lib_size, 1, mean)
  
  scdata$num.cells <- Matrix::rowSums(scdata$counts != 0)
  to.keep <- scdata$num.cells > 0
  
  scdata$data <- log2(counts_norm_lib_size[to.keep, ] + 1)
  
  
  
  
  
  scdata$hvg <- .getMeanVarTrend(scdata$data)
  
  ##select the top 1000 highly variable genes for the PCA
  hvggenes <- rownames(head(scdata$hvg,1000))
  scdata$pca <- prcomp(t(scdata$data[rownames(scdata$data)%in%hvggenes, , drop = FALSE]), rank. = 10)
  
  design <- model.matrix( ~ scdata$pca$x[, 1])
  fit <- lmFit(scdata$data, design)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  
  ##select the top 500 pc1 genes
  scdata$pc1genes <- topTable(fit, coef = 2, n = 500)
## enriched pathways in top 1000 hvgs and 500 pc1 genes
  if(!missing(organism)){
    scdata$hvgPathway <- .getIndividualPathway(head(scdata$hvg,1000),organism)
    scdata$pc1Pathway <- .getIndividualPathway(scdata$pc1genes, organism)
  }

  return(scdata)
}
