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
prepareSCRNAData <- function(inputfile, organism) {
  
  rawdata<-data.frame(fread(inputfile),row.names=1)
  counts<-.tosparse(rawdata)
  
  scdata <- list()
  
  scdata$rawdata <- counts
  
  scdata$total_counts <- Matrix::colSums(counts)
  scdata$total_features <- Matrix::colSums(counts != 0)
  scdata$log10_total_counts<-log10(scdata$total_counts)
  scdata$log10_total_features <- log10(scdata$total_features)
  
  
  is.mito <- grepl("^mt-|^MT-", rownames(counts))
  
  scdata$total_counts_Mt <- Matrix::colSums(counts[is.mito, ])
  scdata$log10_total_counts_Mt <- log10(scdata$total_counts_Mt+1)
  scdata$pct_counts_Mt <- 100 * scdata$total_counts_Mt/scdata$total_counts
  
  is.rRNA<-grepl("^Rp[sl][[:digit:]]|^RP[SL][[:digit:]]",rownames(counts))
  scdata$total_counts_rRNA <- Matrix::colSums(counts[is.rRNA, ])
  scdata$pct_counts_rRNA<-100*scdata$total_counts_rRNA/scdata$total_counts
  
  scdata$libsize.drop <- .findOutlier(scdata$log10_total_counts,  type = "lower")
  ##filter cells with less than 200 genes or 3 mad lower
  scdata$feature.drop <- .findOutlier(scdata$log10_total_features, type = "lower", lower.limit=2)
  ##filter cells with larger than 20% mtRNA or 3mad higher
  scdata$mito.drop <- .findOutlier(scdata$pct_counts_Mt, type = "higher",upper.limit=20)
  
  is.drop<- (scdata$libsize.drop | scdata$feature.drop | scdata$mito.drop)
  
  scdata$num.cells <- Matrix::rowSums(counts != 0)

  gene.keep <- scdata$num.cells > 0


  scdata$data <- counts[gene.keep, !is.drop]
  
  scdata$ave.counts <- rowMeans(scdata$data)


##normalize to 10000

 lib_size <- 10000/scdata$total_counts[!is.drop]
  

rowind<-scdata$data@i+1  
colind<-findInterval(seq(scdata$data@x)-1,scdata$data@p[-1])+1
scdata$data@x<-log2(scdata$data@x*lib_size[colind]+1)

  

#####
  
  
  
  
  scdata$hvg <- .getMeanVarTrend(scdata$data)
  
  ##select the top 1000 highly variable genes for the PCA
  hvggenes <- rownames(head(scdata$hvg,1000))
  scdata$pca <- prcomp_irlba(t(scdata$data[rownames(scdata$data)%in%hvggenes, , drop = FALSE]), n= 10)
  
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
