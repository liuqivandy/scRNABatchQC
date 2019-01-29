#' prepareSCRNAData
#' @description prepare statistics information from one single cell RNAseq data table

#' @param inputfile; the input file name of gene-by-cell count matrix; each file should be delimited, either in tsv or csv format
#' @param sampleRatio float; the ratio of cells sampled from each experiment to examine the expression similarity(default: 1)
#' @param nHVGs integer; the number of highly variable genes (default: 1000)
#' @param sf integer; Scale factor  (default: 10000)
 
#' @return a named list containing scRNA data, 
#'                                 hvggenes:high variance genes, 
#'                                 pc1genes: genes in first principle component,
#'                                

#' @examples
#' library(scRNABatchQC)
#' sce<-prepareSCRNAData(inputfile="data1.csv")
#' @import R.utils ggplot2 gplots limma data.table irlba Rtsne WebGestaltR rmdformats Matrix statmod DT RCurl
#' @importFrom devtools session_info
#' 
#' @export
prepareSCRNAData <- function(inputfile, organism=c("hsapiens","mmusculus"),sampleRatio=1,nHVGs=1000,sf=10000) {
  organism<-match.arg(organism)
  rawdata<-data.frame(fread(inputfile),row.names=1)
  counts<-.tosparse(rawdata)
  
  scdata <- list()
  
  #scdata$rawdata <- counts
  scdata$ngene<-dim(counts)[1]
  scdata$ncell<-dim(counts)[2]
  scdata$total_counts <- Matrix::colSums(counts)
  scdata$total_features <- Matrix::colSums(counts != 0)
  
  
  
  is.mito <- grepl("^mt-|^MT-", rownames(counts))
  
  scdata$total_counts_Mt <- Matrix::colSums(counts[is.mito, ])
  scdata$pct_counts_Mt <- 100 * scdata$total_counts_Mt/scdata$total_counts
  
  is.rRNA<-grepl("^Rp[sl][[:digit:]]|^RP[SL][[:digit:]]",rownames(counts))
  scdata$total_counts_rRNA <- Matrix::colSums(counts[is.rRNA, ])
  
  scdata$pct_counts_rRNA<-100*scdata$total_counts_rRNA/scdata$total_counts
  
  
  ###
  

  
  #
  
  scdata$libsize.drop <- .findOutlier(scdata$total_counts,  log=TRUE,type = "lower")
  ##filter cells with less than 200 genes or 3 mad lower
  scdata$feature.drop <- .findOutlier(scdata$total_features, log=TRUE,type = "lower", lower.limit=2)
  ##filter cells with larger than 20% mtRNA or 3mad higher
  scdata$mito.drop <- .findOutlier(scdata$pct_counts_Mt, type = "higher",upper.limit=20)
  
  is.drop<- (scdata$libsize.drop | scdata$feature.drop | scdata$mito.drop)
  
  scdata$num.cells <- Matrix::rowSums(counts != 0)

  gene.keep <- scdata$num.cells > 0
  

  scdata$data <- counts[gene.keep, !is.drop]
  scdata$log10_total_counts<-log10(scdata$total_counts)[!is.drop]
  scdata$log10_total_features <- log10(scdata$total_features)[!is.drop]
  scdata$log10_total_counts_rRNA <- log10(scdata$total_counts_rRNA+1)[!is.drop]
  scdata$log10_total_counts_Mt <- log10(scdata$total_counts_Mt+1)[!is.drop]
  
  scdata$ave.counts <- Matrix::rowMeans(scdata$data)
  scdata$num.cells<-scdata$num.cells[gene.keep]


##normalize to scale factor (sf), the default is 10000

 lib_size <- sf/scdata$total_counts[!is.drop]
  

rowind<-scdata$data@i+1  
colind<-findInterval(seq(scdata$data@x)-1,scdata$data@p[-1])+1
scdata$data@x<-log2(scdata$data@x*lib_size[colind]+1)

  

#####
  
  
  
  
  scdata$hvg <- .getMeanVarTrend(scdata$data)
  ##explained by feature###
  scdata$genevar_by_counts<-.getVarExplainedbyFeature(scdata,"log10_total_counts")
  scdata$genevar_by_features<-.getVarExplainedbyFeature(scdata,"log10_total_features")
  scdata$genevar_by_Mt<-.getVarExplainedbyFeature(scdata,"log10_total_counts_Mt")
  scdata$genevar_by_rRNA<-.getVarExplainedbyFeature(scdata,"log10_total_counts_rRNA")
  
  
  ##select the top HVGs highly variable genes for the PCA, default is 1000
  hvggenes <-  rownames(scdata$hvg)[order(scdata$hvg$zval,decreasing=T)][1:nHVGs]
  scdata$pca <- prcomp_irlba(t(scdata$data[rownames(scdata$data)%in%hvggenes, , drop = FALSE]), n= 2)
  
  design <- model.matrix( ~ scdata$pca$x[, 1])
  fit <- lmFit(scdata$data, design)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  
  ##select the top 500 pc1 genes
  scdata$pc1genes <- topTable(fit, coef = 2, n = 500)
## enriched pathways in top 1000 hvgs and 500 pc1 genes
  if(!missing(organism)){
    scdata$hvgPathway <- .getIndividualPathway(hvggenes,organism=organism)
    scdata$pc1Pathway <- .getIndividualPathway(rownames(scdata$pc1genes), organism=organism)
  }
  ###only output partial samples to perform comparison between samples
  if (sampleRatio<1){
       nsample=round(dim(scdata$data)[2]*sampleRatio,0)
       sampleind<-sample(1:dim(scdata$data)[2],nsample)
       scdata$data<-scdata$data[,sampleind]
   }
  #######
  return(scdata)
}
