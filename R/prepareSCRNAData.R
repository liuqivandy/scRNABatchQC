library(SingleCellExperiment)
library(scater)
library(scran)
library(cluster)
library(limma)
library(dynamicTreeCut)
library(Rtsne)

##' prepareSCRNAData
##'
##' The function prepare statistics information from single cell RNAseq data table.
##'
##' @param count the count table with first column as gene name
##' @return a named list containing sce:SingleCellExperiment, 
##'                                 hvg:high variance genes, 
##'                                 pc1genes: genes in first principle component,
##'                                 var.fit: variance trend fit result
##' @importFrom SingleCellExperiment SingleCellExperiment reducedDim
##' @importFrom Scater calculateQCMetrics isOutlier calcAverage nexprs normalize runPCA .get_palette
##' @importFrom Scran quickCluster computeSumFactors trendVar decomposeVar
##' @importFrom limma lmFit eBayes topTable
##' @importFrom dynamicTreeCut cutreeDynamic
##' @importFrom cluster silhouette
##' @export prepareSCRNAData
##' @examples 
##' #count1 <- as.matrix(read.csv("sample1.csv"))
##' #sce1<-prepareSCRNAData(count1)
prepareSCRNAData <- function(count) {
  sce <- SingleCellExperiment(list(counts = count))
  
  #is mitochondrial genes (human 14 mitochondrial genes and mouse 13 mitochondrial genes)
  is.mito <- grepl("^mt-|^MT-", rownames(sce)) 
  
  sce <- calculateQCMetrics(sce, feature_controls = list(Mt = is.mito))

  #remove empty genes and samples
  emptyGene=rowSums(count)==0
  emptySample=colSums(count)==0
  sce <- sce[!emptyGene, !emptySample]
  
  ##remove low quality cells
  libsize.drop <- isOutlier(sce$total_counts, nmads = 3, type = "lower", log = TRUE)
  feature.drop <- isOutlier(sce$total_features, nmads = 3, type = "lower", log = TRUE)
  mito.drop <- isOutlier(sce$pct_counts_Mt, nmads = 3, type = "higher")
  combined.drop <- libsize.drop | feature.drop | mito.drop
  sce <- sce[, !combined.drop]
  
  ##remove genes not expressed in any cells
  ave.counts <- calcAverage(sce)
  rowData(sce)$ave.count <- ave.counts
  num.cells <- nexprs(sce, byrow = TRUE)
  rowData(sce)$num.cells <- num.cells
  lowGene <- num.cells == 0
  sce <- sce[!lowGene, ]
  
  ##quickCluster and normalization
  high.ave <- rowData(sce)$ave.count >= 0.1
  clusters <- quickCluster(sce, subset.row=high.ave, method="igraph")
  sce <- computeSumFactors(sce, cluster=clusters, subset.row=high.ave, min.mean=NULL, positive=TRUE)
  
  sizeFactorZero <- sizeFactors(sce) == 0
  sce <- sce[, !sizeFactorZero]
  
  metadata(sce)$filters<-c(SampleInit=length(emptySample),
                           SampleEmpty=count(emptySample),
                           SampleLibsizeDrop=count(libsize.drop),
                           SampleFeatureDrop=count(feature.drop),
                           SampleMitoDrop=count(mito.drop),
                           SampleCombinedDrop=count(combined.drop),
                           SampleSizeFactorZero=count(sizeFactorZero),
                           GeneInit = length(emptyGene),
                           GeneEmpty=count(emptyGene),
                           GeneMitoCount=count(is.mito),
                           GeneLowExpress=count(lowGene))

  sce <- normalize(sce)
  
  ##modeling the technical noise on normalized data and modeling mean-variance relationship
  var.fit <- trendVar(sce, parametric = TRUE, span = 0.2, use.spikes = FALSE)
  var.out <- decomposeVar(sce, var.fit)
  
  sce <- runPCA(sce,ncomponents = 10) 
  
  ##interpreting heterogeneity across PC1
  pc1 <- reducedDim(sce, "PCA")[,1]
  design <- model.matrix(~pc1)
  
  fit <- lmFit(logcounts(sce), design)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  
  pc1genes <- topTable(fit, coef = 2, n = dim(sce)[1], sort.by = "none")

  return(list(sce = sce, hvg = var.out, pc1genes = pc1genes, var.fit = var.fit))
}

##' preparePCATSNEData
##'
##' The function prepare statistics information from multiple scRNA dataset.
##'
##' @param sces a named list of makeSCRNAdata result
##' @return a sce:SingleCellExperiment data with PCA and TSNE
##' @importFrom SingleCellExperiment SingleCellExperiment reducedDim
##' @importFrom Scater calculateQCMetrics isOutlier calcAverage nexprs normalize runPCA .get_palette
##' @importFrom Scran quickCluster computeSumFactors trendVar decomposeVar
##' @importFrom Rtsne Rtsne
##' @export preparePCATSNEData
##' @examples 
##' #sces <- prepareSCRNADataSet(sampleTable)
##' #sceall <- preparePCATSNEData(sces)
preparePCATSNEData <- function(sces, ncomponents = 10, perplexity = 20) {
  allCount <- counts(sces[[1]]$sce)
  conditions <- rep(names(sces)[1], dim(sces[[1]]$sce)[2])
  colnames(allCount) <- paste0(names(sces)[1], "cell", 1:dim(sces[[1]]$sce)[2])
  
  for (i in 2:length(sces)) {
    allCount <- merge(allCount, counts(sces[[i]]$sce), by = "row.names", all = T)
    rownames(allCount) <- allCount[, 1]
    allCount <- allCount[, -1]
    conditions <- c(conditions, rep(names(sces)[i], dim(sces[[i]]$sce)[2]))
    colnames(allCount)[(ncol(allCount) - dim(sces[[i]]$sce)[2] + 1):ncol(allCount)] <- paste0(names(sces)[i], "cell", 1:dim(sces[[i]]$sce)[2])
  }
  allCount[is.na(allCount)] <- 0
  
  sceall <- SingleCellExperiment(list(counts = (as.matrix(allCount))))
  
  colData(sceall)$condition <- conditions
  
  ave.counts <- calcAverage(sceall)
  high.ave <- ave.counts >= 0.1
  clusters <- quickCluster(sceall, subset.row = high.ave, method = "igraph")
  sceall <- computeSumFactors(sceall, cluster = clusters, subset.row = high.ave, min.mean = NULL)
  sceall <- normalize(sceall)
  
  sceall <- runPCA(sceall, ncomponents = ncomponents)
  
  set.seed(100)

  tsne_out <- Rtsne(reducedDim(sceall, "PCA"), initial_dims = ncol(reducedDim(sceall, "PCA")), pca = FALSE, perplexity = perplexity)
  reducedDim(sceall, "TSNE") <- tsne_out$Y
  
  return(sceall)
}

DEBUG=FALSE
if(DEBUG){
  count<-as.matrix(read.csv("Z:/shengq1/20180214_scRNABatchQC/S1.csv", row.names=1))
  sce<-prepareSCRNAData(count)
}
