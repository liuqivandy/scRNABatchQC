library(Rtsne)

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
