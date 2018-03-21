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
prepareSCRNAData <- function(counts, organism) {
  if(is.data.frame(counts)){
    counts<-as.matrix(counts)
  }
  stopifnot(is.matrix(counts))
  
  totalCount<-sum(counts)
  totalCell<-ncol(counts)
  totalGene<-nrow(counts)

  countInCells<-apply(counts,2,sum)
  rCount<-.getRange(countInCells)

  geneInCells<-apply(counts > 0,2,sum)
  rGene<-.getRange(geneInCells)

  sce <- SingleCellExperiment(list(counts = counts))
  
  #is mitochondrial genes (human 14 mitochondrial genes and mouse 13 mitochondrial genes)
  is.mito <- grepl("^mt-|^MT-", rownames(sce)) 
  
  sce <- calculateQCMetrics(sce, feature_controls = list(Mt = is.mito))
  
  maxMtRNAPercentage<-max(sce$pct_counts_Mt)
  maxRRNAPercentage<-0

  ##remove low quality cells
  libsize.drop <- isOutlier(sce$total_counts, nmads = 3, type = "lower", log = TRUE)
  FCount<-sum(libsize.drop)
  
  feature.drop <- isOutlier(sce$total_features, nmads = 3, type = "lower", log = TRUE)
  FGene<-sum(feature.drop)
  
  FrRNA<-0

  mito.drop <- isOutlier(sce$pct_counts_Mt, nmads = 3, type = "higher")
  FmtRNA<-sum(mito.drop)

  combined.drop <- libsize.drop | feature.drop | mito.drop
  Ftotal<-sum(combined.drop)
  
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
  
  metadata(sce)$filters<-c("Count"=totalCount,
                           "Cell"=totalCell,
                           "Gene"=totalGene,
                           "R-Count"=rCount,
                           "R-Gene"=rGene,
                           "mtRNA"=maxMtRNAPercentage,
                           "rRNA"=maxRRNAPercentage,
                           "F-Count"=FCount,
                           "F-Gene"=FGene,
                           "F-rRNA"=FrRNA,
                           "F-mtRNA"=FmtRNA,
                           "F"=Ftotal)

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
  
  if(!missing(organism)){
    metadata(sce)$hvgPathway <- .getIndividualPathway(var.out, "FDR", organism)
    metadata(sce)$pc1Pathway <- .getIndividualPathway(pc1genes, "adj.P.Val", organism)
  }

  return(list(sce = sce, hvg = var.out, pc1genes = pc1genes, var.fit = var.fit))
}
