library(SingleCellExperiment)
library(scater)
library(scran)
library(pheatmap)
library(cluster)
library(limma)
library(dynamicTreeCut)
library(gplots)
library(Rtsne)

makeSCRNAdata <- function(object) {
  stopifnot(methods::is(object, "data.table"))
  ### object is a data.table object, the first column is gene symbols
  
  ### automatically transpose the matrix
  if (sum(grepl("^mt-|^MT-", as.matrix(object[, 1]))) == 0) {
    object <- object[, data.table(t(.SD), keep.rownames = TRUE), .SDcols=-1]
  }

  sce <- SingleCellExperiment(list(counts = (as.matrix(object[, -1]))))
  rownames(sce) <- as.matrix(object[, 1])

  is.mito <- grepl("^mt-|^MT-", rownames(sce)) #  is mitochondrial genes # human 14 mitochondrial genes #mouse 13 mitochondrial genes

  sce <- calculateQCMetrics(sce, feature_controls = list(Mt = is.mito))
  
  ##remove low quality cells
  ##define outlier
  libsize.drop <- isOutlier(sce$total_counts, nmads = 3, type = "lower", log = TRUE)
  feature.drop <- isOutlier(sce$total_features, nmads = 3, type = "lower", log = TRUE)
  mito.drop <- isOutlier(sce$pct_counts_Mt, nmads = 3, type = "higher")
  ##remove
  sce <- sce[, !(libsize.drop | feature.drop | mito.drop)]
  
  ave.counts <- calcAverage(sce)
  rowData(sce)$ave.count <- ave.counts
  #	demo.keep <- ave.counts >= 0
  #	filtered.sce <- sce[demo.keep,]
  
  ##remove genes not expressed in any cells
  num.cells <- nexprs(sce, byrow = TRUE)
  rowData(sce)$num.cells <- num.cells
  #	to.keep <- num.cells > 0
  #	sce <- sce[to.keep,]
  
  high.ave <- rowData(sce)$ave.count >= 0.1
  clusters <- quickCluster(sce, subset.row = high.ave, method = "igraph")
  sce <- computeSumFactors(sce, cluster = clusters, subset.row = high.ave, min.mean = NULL)
  
  ##Heterogenous populations ##quickCluster first#####
  ##Normalize#####
  
  sce <- normalize(sce)  #add logcounts assay into the object
  
  ###Modeling the technical noise ######on normalized data
  #modeling mean-variance relationship
  var.fit <- trendVar(sce, parametric = TRUE, span = 0.2, use.spikes = FALSE)
  var.out <- decomposeVar(sce, var.fit)
  
  ############
  sce <- runPCA(sce, ncomponents = 10)
  
  ##interpreting heterogeneity across PC1
  pcs <- reducedDim(sce, "PCA")
  design <- model.matrix( ~ pcs[, 1])
  
  fit <- lmFit(logcounts(sce), design)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  
  pc1genes <- topTable(fit, coef = 2, n = dim(sce)[1], sort.by = "none")
  
  ##clustering cells into subpopulations###############

  my.dist <- dist(pcs)
  my.tree <- hclust(my.dist, method = "ward.D2")
  
  my.clusters <- unname(cutreeDynamic(my.tree, distM = as.matrix(my.dist), verbose = 0))
  
  sce$cluster <- factor(my.clusters)

  return(list(sce = sce, var.out = var.out, var.fit = var.fit, pcs = pcs, pc1genes = pc1genes, dist = my.dist, tree = my.tree, cluster = my.clusters))
}

