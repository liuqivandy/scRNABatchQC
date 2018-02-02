library(SingleCellExperiment)
library(scater)
library(scran)
library(pheatmap)
library(cluster)
library(limma)
library(dynamicTreeCut)
library(gplots)
library(Rtsne)

scRNAseqMultiQC <- function(count1, count2, count3, ERCC = "spikein", GeneofInterest = "") {

    sce1 <- scRNAseqQC(count1)
    sce2 <- scRNAseqQC(count2)
    sce3 <- scRNAseqQC(count3)
    
    ## should report how many cells and how many genes are filtered out in each sample############
    ##Technical features########################
    ###total_counts,total_features, pt_counts_Mt distribution in each sample##################
    plot.multi.dens( list(sce1$sce$total_counts, sce2$sce$total_counts, sce3$sce$total_counts), xlab = "total_counts" )
    plot.multi.dens( list(sce1$sce$total_features, sce2$sce$total_features, sce3$sce$total_features), xlab = "total_genes" )
    plot.multi.dens( list(sce1$sce$pct_counts_Mt, sce2$sce$pct_counts_Mt, sce3$sce$pct_counts_Mt), xlab = "pct_counts_Mt")
    
    ### top 500 genes count distribution # cannot define color
    scecombine <- SingleCellExperiment(list(counts = cbind(rowData(sce1$sce)$ave.count, rowData(sce2$sce)$ave.count, rowData(sce3$sce)$ave.count)))
    colData(scecombine)$condition <- c(1, 2, 3)
    plotScater(scecombine, colour_by = "condition")  ##scater
    ##  averge count vs. detection rate
    #   plot(log10(rowData(sce1$sce)$ave.count),rowData(sce1$sce)$num.cells,pch=20, xlab="Ave counts",ylab="Num of cells",col=1 )
    #   points(log10(rowData(sce2$sce)$ave.count),rowData(sce2$sce)$num.cells,pch=20, xlab="Ave counts",ylab="Num of cells",col=2)
    #   points(log10(rowData(sce3$sce)$ave.count),rowData(sce3$sce)$num.cells,pch=20, xlab="Ave counts",ylab="Num of cells",col=3)
    
    plot(log10(rowData(sce1$sce)$ave.count), rowData(sce1$sce)$num.cells / dim(sce1$sce)[2], 
      pch = 20, xlab = "Ave counts of genes", ylab = "Detecting Rate", col = 1)
    points(log10(rowData(sce2$sce)$ave.count), rowData(sce2$sce)$num.cells / dim(sce2$sce)[2], 
      pch = 20, xlab = "Ave counts of genes", ylab = "Detecting Rate", col = 2)
    points(log10(rowData(sce3$sce)$ave.count), rowData(sce3$sce)$num.cells / dim(sce3$sce)[2],
      pch = 20, xlab = "Ave counts of genes", ylab = "Detecting Rate", col = 3)
    
    ##variance trend
    plot(sce1$hvg$mean, sce1$hvg$total, pch = 16, cex = 0.6,
      xlab = "Mean log-expression", ylab = "Variance of log-expression", col = 1)
    curve(sce1$var.fit$trend(x), col = 1, lwd = 2, add = TRUE)
    
    points(sce2$hvg$mean, sce2$hvg$total, pch = 16, cex = 0.6, 
      xlab = "Mean log-expression", ylab = "Variance of log-expression", col = 2)
    curve(sce2$var.fit$trend(x), col = 2, lwd = 2, add = TRUE)
    points(sce3$hvg$mean, sce3$hvg$total, pch = 16, cex = 0.6,
      xlab = "Mean log-expression", ylab = "Variance of log-expression", col = 3)
    curve(sce3$var.fit$trend(x), col = 3, lwd = 2, add = TRUE)
    
    ###explain variables ###combine different sample into one figure###########
    plotExplanatoryVariables(sce1$sce, variables = c("log10_total_counts")) ##scater
    plotExplanatoryVariables(sce2$sce, variables = c("log10_total_counts"))
    plotExplanatoryVariables(sce3$sce, variables = c("log10_total_counts"))
    
    plotExplanatoryVariables(sce1$sce, variables = c("log10_total_features"))
    plotExplanatoryVariables(sce1$sce, variables = c("log10_total_counts_Mt"))
    
    ###########Global similarity########################################ave count similarity, pca plot, tsne plot ################
    ave.count <- assay(scecombine)
    pairs(log10(ave.count), xaxt = "n", yaxt = "n", upper.panel = panel.cor, gap = 0, lower.panel = panel.dot)
    
    ####construct a combined sce; sceall
    sceall <- SingleCellExperiment(list(counts = cbind(counts(sce1$sce), counts(sce2$sce), counts(sce3$sce))))
    
    colData(sceall)$condition <- c(rep(1, dim(sce1$sce)[2]), rep(2, dim(sce2$sce)[2]), rep(3, dim(sce3$sce)[2]))
    
    ave.counts <- calcAverage(sceall)
    high.ave <- ave.counts >= 0.1
    clusters <- quickCluster(sceall, subset.row = high.ave, method = "igraph")
    sceall <- computeSumFactors(sceall, cluster = clusters, subset.row = high.ave, min.mean = NULL)
    sceall <- normalize(sceall)  #add logcounts assay into the object
    
    ####
    sceall <- runPCA(sceall, ncomponents = 10)
    #plot(reducedDim(sceall,"PCA")[,1:2]),col=sceall$condition,pch=20)
    plotPCA(sceall, colour_by = "condition")  ##scater
    plotTSNE(sceall, use_dimred = "PCA", colour_by = "condition", perplexity = 20, rand_seed = 100)  ###scater
    
    
    ###Biological features similarity### high variable genes similarity ####pc1 genes similarity########################################ave
    ###highly variable genes### select the top 50 genes (adjustable) with FDR<0.01
    hvgout1 <- sce1$hvg[sce1$hvg$FDR < 0.01, ]
    hvg1 <- rownames(hvgout1)[order(hvgout1$bio, decreasing = TRUE)][1:min(50, dim(hvgout1)[1])]
    hvgout2 <- sce2$hvg[sce2$hvg$FDR < 0.01, ]
    hvg2 <- rownames(hvgout2)[order(hvgout2$bio, decreasing = TRUE)][1:min(50, dim(hvgout2)[1])]
    hvgout3 <- sce3$hvg[sce3$hvg$FDR < 0.01, ]
    hvg3 <- rownames(hvgout3)[order(hvgout3$bio, decreasing = TRUE)][1:min(50, dim(hvgout3)[1])]
    hvglist <- unique(c(hvg1, hvg2, hvg3))
    matchid <- rownames(sce1$sce) %in% hvglist
    hvgbio <- cbind(sce1$hvg$bio[matchid], sce2$hvg$bio[matchid], sce3$hvg$bio[matchid])
    rownames(hvgbio) <- rownames(sce1$sce)[matchid]
    heatmap.2(hvgbio, margins = c(2, 10), cexRow = 0.5)
    ###add a similar KEGG pathway enrichment heatmap
    
    #####pc1 genes ## select the top50 genes (adjustable) with FDR<0.01
    pcout1 <- sce1$pc1genes[sce1$pc1genes$adj.P.Val < 0.01, ]
    pc1g <- rownames(pcout1)[order(abs(pcout1$logFC), decreasing = TRUE)][1:min(50, dim(pcout1)[1])]
    pcout2 <- sce2$pc1genes[sce2$pc1genes$adj.P.Val < 0.01, ]
    pc2g <- rownames(pcout2)[order(abs(pcout2$logFC), decreasing = TRUE)][1:min(50, dim(pcout2)[1])]
    pcout3 <- sce3$pc1genes[sce3$pc1genes$adj.P.Val < 0.01, ]
    pc3g <- rownames(pcout3)[order(abs(pcout3$logFC), decreasing = TRUE)][1:min(50, dim(pcout3)[1])]
    pcglist <- unique(c(pc1g, pc2g, pc3g))
    matchid <- rownames(sce1$sce) %in% pcglist
    pcgFC <- cbind(sce1$pc1genes$logFC[matchid], sce2$pc1genes$logFC[matchid], sce3$pc1genes$logFC[matchid])
    rownames(pcgFC) <- rownames(sce1$sce)[matchid]
    heatmap.2(abs(pcgFC), margins = c(2, 10), cexRow = 0.5)
    #########Pairwise comparsion########################find the differentially expressed genes between two scRNA-seq data
    
    design <- model.matrix( ~ 0 + as.factor(sceall$condition))
    colnames(design) <- c("s1", "s2", "s3")
    
    fit <- lmFit(logcounts(sceall), design)
    contrast.matrix <- makeContrasts(s1 - s2, s1 - s3, s2 - s3, levels = design)
    
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)
    s1vss2 <- topTable(fit2, coef = 1, num = dim(sceall)[1], sort.by = "none")
    s1vss3 <- topTable(fit2, coef = 2, num = dim(sceall)[1], sort.by = "none")
    s2vss3 <- topTable(fit2, coef = 3, num = dim(sceall)[1], sort.by = "none")
    diff_s1vss2 <- s1vss2[abs(s1vss2$logFC) > 1 & s1vss2$adj.P.Val < 0.01, ]
    diffg1 <- rownames(diff_s1vss2)[order(abs(diff_s1vss2$logFC), decreasing = TRUE)][1:min(50, dim(diff_s1vss2)[1])]
    diff_s1vss3 <- s1vss3[abs(s1vss3$logFC) > 1 & s1vss3$adj.P.Val < 0.01, ]
    diffg2 <- rownames(diff_s1vss3)[order(abs(diff_s1vss3$logFC), decreasing = TRUE)][1:min(50, dim(diff_s1vss3)[1])]
    diff_s2vss3 <- s2vss3[abs(s2vss3$logFC) > 1 & s2vss3$adj.P.Val < 0.01, ]
    diffg3 <- rownames(diff_s2vss3)[order(abs(diff_s2vss3$logFC), decreasing = TRUE)][1:min(50, dim(diff_s2vss3)[1])]
    diffg <- unique(c(diffg1, diffg2, diffg3))
    matchid <- rownames(s1vss2) %in% diffg
    diffFC <- cbind(s1vss2$logFC[matchid], s1vss3$logFC[matchid], s2vss3$logFC[matchid])
    rownames(diffFC) <- rownames(s1vss2)[matchid]
    colnames(diffFC) <- c("s1vss2", "s1vss3", "s2vss3")
    heatmap.2(diffFC, cexRow = 0.6, cexCol = 0.6, margins = c(5, 10))
  }


####subfunctions#########
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, method = "spearman"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste(prefix, txt, sep = "")

  if (missing(cex.cor)) {
    cex <- 0.8 / strwidth(txt)
  }

  col <- rgb((1:100) / 100, 0, 0)
  rect(0, 0, 1, 1, col = col[round(r * 100)])
  
  text(0.5, 0.5, txt, cex = cex * r, col = "white")
}


panel.dot <- function(x, y, ...) {
  points(x, y, pch = 16)
}



########
plot.multi.dens <- function(s, xlab = "")
{
  junk.x = NULL
  junk.y = NULL
  for (i in 1:length(s)) {
    junk.x = c(junk.x, density(s[[i]])$x)
    junk.y = c(junk.y, density(s[[i]])$y)
  }
  xr <- range(junk.x)
  yr <- range(junk.y)
  plot( density(s[[1]]), xlim = xr, ylim = yr, main = "", xlab = xlab)
  for (i in 1:length(s)) {
    lines(density(s[[i]]), xlim = xr, ylim = yr, col = i)
  }
}

###
###deal with one single RNAseq####################################
scRNAseqQC <- function(count) {
  sce <- SingleCellExperiment(list(counts = (count)))
  #  is mitochondrial genes # human 14 mitochondrial genes #mouse 13 mitochondrial genes
  is.mito <- grepl("^mt-", rownames(sce))
  ##
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
  pc1 <- reducedDim(sce, "PCA")[, 1]
  design <- model.matrix( ~ pc1)
  
  fit <- lmFit(logcounts(sce), design)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  
  pc1genes <- topTable(fit, coef = 2, n = dim(sce)[1], sort.by = "none")

  ##clustering cells into subpopulations###############
  pcs <- reducedDim(sce, "PCA")
  my.dist <- dist(pcs)
  my.tree <- hclust(my.dist, method = "ward.D2")
  
  my.clusters <- unname(cutreeDynamic(my.tree, distM = as.matrix(my.dist), verbose = 0))
  
  sce$cluster <- factor(my.clusters)

  ####check the separatedness of clusters using the silhouette width
  
  clust.col <- scater:::.get_palette("tableau10medium") # hidden scater colours
  sil <- silhouette(my.clusters, dist = my.dist)
  sil.cols <- clust.col[ifelse(sil[, 3] > 0, sil[, 1], sil[, 2])]
  sil.cols <- sil.cols[order(-sil[, 1], sil[, 3])]
  plot(sil, main = paste(length(unique(my.clusters)), "clusters"), border = sil.cols, col = sil.cols, do.col.sort = FALSE)
  return(list(sce = sce, hvg = var.out, pc1genes = pc1genes, var.fit = var.fit))

}
