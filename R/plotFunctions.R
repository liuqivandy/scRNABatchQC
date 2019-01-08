
DEFAULT_LINE_SIZE <- 1.5
DEFAULT_POINT_SIZE <- 1

### package cluster, dynamicTreeCut
#plotClusterSeparateness <- function(sce, ...) {
  ####check the separateness of clusters using the silhouette width
 # pcs <- sce$pca$x
 # my.dist <- dist(pcs)
  #my.tree <- hclust(my.dist, method = "ward.D2")
  
  #my.clusters <- unname(cutreeDynamic(my.tree, distM = as.matrix(my.dist), verbose = 0))
  
  #sce$cluster <- factor(my.clusters)
  
  ####check the separatedness of clusters using the silhouette width
 # clust.col <- rainbow(10)
 # sil <- silhouette(my.clusters, dist = my.dist)
 # sil.cols <- clust.col[ifelse(sil[, 3] > 0, sil[, 1], sil[, 2])]
 # sil.cols <- sil.cols[order(-sil[, 1], sil[, 3])]
  #plot(sil, main = paste(length(unique(my.clusters)), "clusters"), 
 #      border = sil.cols, col = sil.cols, do.col.sort = FALSE, ...)
#}

plotDensity <- function(sces, feature, featureLabel = "", 
                        scolors = 1:length(sces), lineSize = DEFAULT_LINE_SIZE ) {
  featureData <- .getColData(sces, feature)
  featureLabel <- ifelse(featureLabel == "", feature, featureLabel)

  p <- ggplot(featureData, aes(x = Value)) + 
    stat_density(aes(color = Sample), size = lineSize,geom="line",position="identity") + 
    scale_colour_manual(values = scolors) +
    xlab(featureLabel) + theme_classic()
  
  return(p)
}

### top 500 genes count distribution
plotGeneCountDistribution <- function(sces, scolors = 1:length(sces), 
                                      nfeatures = 500, lineSize = DEFAULT_LINE_SIZE) {
  prop_mat <- c()
  
  for (i in 1:length(sces)) {
    aveCountSum <- sum(sces[[i]]$ave.counts)
    aveCountProp <- cumsum(sort(sces[[i]]$ave.counts, decreasing = TRUE)) / aveCountSum
    prop_mat <- cbind(prop_mat, aveCountProp)
    colnames(prop_mat)[i] <- names(sces)[i]
  }
  rownames(prop_mat) <- seq_len(nrow(prop_mat))
  
  prop_to_plot <- reshape2::melt(prop_mat[seq_len(nfeatures),,drop=FALSE ])
  colnames(prop_to_plot) <- c("Feature", "Sample", "Proportion_Library")
  
  p <- ggplot(prop_to_plot, 
              aes_string(x = "Feature", y = "Proportion_Library", 
                         group = "Sample", colour = "Sample")) +
    geom_line(size = lineSize) + 
    xlab("Number of features") + ylab("Cumulative proportion of library") +
    scale_color_manual(values = scolors) +
    theme_classic()
  
  return(p)
}


#plotAveCountVSNumberOfCells <- function(sces, scolors = 1:length(sces), size = DEFAULT_POINT_SIZE) {
#  avedetect <- data.frame()
#  for (i in 1:length(sces)) {
#    tmpavedec <- data.frame(avecount = log10(sces[[i]]$ave.counts), 
#                            numberOfCells = sces[[i]]$num.cells, 
#                            Sample = rep(names(sces)[i], length(sces[[i]]$ave.counts)))
#    avedetect <- rbind(avedetect, tmpavedec)
 # }
  
 # p <- ggplot(avedetect, aes_string(x = "avecount", y = "numberOfCells", 
 #                                   group = "Sample", colour = "Sample")) + 
 #   geom_point(size=pointSize)  + xlab("log10(Average count of genes)") + ylab("Number of cells") +
 #   scale_color_manual(values = scolors) +
 #   theme_classic()
 # 
 # return(p)
#}

####averge count vs. detection rate

plotAveCountVSdetectRate <- function(sces, scolors = 1:length(sces), lineSize = DEFAULT_LINE_SIZE) {
  avedetect <- data.frame()
  for (i in 1:length(sces)) {
    tmpavedec <- data.frame(avecount = log10(sces[[i]]$ave.counts), 
                            detectrate = sces[[i]]$num.cells / sces[[i]]$ncell, 
                            Sample = rep(names(sces)[i], length(log10(sces[[i]]$ave.counts))))
    avedetect <- rbind(avedetect, tmpavedec)
  }
  
  p <- ggplot(avedetect, aes_string(x = "avecount", y = "detectrate", 
                                    group = "Sample", colour = "Sample")) + 
    geom_smooth(size=lineSize) + xlab("Average count of genes") + ylab("Detection rate") +
    scale_color_manual(values = scolors) +
    theme_classic()
  
  return(p)
}

##variance trend
plotVarianceTrend <- function(sces, scolors = 1:length(sces), 
                              pointSize=DEFAULT_POINT_SIZE, lineSize=DEFAULT_LINE_SIZE) {
  meanvar_dat <- data.frame()

  for (i in 1:length(sces)) {
    tmpmeanvar <- data.frame(mean = sces[[i]]$hvg$mean, 
                              var= sces[[i]]$hvg$var, 
                              
                              Sample = rep(names(sces)[i], length(sces[[i]]$hvg$mean)))
        
    
    meanvar_dat<- rbind(meanvar_dat, tmpmeanvar)
    
  }
  
  #pp <- ggplot(meanvar_dat, aes_string(x = "mean", y = "var", group = "Sample", colour = "Sample")) + geom_point(size=pointSize)+geom_smooth()
  #pl <- ggplot(trend_dat, aes_string(x = "mean", y = "var", group = "Sample", colour = "Sample")) + geom_line(alpha = 0.3, size = 1.5)
  p <- ggplot(meanvar_dat, aes_string(x = "mean", y = "var", group = "Sample", colour = "Sample"))  + 
    geom_point(size=pointSize) + 
    geom_smooth(size=lineSize) + 
    scale_color_manual(values = scolors) + 
    xlab("Mean log-expression") + 
    ylab("Variance of log-expression") +
    theme_classic()
  
  return(p)
}

plotMultiSamplesOneExplanatoryVariables <- function(sces, scolors = 1:length(sces), 
                                                    feature, lineSize = DEFAULT_LINE_SIZE) {
  if(missing(feature)){
    stop("Need to specify feature of plotMultiSamplesOneExplanatoryVariables")
  }

  pct_var_explained <- c()
  sample <- c()
  for (i in 1:length(sces)) {
    tmp_pct_VE <- sces[[i]][[feature]]
    pct_var_explained <- c(pct_var_explained, tmp_pct_VE)
    sample <- c(sample, rep(names(sces)[i], length(tmp_pct_VE)))
  }
  
  dat <- data.frame(Pct_Var_Explained = pct_var_explained,  Sample = sample)
  
  p <- ggplot(dat, aes(x = Pct_Var_Explained, colour = Sample)) + 
    geom_line(stat = "density", size = lineSize, trim = T) + 
    geom_vline(xintercept = 1, linetype = 2) + 
    scale_x_log10(breaks = 10 ^ (-3:2), labels = c(0.001, 0.01, 0.1, 1, 10, 100)) + 
    xlab(paste0("% variance explained (log10-scale)")) + 
    ylab("Density") + 
    scale_color_manual(values = scolors) +
    coord_cartesian(xlim = c(10 ^ (-3), 100)) + theme_classic()
  
  return(p)
}

###########Global similarity#####################################
###ave count similarity, pca plot, tsne plot ################
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
  points(x, y, pch = 16, ...)
}

plotSampleSimilarity <- function(sces, ...) {
  aveCount <- sces[[1]]$ave.counts
  for (i in 2:length(sces)) {
    aveCount <- merge(aveCount, sces[[i]]$ave.counts, by = "row.names")
    rownames(aveCount) <- aveCount[, 1]
    aveCount <- aveCount[, -1]
  }
  colnames(aveCount) <- names(sces)
  
  pairs(log10(aveCount), xaxt = "n", yaxt = "n", upper.panel = panel.cor, gap = 0, lower.panel = panel.dot, ...)
}

####################### PCA ##############

plotAllPCA <- function(pca_tsne_data, scolors = 1:length(sces), pointSize = DEFAULT_POINT_SIZE) {
  pcadata <- data.frame(PC1 = pca_tsne_data$pca$x[, 1], PC2 = pca_tsne_data$pca$x[, 2], Sample = (pca_tsne_data$condition))
  
  eigs <- pca_tsne_data$pca$sdev ^ 2
  pc1pct <- eigs[1] / sum(eigs)
  pc2pct <- eigs[2] / sum(eigs)
  
  p_pca <- ggplot(pcadata, aes(x = PC1, y = PC2, label = Sample)) + 
    geom_point(aes(col = Sample), size = pointSize) + 
    xlab(paste0("PC1(", round(pc1pct * 100), "%)")) + 
    ylab(paste0("PC2(", round(pc2pct * 100), "%)")) + 
    scale_colour_manual(values = scolors) + theme_classic()
  
  return(p_pca)
}

####################### TSNE ##############

plotAllTSNE <- function(pca_tsne_data, scolors = 1:length(sces), pointSize = DEFAULT_POINT_SIZE) {
  tsnedata <- data.frame(D1 = pca_tsne_data$tsne[, 1], D2 = pca_tsne_data$tsne[, 2], Sample = (pca_tsne_data$condition))
  
  p_tsne <- ggplot(tsnedata, aes(x = D1, y = D2, label = Sample)) + 
    geom_point(aes(col = Sample), size = pointSize) + 
    xlab("Dimension 1") + ylab("Dimension 2") + 
    scale_colour_manual(values = scolors) + theme_classic()
  
  return(p_tsne)
}

plotPairwiseDifference <- function(scesall, FDR = 0.01, geneNo = 50, ...) {
  diffFC <- .getDiffGenes(scesall, FDR = FDR, geneNo = geneNo)
  heatmap.2(as.matrix(diffFC$genes), cexRow = 0.6, cexCol = 0.6, ...)
}
