library(ggplot2)

DEFAULT_LINE_SIZE <- 1.5
DEFAULT_POINT_SIZE <- 1

.getColData<-function(sces, feature){
  if(missing(feature)){
    stop("Need to specify feature of .getColData")
  }
  
  fNo <- which(colnames(colData(sces[[1]]$sce)) == feature)
  if(fNo == 0){
    stop(paste0("Feature ", feature, " is not exists in object sces"))
  }
  
  result<-NULL
  for (i in 1:length(sces)) {
    result<-rbind(result, data.frame(Sample=names(sces)[i], Value=colData(sces[[i]]$sce)[, fNo]))
  }
  return(result)
}

.getCbindRowData<-function(sces, feature){
  if(missing(feature)){
    stop("Need to specify feature of .getRowData")
  }
  
  fNo <- which(colnames(rowData(sces[[1]]$sce)) == feature)
  if(fNo == 0){
    stop(paste0("Feature ", feature, " is not exists in object sces"))
  }
  
  result<-NULL
  for (i in 1:length(sces)) {
    result<-cbind(result, rowData(sces[[i]]$sce)[, fNo])
  }
  colnames(result)<-names(sces)
  return(result)
}

plotClusterSeparateness <- function(sce) {
  ####check the separateness of clusters using the silhouette width
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
}


plotDensity <- function(sces, feature, featureLabel="", scolors, size = DEFAULT_LINE_SIZE ) {
  featureData<-.getColData(sces, feature)
  featureLabel=ifelse(featureLabel=="", feature, featureLabel)
  
  if(missing(scolors)){
    scolors =  1:length(sces)
  }
  
  p<-ggplot(featureData, aes(x=Value)) + 
    geom_density(aes(color=Sample), size=size) + 
    scale_colour_manual(values=scolors) +
    xlab(featureLabel) +
    theme_bw()
  return(p)
}

### top 500 genes count distribution
plotGeneCountDistribution <- function(sces, scolors, nfeatures = 500, size = DEFAULT_LINE_SIZE) {
  prop_mat <- c()
  
  for (i in 1:length(sces)) {
    aveCountSum <- sum(rowData(sces[[i]]$sce)$ave.count)
    aveCountProp <- cumsum(sort(rowData(sces[[i]]$sce)$ave.count, decreasing = TRUE)) / aveCountSum
    prop_mat <- cbind(prop_mat, aveCountProp)
    colnames(prop_mat)[i] <- names(sces)[i]
  }
  rownames(prop_mat) <- seq_len(nrow(prop_mat))
  
  prop_to_plot <- reshape2::melt(prop_mat[seq_len(nfeatures), ])
  colnames(prop_to_plot) <- c("Feature", "Sample", "Proportion_Library")
  
  p <- ggplot(prop_to_plot, 
              aes_string(x = "Feature", y = "Proportion_Library", 
                         group = "Sample", colour = "Sample")) +
    geom_line(size=size) + 
    xlab("Number of features") + ylab("Cumulative proportion of library") +
    scale_color_manual(values = scolors)
  
  return(p)
}

####averge count vs. detection rate

plotAveCountVSdetectRate <- function(sces, scolors, size = DEFAULT_POINT_SIZE) {
  avedetect <- data.frame()
  for (i in 1:length(sces)) {
    tmpavedec <- data.frame(avecount = log10(rowData(sces[[i]]$sce)$ave.count), 
                            detectrate = rowData(sces[[i]]$sce)$num.cells / dim(sces[[i]]$sce)[2], 
                            Sample = rep(names(sces)[i], length(log10(rowData(sces[[i]]$sce)$ave.count))))
    avedetect <- rbind(avedetect, tmpavedec)
  }
  
  p <- ggplot(avedetect, aes_string(x = "avecount", y = "detectrate", 
                                    group = "Sample", colour = "Sample")) + 
    geom_point() + xlab("Average count of genes") + ylab("Detecting rate") +
    scale_color_manual(values = scolors)
  return(p)
}

##variance trend
plotVarianceTrend <- function(sces, scolors, pointSize=DEFAULT_POINT_SIZE, lineSize=DEFAULT_LINE_SIZE) {
  vartrend_dat <- data.frame()
  for (i in 1:length(sces)) {
    tmpvartrend <- data.frame(mean = sces[[i]]$hvg$mean, 
                              total = sces[[i]]$hvg$total, 
                              trend = sces[[i]]$var.fit$trend(sces[[i]]$hvg$mean),
                              Sample = rep(names(sces)[i], length(sces[[i]]$hvg$mean)))
    vartrend_dat <- rbind(vartrend_dat, tmpvartrend)
  }
  
  pp <- ggplot(vartrend_dat, aes_string(x = "mean", y = "total", group = "Sample", colour = "Sample")) + geom_point()
  pl <- ggplot(vartrend_dat, aes_string(x = "mean", y = "trend", group = "Sample", colour = "Sample")) + geom_line(alpha = 0.3, size = 1.5)
  p <- ggplot(vartrend_dat) + 
    geom_point(pp$mapping, size=pointSize) + 
    geom_line(pl$mapping, size=lineSize) + 
    scale_color_manual(values = scolors) + 
    xlab("Mean log-expression") + 
    ylab("Variance of log-expression")
  
  return(p)
}

plotMultiSamplesOneExplanatoryVariables <- function(sces, scolors, feature, featureLabel="", size = DEFAULT_LINE_SIZE) {
  if(missing(scolors)){
    scolors =  1:length(sces)
  }
  
  if(missing(feature)){
    stop("Need to specify feature of plotMultiSamplesOneExplanatoryVariables")
  }
  featureLabel=ifelse(featureLabel=="", feature, featureLabel)
  
  pct_var_explained <- c()
  sample <-c()
  for (i in 1:length(sces)) {
    pev <- plotExplanatoryVariables(sces[[i]]$sce, variables = c(feature)) ##scater
    pct_var_explained <- c(pct_var_explained, pev$data$Pct_Var_Explained)
    sample <-c(sample, rep(names(sces)[i], length(pev$data$Pct_Var_Explained)))
  }
  
  dat <- data.frame(Pct_Var_Explained = pct_var_explained,  Sample = sample)
  
  p <- ggplot(dat, aes(x = Pct_Var_Explained, colour = Sample)) + 
    geom_line(stat = "density", size = size, trim = T) + 
    geom_vline(xintercept = 1, linetype = 2) + 
    scale_x_log10(breaks = 10 ^ (-3:2), labels = c(0.001, 0.01, 0.1, 1, 10, 100)) + 
    xlab(paste0("% variance explained (log10-scale)")) + 
    ylab("Density") + 
    ggtitle(featureLabel) + 
    scale_color_manual(values = scolors) +
    coord_cartesian(xlim = c(10 ^ (-3), 100))
  return(p)
}

###########Global similarity########################################ave count similarity, pca plot, tsne plot ################
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

plotSampleSimilarity <- function(sces) {
  aveCount <- rowData(sces[[1]]$sce)$ave.count
  for (i in 2:length(sces)) {
    aveCount <- cbind(aveCount, rowData(sces[[i]]$sce)$ave.count)
  }
  colnames(aveCount) <- names(sces)
  
  pairs(log10(aveCount), xaxt = "n", yaxt = "n", upper.panel = panel.cor, gap = 0, lower.panel = panel.dot)
}

####################### PCA ##############

plotAllPCA <- function(sceall, scolors, size = 2) {
  pcadata <- data.frame(PC1 = reducedDim(sceall, "PCA")[, 1],
                        PC2 = reducedDim(sceall, "PCA")[, 2], 
                        Sample = as.factor(colData(sceall)["condition"][, 1]))
  pcalabs <- attr(reducedDim(sceall), "percentVar")
  
  p_pca <- ggplot(pcadata, aes(x = PC1, y = PC2, label = Sample)) + 
    geom_point(aes(col = Sample), size = size) + 
    xlab(paste0("Componet 1:", round(pcalabs[1] * 100), "% Variance")) + 
    ylab(paste0("Componet 2:", round(pcalabs[2] * 100), "% Variance")) + 
    scale_colour_manual(values = scolors)
  
  return(p_pca)
}

plotAllTSNE <- function(sceall, scolors, size = 2) {
  tsnedata <- data.frame(D1 = tsne_out$Y[, 1],
                         D2 = tsne_out$Y[, 2], 
                         Sample = as.factor(colData(sceall)["condition"][, 1]))
  
  p_tsne <- ggplot(tsnedata, aes(x = D1, y = D2, label = Sample)) + 
    geom_point(aes(col = Sample), size = size) + 
    xlab("Dimension 1") + ylab("Dimension 2") + 
    scale_colour_manual(values = scolors)
  return(p_tsne)
}
