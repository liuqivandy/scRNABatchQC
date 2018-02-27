library(ggplot2)
library(reshape2)
library(WebGestaltR)

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
    scale_color_manual(values = scolors) +
    theme_classic()
  
  return(p)
}


plotAveCountVSNumberOfCells <- function(sces, scolors, size = DEFAULT_POINT_SIZE) {
  avedetect <- data.frame()
  for (i in 1:length(sces)) {
    tmpavedec <- data.frame(avecount = log10(rowData(sces[[i]]$sce)$ave.count), 
                            numberOfCells = rowData(sces[[i]]$sce)$num.cells, 
                            Sample = rep(names(sces)[i], length(rowData(sces[[i]]$sce)$ave.count)))
    avedetect <- rbind(avedetect, tmpavedec)
  }
  
  p <- ggplot(avedetect, aes_string(x = "avecount", y = "numberOfCells", 
                                    group = "Sample", colour = "Sample")) + 
    geom_point() + xlab("log10(Average count of genes)") + ylab("Number of cells") +
    scale_color_manual(values = scolors) +
    theme_classic()
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
    scale_color_manual(values = scolors) +
    theme_classic()
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
    ylab("Variance of log-expression") +
    theme_classic()
  
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
    coord_cartesian(xlim = c(10 ^ (-3), 100)) + 
    theme_classic()
  return(p)
}

### Biological features similarity
### select the top 50 genes (adjustable) with FDR<0.01
### plotBiologicalSimilarity(sces, objectName="hvg", filterName="FDR", valueName="bio")
### plotBiologicalSimilarity(sces, objectName="pc1genes", filterName="adj.P.Val", valueName="logFC")
plotBiologicalSimilarity <- function(sces, objectName, filterName, valueName, defaultValue=0) {
  objIndex <- which(names(sces[[1]]) == objectName)
  sobj <- sces[[1]][[objIndex]]
  filterIndex  <- which(colnames(sobj) == filterName)
  valueIndex  <- which(colnames(sobj) == valueName)
  
  genelist<-c()
  for (i in 1:length(sces)) {
    sce<-sces[[i]]
    sobj <- sce[[objIndex]]
    sobj<-sobj[sobj[, filterIndex] < 0.01, ]
    sgene <- rownames(sobj)[order(abs(sobj[,valueName]), decreasing = TRUE)][1:min(50, dim(sobj)[1])]
    genelist<-c(genelist, sgene)
  }
  genelist<-unique(genelist)
  
  sdata <- NULL
  for (i in 1:length(sces)) {
    sce<-sces[[i]]
    sobj <- sce[[objIndex]]
    matchid <- rownames(sobj) %in% genelist
    filtered<-sobj[matchid,]
    sdata<-rbind(sdata, data.frame(Sample=names(sces)[i], Feature=rownames(filtered), Value=filtered[,valueIndex]))
  }
  
  mdata<-dcast(sdata, Feature~Sample, value.var="Value", fill=defaultValue)
  rownames(mdata)<-mdata$Feature
  mdata<-as.matrix(mdata[,c(2:ncol(mdata))])
  
  heatmap.2(mdata, margins = c(5, 10), cexRow = 0.5)
}

plotPathwaySimilarity <- function(sces, objectName, filterName, organism) {
  objIndex <- which(names(sces[[1]]) == objectName)
  sobj <- sces[[1]][[objIndex]]
  filterIndex  <- which(colnames(sobj) == filterName)

  sdata<-NULL
  for (i in 1:length(sces)) {
    sce<-sces[[i]]
    sobj <- sce[[objIndex]]
    sobj<-sobj[sobj[, filterIndex] < 0.01, ]
    sgenes <- rownames(sobj)
    spathway<-WebGestaltR(enrichMethod="ORA",organism=organism,
                    enrichDatabase="pathway_KEGG",interestGene=sgenes,
                    interestGeneType="genesymbol",referenceSet="genome",
                    is.output=FALSE)
    sdata<-rbind(sdata, data.frame(Sample=names(sces)[i],
                               Pathway=spathway$description,
                               FDR=-log10(spathway$FDR)))
  }
  
  mdata<-dcast(sdata, Feature~Sample, value.var="Value", fill=defaultValue)
  rownames(mdata)<-mdata$Feature
  mdata<-as.matrix(mdata[,c(2:ncol(mdata))])
  
  heatmap.2(mdata, margins = c(5, 10), cexRow = 0.5)
}
