
#' plotDensity 
#' @description plot the distribution of total_counts, total_features, pct_counts_Mt (percentage of mtRNA counts) or pct_counts_rRNA (percentage of rRNA counts) for multiple single cell RNAseq datasets
#' @param sces a list of SingleCellExperiment objects; each object containing QC metadata for each dataset; (results from \code{\link{Process_scRNAseq}} )
#' @param feature string; which features to plot; features can be total_counts, total_features, pct_counts_Mt, pct_counts_rRNA; (default: total_counts)
#' @param featureLabel string; label of feature
#' @param scolors vector of color
#' @param lineSize integer; size of line
#' @import ggplot2 SingleCellExperiment
#' @export
#' @examples 
#' library(scRNABatchQC)
#' sces<-Process_scRNAseq(inputfiles=c("https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar1.csv.gz",
#'                                     "https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar5.csv.gz"))
#'                 
#' plotDensity(sces,"total_counts")
#' # plot the total counts distribution of the first dataset
#' plotDensity(sces[1],"total_counts")
#' @seealso \code{\link{Process_scRNAseq}} 
plotDensity <- function(sces, feature=c("total_counts","total_features","pct_counts_Mt","pct_counts_rRNA"), 
                        featureLabel=NULL, scolors = 1:length(sces), lineSize = 1 ) {
  
  feature<-match.arg(feature)
  featureData <- .getRawColData(sces, feature)
  featureLabel<- ifelse(is.null(featureLabel),feature, featureLabel)
  
  p <- ggplot(featureData, aes(x = Value)) + 
    stat_density(aes(color = Sample), size = lineSize,geom="line",position="identity") + 
    scale_colour_manual(values = scolors) +
    xlab(featureLabel) + ylab("Density")+ theme_classic()+ 
    guides(col = guide_legend(ncol=ceiling(length(sces)/10)))+
    ggtitle(feature)+
    theme(plot.title = element_text(hjust = 0.5))

  return(p)
}


#' plot the gene count distribution
#' @description plot the gene count distribution of top n (default:500) highly expressed genes 
#' @param sces list; a list of SingleCellExperiment objects; each object containing QC metadata for each dataset; (results from \code{\link{Process_scRNAseq}} )
#' @param ngenes integer; the number of highly expressed genes to plot the distribution (default:500)
#' @param scolors a vector of integer; the color of each dataset (default: 1:length(sces))
#' @param lineSize integer; the line size of the plot  (default: 1)
#' @import ggplot2 SingleCellExperiment
#' @export
#' @examples 
#' library(scRNABatchQC)
#' sces<-Process_scRNAseq(inputfiles=c("https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar1.csv.gz",
#'                                     "https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar5.csv.gz"))
#'   
#' #plot the average count distribution for all datasets
#' plotGeneCountDistribution(sces)
#' #plot the average count distribution for the first dataset
#' plotGeneCountDistribution(sces[1])
#' @seealso \code{\link{Process_scRNAseq}} 
plotGeneCountDistribution <- function(sces, scolors = 1:length(sces), 
                                      ngenes = 500, lineSize = 1) {
  prop_mat <- c()
  if (is.null(names(sces))) names(sces)<-1:length(sces)
  for (i in 1:length(sces)) {
    aveCountSum <- sum(sces[[i]]@rowRanges@elementMetadata$ave.counts)
    aveCountProp <- cumsum(sort(sces[[i]]@rowRanges@elementMetadata$ave.counts, decreasing = TRUE)[1:ngenes]) / aveCountSum
    prop_mat <- cbind(prop_mat, aveCountProp)
    colnames(prop_mat)[i] <- names(sces)[i]
  }
  rownames(prop_mat) <- seq_len(nrow(prop_mat))
  
  prop_to_plot <- reshape2::melt(prop_mat[,,drop=FALSE ])
  colnames(prop_to_plot) <- c("Feature", "Sample", "Proportion_Library")
  
  p <- ggplot(prop_to_plot, 
              aes_string(x = "Feature", y = "Proportion_Library", 
                         group = "Sample", colour = "Sample")) +
    geom_line(size = lineSize) + 
    xlab("Number of genes") + ylab("Cumulative proportion of library") +
    scale_color_manual(values = scolors) +
    theme_classic()+guides(col = guide_legend(ncol=ceiling(length(sces)/10)))+
    ggtitle("Cumulative percentage of highly expressed genes")+
    theme(plot.title = element_text(hjust = 0.5))
  
  return(p)
}

####averge count vs. detection rate
#' plotAveCountVSdetectRate 
#' @description plot the gene detection rates vs. the average gene expression level in the datasets  
#' @param sces list; a list of SingleCellExperiment objects; each object containing QC metadata for each dataset; (results from \code{\link{Process_scRNAseq}} )
#' @param scolors a vector of integer; the color of each dataset (default: 1:length(sces))
#' @param lineSize integer; the line size of the plot  (default: 1)
#' @import ggplot2 SingleCellExperiment
#' @export
#'
#' @examples 
#' library(scRNABatchQC)
#' sces<-Process_scRNAseq(inputfiles=c("https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar1.csv.gz",
#'                                     "https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar5.csv.gz"))
#'   
#' #plot the gene detection rates for all datasets
#' plotAveCountVSdetectRate(sces)
#' #plot the gene detection rates for the first dataset
#' plotAveCountVSdetectRate(sces[1])
#' @seealso \code{\link{Process_scRNAseq}} 
plotAveCountVSdetectRate <- function(sces, scolors = 1:length(sces), lineSize = 1) {
  avedetect <- data.frame()
  if (is.null(names(sces))) names(sces)<-1:length(sces)
  for (i in 1:length(sces)) {
    tmpavedec <- data.frame(avecount = log10(sces[[i]]@rowRanges@elementMetadata$ave.counts), 
                            detectrate = sces[[i]]@rowRanges@elementMetadata$num.cells / sces[[i]]@metadata$rawmeta$ncell, 
                            Sample = rep(names(sces)[i], length(log10(sces[[i]]@rowRanges@elementMetadata$ave.counts))))
    avedetect <- rbind(avedetect, tmpavedec)
  }
  
  p <- ggplot(avedetect, aes_string(x = "avecount", y = "detectrate", 
                                    group = "Sample", colour = "Sample")) + 
    geom_smooth(size=lineSize) + xlab("Average count of genes (log10)") + ylab("Detection rate") +
    scale_color_manual(values = scolors) +
    theme_classic()+guides(col = guide_legend(ncol=ceiling(length(sces)/10)))+
    ggtitle("Detection rate")+
    theme(plot.title = element_text(hjust = 0.5))
  
  return(p)
}

#' plot Mean~variance trend
#' @description plot the gene mean vs. variance trend  
#' @param sces list; a list of SingleCellExperiment objects; each object containing QC metadata for each dataset; (results from \code{\link{Process_scRNAseq}} )
#' @param scolors a vector of integer; the color of each dataset (default: 1:length(sces))
#' @param pointSize; the point size of each dot representing each gene
#' @param lineSize integer; the line size of the fitted mean-variance trend curve  (default: 1)
#' @import ggplot2 SingleCellExperiment
#' @export
#'
#' @examples 
#' library(scRNABatchQC)
#' sces<-Process_scRNAseq(inputfiles=c("https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar1.csv.gz",
#'                                     "https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar5.csv.gz"))
#'   
#' #plot the mean~variance trend for all datasets
#' plotVarianceTrend(sces)
#' #plot the mean~variance trend for the first dataset
#' plotVarianceTrend(sces[1])
#' @seealso \code{\link{Process_scRNAseq}} 
plotVarianceTrend <- function(sces, scolors = 1:length(sces), 
                              pointSize=0.8, lineSize=1) {
  meanvar_dat <- data.frame()
  if (is.null(names(sces))) names(sces)<-1:length(sces)
  for (i in 1:length(sces)) {
    tmpmeanvar <- data.frame(mean = sces[[i]]@rowRanges@metadata$hvg$mean, 
                             var= sces[[i]]@rowRanges@metadata$hvg$var, 
                             Sample = rep(names(sces)[i], length(sces[[i]]@rowRanges@metadata$hvg$mean)))

    meanvar_dat<- rbind(meanvar_dat, tmpmeanvar)
    
  }
  
  p <- ggplot(meanvar_dat, aes_string(x = "mean", y = "var", group = "Sample", colour = "Sample"))  + 
    geom_point(size=pointSize) + 
    geom_smooth(size=lineSize) + 
    scale_color_manual(values = scolors) + 
    xlab("Mean log-expression") + 
    ylab("Variance of log-expression") +
    theme_classic()+guides(col = guide_legend(ncol=ceiling(length(sces)/10)))+
    ggtitle("Mean-Variance trend")+
    theme(plot.title = element_text(hjust = 0.5))
  return(p)
}

#' plot variance explained by one specific feature
#' @description plot the variances explained by the explanatory variable 
#' @param sces list; a list of SingleCellExperiment objects; each object containing QC metadata for each dataset; (results from \code{\link{Process_scRNAseq}} )
#' @param feature string; the explanatory variable to plot; feature should be genevar_by_counts, genevar_by_features, genevar_by_Mt, or genevar_by_rRNA; (default: genevar_by_counts)
#' @param scolors a vector of color; the color of each dataset (default: 1:length(sces))
#' @param lineSize integer; the line size of the fitted mean-variance trend curve  (default: 1)
#' @import ggplot2 SingleCellExperiment
#' @export
#'
#' @examples 
#' library(scRNABatchQC)
#' sces<-Process_scRNAseq(inputfiles=c("https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar1.csv.gz",
#'                                     "https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar5.csv.gz"))
#'   
#' #plot the variance explained by the total_counts 
#' plotVarianceExplained(sces)
#' #plot the variance explained by the total number of genes
#' plotVarianceExplained(sces,feature="genevar_by_features")
#' @seealso \code{\link{Process_scRNAseq}} 
plotVarianceExplained <- function(sces, feature=c("genevar_by_counts","genevar_by_features","genevar_by_Mt","genevar_by_rRNA"), scolors = 1:length(sces), 
                                  lineSize = 1) {
  if (is.null(names(sces))) names(sces)<-1:length(sces)
  feature<-match.arg(feature)
  featureData <- .getGeneData(sces, feature)
  featureLabel<- switch(feature,"genevar_by_counts"="by total counts","genevar_by_features"="by total number of genes","genevar_by_Mt"="by total number of mtRNA", "genevar_by_rRNA"="by total number of rRNA")
  p <- ggplot(featureData, aes(x = Value)) + 
    stat_density(aes(color = Sample), size = lineSize,geom="line",position="identity") + 
    geom_vline(xintercept = 1, linetype = 2) + 
    scale_x_log10(breaks = 10 ^ (-3:2), labels = c(0.001, 0.01, 0.1, 1, 10, 100)) + 
    xlab(paste0("% variance explained ",featureLabel, "(log10)")) + 
    ylab("Density") + 
    scale_color_manual(values = scolors) +
    coord_cartesian(xlim = c(10 ^ (-3), 100)) + theme_classic()+guides(col = guide_legend(ncol=ceiling(length(sces)/10)))+
    ggtitle(paste0("Var explained ",featureLabel))+
    theme(plot.title = element_text(hjust = 0.5))
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

#' plot the similarity across sample sets
#' @description plot the average expression similarity across single cell RNAseq datasets
#' @param sces list; a list of SingleCellExperiment objects; each object containing QC metadata for each dataset; (results from \code{\link{Process_scRNAseq}} )
#' @param ... parameters passing to function pairs
#' @import SingleCellExperiment 
#' @export
#'
#' @examples 
#' library(scRNABatchQC)
#' sces<-Process_scRNAseq(inputfiles=c("https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar1.csv.gz",
#'                                     "https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar5.csv.gz"))
#'   
#' plotSampleSimilarity(sces)

#' @seealso \code{\link{Process_scRNAseq}} 

plotSampleSimilarity <- function(sces, ...) {
  
  if (length(sces)<=1) stop("there should be more than one dataset")
  if (is.null(names(sces))) names(sces)<-1:length(sces)
  
  aveCount <- sces[[1]]@rowRanges@elementMetadata$ave.counts
  
  for (i in 2:length(sces)) {
    aveCount <- merge(aveCount, sces[[i]]@rowRanges@elementMetadata$ave.counts, by = "row.names")
    rownames(aveCount) <- aveCount[, 1]
    aveCount <- aveCount[, -1]
  }
  colnames(aveCount) <- names(sces)
  
  pairs(log10(aveCount), xaxt = "n", yaxt = "n", upper.panel = panel.cor, gap = 0, lower.panel = panel.dot, ...)
}

##plot the highly variable genes
#' plot HVGs
#' @description plot the highly variable genes across single cell RNAseq datasets
#' @param sces list; a list of SingleCellExperiment objects; each object containing QC metadata for each dataset; (results from \code{\link{Process_scRNAseq}} )
#' @param margins margins for heatmap.2 (default: c(5,5))
#' @param keysize integer for heatmap.2 (default:1)
#' @param col color for heatmap.2    (default: bluered(75))
#' @param ... parameters passing to heatmap.2
#' @return a matrix containing the z-score of dispersion of highly variable genes in each dataset
#' @import SingleCellExperiment 
#' @export
#'
#' @examples 
#' library(scRNABatchQC)
#' sces<-Process_scRNAseq(inputfiles=c("https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar1.csv.gz",
#'                                     "https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar5.csv.gz"))
#' plotHVGs(sces)
#' @seealso \code{\link{Process_scRNAseq}} 

plotHVGs<-function(sces,margins=c(5,5),keysize=1,col=bluered(75),...){
  
  hvgBiologicalSimilarity <- .getBiologicalSimilarity(sces, objectName = "hvg", valueName = "zval")
  if (nrow(hvgBiologicalSimilarity)<2 | ncol(hvgBiologicalSimilarity)<2 ) {
      warning("there are less than two HVGs or less than two datasets, no heatmap will be generated")
    }	else{
      heatmap.2(hvgBiologicalSimilarity, cexRow = 0.5, margins = margins, keysize=keysize, col=col, main="HVGs",key.title="zval",key.xlab="",key.ylab="",...)
    }
   return(hvgBiologicalSimilarity)
}

##plot the pathways enriched in highly variable genes
#' plot pathways enriched in HVGs
#' @description plot the pathways enriched in the highly variable genes across single cell RNAseq datasets
#' @param sces list; a list of SingleCellExperiment objects; each object containing QC metadata for each dataset; (results from \code{\link{Process_scRNAseq}} )
#' @param margins margins for heatmap.2 (default: c(5,10))
#' @param keysize integer for heatmap.2  (default: 1)
#' @param col color for heatmap.2  (default: colorpanel(75,low="white",high="red"))
#' @param ... parameters passing to heatmap.2
#' @return a matrix containing the -log10 FDR of enriched pathways in each dataset
#' @import SingleCellExperiment 
#' @export
#'
#' @examples 
#' library(scRNABatchQC)
#' sces<-Process_scRNAseq(inputfiles=c("https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar1.csv.gz",
#'                                     "https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar5.csv.gz"))
#' plotHVGsPathways(sces)
#' @seealso \code{\link{Process_scRNAseq}} 

plotHVGsPathways<-function(sces,margins=c(5,10),keysize=1,col=colorpanel(75,low="white",high="red"),...){
  
  hvgPathways <- .getMultiplePathways(sces, metaObjectName = "hvgPathway")
  if (nrow(hvgPathways)<2 | ncol(hvgPathways)<2 ) {
    warning("there are less than two enriched pathways or less than two datasets, no heatmap will be generated")
    }	else{
    heatmap.2(hvgPathways, cexRow = 0.5, margins = margins, keysize=keysize, col=col, main="Pathways enriched in HVGs",key.title="-logFDR",key.xlab="",key.ylab="", ...)
  }
  return(hvgPathways)
}

##plot  genes highly associated with a certain principle component
#' plot PC-related genes
#' @description plot the genes highly associated with a certain principle component
#' @param sces list; a list of SingleCellExperiment objects; each object containing QC metadata for each dataset;(results from \code{\link{Process_scRNAseq}} )
#' @param margins margins for heatmap.2 (default: c(5,5))
#' @param keysize integer for heatmap.2   (default: 1)
#' @param col color for heatmap.2    (default: bluered(75))
#' @param ... parameters passing to heatmap.2
#' @return a matrix containing the log2 fold change of genes highly associated with a certain principle component
#' @import SingleCellExperiment 
#' @export
#'
#' @examples 
#' library(scRNABatchQC)
#' sces<-Process_scRNAseq(inputfiles=c("https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar1.csv.gz","https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar5.csv.gz"))
#' plotPCgenes(sces)
#' @seealso \code{\link{Process_scRNAseq}}  

plotPCgenes<-function(sces,margins=c(5,5),keysize=1,col=bluered(75),...){

   pc1geneBiologicalSimilarity <- .getBiologicalSimilarity(sces, objectName = "pc1genes", valueName = "logFC")
  if (nrow(pc1geneBiologicalSimilarity)<2 | ncol(pc1geneBiologicalSimilarity)<2)  {
     warning("there are less than two genes or less than two datasets, no heatmap will be generated") 
	 } else{
    heatmap.2(pc1geneBiologicalSimilarity, cexRow = 0.5, margins = margins, keysize=keysize, col=col, main="PCgenes",key.title="logFC",key.xlab="",key.ylab="", ...)
  }
  return(pc1geneBiologicalSimilarity)
}


##plot the pathways enriched in PC genes
#' plot pathways enriched in PC-related genes 
#' @description plot the pathways enriched in the PC genes across single cell RNAseq datasets
#' @param sces list; a list of SingleCellExperiment objects; each object containing QC metadata for each dataset; (results from \code{\link{Process_scRNAseq}} )
#' @param margins margins for heatmap.2 (default: c(5,10))
#' @param keysize integer for heatmap.2  (default: 1)
#' @param col color for heatmap.2  (default: colorpanel(75,low="white",high="red"))
#' @param ... parameters passing to heatmap.2
#' @return a matrix containing the -log10 FDR of enriched pathways in each dataset
#' @import SingleCellExperiment 
#' @export
#'
#' @examples 
#' library(scRNABatchQC)
#' sces<-Process_scRNAseq(inputfiles=c("https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar1.csv.gz",
#'                                     "https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar5.csv.gz"))
#' plotPCPathways(sces)
#' @seealso \code{\link{Process_scRNAseq}} 

plotPCPathways<-function(sces,margins=c(5,10),keysize=1,col=colorpanel(75,low="white",high="red"),...){
  
  pc1Pathways <- .getMultiplePathways(sces, metaObjectName = "pc1Pathway")
  
  if (nrow(pc1Pathways)<2 | ncol(pc1Pathways)<2 ) {
    warning("there are less than two enriched pathways or less than two datasets, no heatmap will be generated")
    }	else{
    
    heatmap.2(pc1Pathways, cexRow = 0.5, margins = margins, keysize=keysize, col=col, main="Pathways enriched in PCgenes",key.title="-logFDR",key.xlab="",key.ylab="", ...)
  }
  return(pc1Pathways)
}

##plot functions on the combined datasets
####################### PCA ##############

#' plot PCA for the combined dataset
#' @description  PCA plot for the combined dataset
#' @param scesMerge SingleCellExperiment object; this object contains the combined data and reduced dimensions using PCA and tSNE; (results from \code{\link{Combine_scRNAseq}})
#' @param scolors a vector of color
#' @param pointSize integer; size of point
#' @import SingleCellExperiment ggplot2
#' @export
#'
#' @examples 
#' library(scRNABatchQC)
#' sces<-Process_scRNAseq(inputfile=c("https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar1.csv.gz",
#'                                    "https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar5.csv.gz"))
#' scesMerge <- Combine_scRNAseq(sces)
#' plotAllPCA(scesMerge)

plotAllPCA <- function(scesMerge, scolors = NULL, pointSize = 0.8) {
  pcadata <- data.frame(Sample = (scesMerge@colData$condition), PC1 = scesMerge@metadata$reducedDims$PCA$x[, 1], PC2 = scesMerge@metadata$reducedDims$PCA$x[, 2],stringsAsFactors=FALSE )
  nsample<-length(unique(scesMerge@colData$condition))
  if (is.null(scolors)) scolors=1:nsample
  eigs <- scesMerge@metadata$reducedDims$PCA$sdev ^ 2
  pc1pct <- eigs[1] / sum(eigs)
  pc2pct <- eigs[2] / sum(eigs)
  
  p_pca <- ggplot(pcadata, aes(x = PC1, y = PC2)) + 
    geom_point(aes(col = Sample), size = pointSize) + 
    xlab(paste0("PC1(", round(pc1pct * 100), "%)")) + 
    ylab(paste0("PC2(", round(pc2pct * 100), "%)")) + 
    scale_colour_manual(values = scolors,breaks=unique(pcadata$Sample)) +
    theme_classic()+
    guides(col = guide_legend(ncol=ceiling(nsample/10)))+
    ggtitle("PCA")+
    theme(plot.title = element_text(hjust = 0.5))
  
  return(p_pca)
}

####################### tSNE ##############
#' plot tSNE for the combined dataset
#' @description  tSNE plot for the combined dataset
#' @param scesMerge SingleCellExperiment object; this object contains the combined data and reduced dimensions using PCA and tSNE;(results from \code{\link{Combine_scRNAseq}}) 
#' @param scolors a vector of color
#' @param pointSize integer; size of point
#' @import SingleCellExperiment ggplot2
#' @export
#'
#' @examples 
#' library(scRNABatchQC)
#' sces<-Process_scRNAseq(inputfile=c("https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar1.csv.gz",
#'                                    "https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar5.csv.gz"))
#' scesMerge<-Combine_scRNAseq(sces)
#' plotAlltSNE(scesMerge)

plotAlltSNE <- function(scesMerge, scolors = NULL, pointSize = 0.8) {
  tsnedata <- data.frame(D1 = scesMerge@metadata$reducedDims$tSNE[, 1], D2 = scesMerge@metadata$reducedDims$tSNE[, 2], Sample = (scesMerge@colData$condition),stringsAsFactors=FALSE)
  nsample<-length(unique(scesMerge@colData$condition))
  if (is.null(scolors)) scolors=1:nsample
  p_tsne <- ggplot(tsnedata, aes(x = D1, y = D2, label = Sample)) + 
    geom_point(aes(col = Sample), size = pointSize) + labs(title="tSNE",x="t-SNE1",y="t-SNE2")+
    scale_colour_manual(values = scolors,breaks=unique(tsnedata$Sample)) + theme_classic()+guides(col = guide_legend(ncol=ceiling(nsample/10)))+
    theme(plot.title = element_text(hjust = 0.5))
  
  return(p_tsne)
}

####pairwise comparison between any two datasets
#' plot differentially expressed genes between any two datasets
#' @description  plot the differentially expressed genes in any pairwise comparison
#' @param scesMerge a SingleCellExperiment object; this object contains the combined datasets, pairwise comparison results and reduced dimensions using PCA and tSNE 
#' @param margins margins for heatmap.2
#' @param keysize integer for heatmap.2
#' @param col color for heatmap.2
#' @param ... parameters passing to heatmap.2
#' @return a matrix containing the log2 fold change of differnetially expressed genes in each comparison
#' @import SingleCellExperiment ggplot2
#' @export
#'
#' @examples 
#' library(scRNABatchQC)
#' sces<-Process_scRNAseq(inputfile=c("https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar1.csv.gz",
#'                                    "https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar5.csv.gz"))
#' scesMerge<-Combine_scRNAseq(sces)
#' plotDiffgenes(scesMerge)
#' @seealso \code{\link{Process_scRNAseq}}, \code{\link{Combine_scRNAseq}} 

plotDiffgenes<-function(scesMerge,margins=c(5,5),keysize=1,col=bluered(75), ...){
  if (is.null(scesMerge@metadata$diffFC$genes)) stop("no differentially expressed genes detected in pairwise comparisons")
  if (nrow(scesMerge@metadata$diffFC$genes)<2 |ncol(scesMerge@metadata$diffFC$genes)<2 ) {
   warning("less than two genes detected or less than two columns, no heatmap will be generated") 
    } else{
   heatmap.2(as.matrix(scesMerge@metadata$diffFC$genes), cexRow = 0.6, margins=margins, keysize=keysize, col=col,key.title="logFC",key.xlab="",key.ylab="", ...)
 }
  return(scesMerge@metadata$diffFC$genes)
}

#' plot pathways enriched in differentially expressed genes
#' @description  plot the differentially expressed genes in any pairwise comparison
#' @param scesMerge a SingleCellExperiment object; this object contains the combined datasets, pairwise comparison results and reduced dimensions using PCA and tSNE; (results from \code{\link{Combine_scRNAseq}}) 
#' @param margins margins for heatmap.2
#' @param keysize integer for heatmap.2
#' @param col color for heatmap.2
#' @param ... parameters passing to heatmap.2
#' @return a matrix containing the -log10 FDR of enriched pathways 
#' @import SingleCellExperiment ggplot2
#' @export
#'
#' @examples 
#' library(scRNABatchQC)
#' sces<-Process_scRNAseq(inputfile=c("https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar1.csv.gz",
#'                                    "https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar5.csv.gz"))
#' scesMerge<-Combine_scRNAseq(sces)
#' plotDiffPathways(scesMerge)
#' @seealso \code{\link{Process_scRNAseq}}, \code{\link{Combine_scRNAseq}}


plotDiffPathways <-function(scesMerge,margins=c(5,5),keysize=1,col=colorpanel(75,low="white",high="red"), ...){
  if (is.null(scesMerge@metadata$diffFC$pathways)) stop("no pathways detected") 
  if (nrow(scesMerge@metadata$diffFC$pathways)<2 | ncol(scesMerge@metadata$diffFC$pathways)<2) {
    warning("less than two pathways detected or less than two columns, no heatmap will be generated")
	} else{
    heatmap.2(as.matrix(scesMerge@metadata$diffFC$pathways), cexRow = 0.6, margins=margins, keysize=keysize, col=col,key.title="-logFDR",key.xlab="",key.ylab="", ...)
  }
  return(scesMerge@metadata$diffFC$pathways)
}
