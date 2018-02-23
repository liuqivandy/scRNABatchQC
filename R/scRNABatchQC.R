
library(RColorBrewer)
library(scater)

source("prepareSCRNADataSet.R")

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

plotDensity <- function(sces, feature, featureLabel="", scolors ) {
  featureData<-.getColData(sces, feature)
  featureLabel=ifelse(featureLabel=="", feature, featureLabel)
  
  g<-ggplot(featureData, aes(x=Value)) + 
    geom_density(aes(color=Sample)) + 
    scale_colour_manual(values=scolors) +
    xlab(featureLabel) +
    theme_bw()
  print(g)
}

plotCumulativeProportion<-function (sces, scolors, nfeatures = 500, exprs_values = "counts") 
{
  exprs_mat = .getCbindRowData(sces, "ave.count")
  nfeatures_total <- nrow(exprs_mat)
  seq_real_estate <- t(plyr::aaply(exprs_mat, 2, .fun = function(x) {
    cumsum(sort(x, decreasing = TRUE))
  }))
  rownames(seq_real_estate) <- seq_len(nfeatures_total)
  nfeatures_to_plot <- nfeatures
  to_plot <- seq_len(nfeatures_to_plot)
  seq_real_estate_long <- reshape2::melt(seq_real_estate[to_plot, ], value.name = exprs_values)
  prop_library <- reshape2::melt(t(t(seq_real_estate[to_plot, ])/scater:::.general_colSums(exprs_mat)), value.name = "prop_library")
  colnames(seq_real_estate_long) <- c("Feature", "Sample", exprs_values)
  seq_real_estate_long$Sample<-as.factor(seq_real_estate_long$Sample)
  seq_real_estate_long$Proportion_Library <- prop_library$prop_library
  g <- ggplot(seq_real_estate_long) + 
    geom_line(aes(x = Feature, y = Proportion_Library, color=Sample)) + 
    scale_colour_manual(values = setNames(scolors, colnames(x)) ) +
    xlab("Number of features") + 
    ylab("Cumulative proportion of library") +
    theme_bw()
  print(g)
}

sampleTable<-data.frame(Sample=c("S1", "S2", "S3"),
                        File=paste0("Z:/shengq1/20180214_scRNABatchQC/", c("qi_m1.csv", "qi_m2.csv", "s1_pan_qi_1.csv")),
                        Transform=c(1, 1, 1))

sces<-prepareSCRNADataSet(sampleTable)

if(length(sces) < 9){
  scolors=brewer.pal(length(sces), "Set1")
}else{
  scolors=rainbow(length(sces))
}
names(scolors)<-names(sces)

plotDensity(sces,"total_counts", "Total count", scolors)
plotDensity(sces,"total_features", "Total feature", scolors)
plotDensity(sces,"pct_counts_Mt", "pct_counts_Mt", scolors)

### top 500 genes count distribution
plotCumulativeProportion(sces, scolors)

