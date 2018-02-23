
library(RColorBrewer)
library(scater)

source("prepareSCRNADataSet.R")
source("plotFunctions.R")

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

DEBUG<-TRUE

if(DEBUG){
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
  
  plotGeneCountDistribution(sces, scolors)

  plotAveCountVSdetectRate(sces, scolors)
  
  plotVarianceTrend(sces, scolors)
}