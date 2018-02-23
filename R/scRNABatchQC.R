
library(RColorBrewer)
library(scater)

source("prepareSCRNAData.R")
source("prepareSCRNADataSet.R")
source("plotFunctions.R")

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
  
  
  pdf("Z:/shengq1/20180214_scRNABatchQC/scRNABatchQC.pdf", onefile=TRUE)
  
  checkClusterSeparateness(sces[[1]]$sce)
  
  print(plotDensity(sces,"total_counts", "Total count", scolors))
  print(plotDensity(sces,"total_features", "Total feature", scolors))
  print(plotDensity(sces,"pct_counts_Mt", "pct_counts_Mt", scolors))
  print(plotGeneCountDistribution(sces, scolors))
  print(plotAveCountVSdetectRate(sces, scolors))
  print(plotVarianceTrend(sces, scolors))
  print(plotMultiSamplesOneExplanatoryVariables(sces, scolors, "log10_total_counts", "log10(total counts)"))
  dev.off()
}