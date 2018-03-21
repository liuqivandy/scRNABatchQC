library(ggplot2)
library(gplots)
library(reshape2)
library(SingleCellExperiment)
library(scater)
library(scran)
library(cluster)
library(limma)
library(dynamicTreeCut)
library(Rtsne)
library(data.table)
library(WebGestaltRsqh)
library(rmdformats)

source("./R/dataFunctions.R")
source("./R/plotFunctions.R")
source("./R/prepareSCRNAData.R")
source("./R/prepareSCRNADataSet.R")
source("./R/prepareReportData.R")

rdatafile<-"Z:/JiePing/scRNABatchQC/sces.rdata"
organism = "mmusculus"

if(!file.exists(rdatafile)){
  sampleTable <- data.frame(Sample = c("S1", "S2", "S3"),
                            File = file.path("Z:/JiePing/scRNABatchQC", c("count1.csv", "count2.csv", "count3.csv")))
  
  plotData<-prepareReportData(sampleTable, organism)
  
  save(file=rdatafile, plotData)
}else{
  load(rdatafile)
}

rmarkdown::render("./R/scRNABatchQCreport.Rmd",
                  output_dir = "Z:/JiePing/scRNABatchQC",
                  output_file = "scRNABatchQCreport.html",
                  params = list(data = plotData, displayFigure = TRUE))
