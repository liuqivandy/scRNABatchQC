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
library(WebGestaltR)

source("./R/dataFunctions.R")
source("./R/plotFunctions.R")
source("./R/prepareSCRNAData.R")
source("./R/prepareSCRNADataSet.R")

rdatafile<-"Z:/JiePing/scRNABatchQC/sces.rdata"
organism = "mmusculus"

if(!file.exists(rdatafile)){
  sampleTable <- data.frame(Sample = c("S1", "S2", "S3"),
                            File = file.path("Z:/JiePing/scRNABatchQC", c("count1.csv", "count2.csv", "count3.csv")))
  sces <- prepareSCRNADataSet(sampleTable, organism)
  scesall <- preparePCATSNEData(sces)
  save(file=rdatafile, sces, scesall)
}else{
  load(rdatafile)
}

rmarkdown::render("./R/scRNABatchQCreport.Rmd", 
                  output_dir = "Z:/JiePing/scRNABatchQC",
                  params = list(data = sces, all = scesall, organ = organism))
