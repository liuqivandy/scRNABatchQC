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
library(knitr)

source("./R/dataFunctions.R")
source("./R/plotFunctions.R")
source("./R/prepareSCRNAData.R")
source("./R/prepareSCRNADataSet.R")
source("./R/prepareReportData.R")

organism = "hsapiens"

sampleNames<-c("A0_1cell", "A0_CQS", "A1_1cell", "A1_CQS", "A2_1cell", "A2_CQS", "A5_1cell", "A5_CQS")
sampleTable <- data.frame(Sample = sampleNames,
                          File = file.path("Z:/JiePing/scRNASeq/1cell_compare/", paste0(sampleNames, ".filtered.csv")))
plotData<-prepareReportData(sampleTable, organism, "qc_a")
rmarkdown::render("./R/scRNABatchQCreport.Rmd",
                  output_dir = "Z:/JiePing/scRNASeq/1cell_compare",
                  output_file = "qc_A.html",
                  params = list(data = plotData))

sampleNames<-c("H0_1cell", "H0_CQS", "H1_1cell", "H1_CQS", "H2_1cell", "H2_CQS", "H5_1cell", "H5_CQS")
sampleTable <- data.frame(Sample = sampleNames,
                          File = file.path("Z:/JiePing/scRNASeq/1cell_compare/", paste0(sampleNames, ".filtered.csv")))
plotData<-prepareReportData(sampleTable, organism, "qc_h")
rmarkdown::render("./R/scRNABatchQCreport.Rmd",
                  output_dir = "Z:/JiePing/scRNASeq/1cell_compare",
                  output_file = "qc_H.html",
                  params = list(data = plotData))
