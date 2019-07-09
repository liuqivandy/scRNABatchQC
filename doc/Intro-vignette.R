## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
knitr::opts_chunk$set(fig.width=7, fig.height=7) 

## ---- eval=F, echo=T-----------------------------------------------------
#  library(scRNABatchQC)
#  output<-scRNABatchQC(inputfiles=c("ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2883nnn/GSM2883184/suppl/GSM2883184_E12_5_wholeThy_venus_1.dge.txt.gz",
#                                    "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2883nnn/GSM2883185/suppl/GSM2883185_E12_5_wholeThy_venus_2.dge.txt.gz",
#  			                            "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2883nnn/GSM2883186/suppl/GSM2883186_E12_5_wholeThy_venus_3.dge.txt.gz"),
#                       logFC=0.585,
#                       organism="mmusculus")

## ----message=FALSE,warning=FALSE-----------------------------------------
# list the organism supported by WebGestalt
library(scRNABatchQC)
listOrganism()

## ----echo=T,message=FALSE,warning=FALSE,error=FALSE,results='hide'-------
library(scRNABatchQC)
rdsFile<-"output.rds"
if(!file.exists(rdsFile)){
  output<-scRNABatchQC(inputfiles=c("ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2883nnn/GSM2883184/suppl/GSM2883184_E12_5_wholeThy_venus_1.dge.txt.gz",
                                    "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2883nnn/GSM2883185/suppl/GSM2883185_E12_5_wholeThy_venus_2.dge.txt.gz",
		  	                            "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2883nnn/GSM2883186/suppl/GSM2883186_E12_5_wholeThy_venus_3.dge.txt.gz"),
                       logFC=0.585,
                       organism="mmusculus",
                       createReport=F)
  saveRDS(output, rdsFile)
}else{
  output<-readRDS(rdsFile)
}

## ----results='hide'------------------------------------------------------
# a list of SingleCellExperiment objects, each of which contains raw and normalized count, technical and biological metadata from one scRNAseq dataset
output$sces
# the first SingleCellExperiment object for the first dataset (here bioplar1.csv.gz)
output$sces[[1]]
# the metadata of the first dataset
names(output$sces[[1]]@metadata)
#raw counts of the first dataset
counts(output$sces[[1]])[1:5,1:5]
#normalized and log-transformed data of the first dataset
logcounts(output$sces[[1]])[1:5,1:5]
#one SingleCellExperiment object,which contains normalized count, technical and biological metadata for the merged dataset (all datasets merged into one)
output$scesMerge
#metadata for the merged dataset
names(output$scesMerge@metadata)
# normalized and log-transformed data of the merged dataset. Note: the merged dataset doesn't contain the raw count
dim(logcounts(output$scesMerge))

## ------------------------------------------------------------------------
# plot the distribution of the total number of counts 
plotDensity(output$sces,"total_counts")
# plot the distribution of mtRNA fraction
plotDensity(output$sces,"pct_counts_Mt")
# plot the distribution of the total counts for the first sample
plotDensity(output$sces[1],"total_counts")
#plot the variance-expression trend
plotVarianceTrend(output$sces)
#plot the variance explained by the total counts
plotVarianceExplained(output$sces,feature="genevar_by_counts")

## ----results='hide'------------------------------------------------------
# plot the HVGs 
plotHVGs(output$sces)

## ----results='hide'------------------------------------------------------
# plot the pathways enriched in the HVGs
plotHVGsPathways(output$sces)

## ----results='hide'------------------------------------------------------
# plot the PC-related genes
plotPCgenes(output$sces)

## ----results='hide'------------------------------------------------------
# plot the pathways enriched in pc-related genes
plotPCPathways(output$sces)

## ------------------------------------------------------------------------
# plot the pairwise similarity of average expression 
plotSampleSimilarity(output$sces)

## ----results='hide'------------------------------------------------------
# PCA plots of all samples
plotAllPCA(output$scesMerge)

## ----results='hide'------------------------------------------------------
# tSNE plots of all samples
plotAlltSNE(output$scesMerge)

## ----results='hide'------------------------------------------------------
# plot the differentially expressed genes in any pairwise comparison
plotDiffgenes(output$scesMerge)

## ----results='hide'------------------------------------------------------
# plot the pathways enriched in DEGs
plotDiffPathways(output$scesMerge)

## ---- eval=F, echo=T-----------------------------------------------------
#  # change the chunk size to save memory
#  result<-scRNABatchQC(inputfiles=c("https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar1.csv.gz",
#  	                          "https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar5.csv.gz"), chunk=2000)
#  gc()

## ---- eval=F, echo=T-----------------------------------------------------
#  # change the sampleRatio to save memory
#  result<-scRNABatchQC(inputfiles=c("https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar1.csv.gz",
#  	                          "https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar5.csv.gz"), sampleRatio = 0.1)
#  gc()

