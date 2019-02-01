#' scRNABatchQC
#' 
#' @description Compare multiple scRNA-seq datasets simultaneously on numerous technical and biological features 

#' @param inputfiles; a character vector giving the input file names (or a URL starting http://, file://, etc.) of gene-by-cell count matrices; each file should be regular delimited file;  Compressed files ending .gz and .bz2 are supported
#' @param samplenames; the names of experiments, if it is NULL(default), the names are set to S1, S2....
#' @param organism; a character string giving the organism, either hsapiens or mmusculus; (default: "hsapiens")
#' @param outputFile; a character string giving the output File name (default: report.html) 
#' @param sampleRatio float; the ratio of cells sampled from each experiment to examine the expression similarity(default: 1)
#' @param nHVGs integer; the number of highly variable genes (default: 1000)
#' @param nPC integer; the number of principle components (default: 10)
#' @param sf integer; Scale factor  (default: 10000)
#' @param logFC float; log2 fold change cutoff to select differentially expressed genes  (default: 1)
#' @param FDR float; FDR cutoff to select differentially expressed genes (default:0.01)
#' 

#' @examples
#' library(scRNABatchQC)  


#' #Check the quality of two single cell RNA-seq datasets from human
#' scRNABatchQC(inputfiles=c("data1.csv","data2.csv"))  
#' 
#' @import R.utils ggplot2 gplots limma data.table irlba Rtsne WebGestaltR rmdformats Matrix statmod DT RCurl
#' @importFrom devtools session_info
#' 
#' @export
scRNABatchQC<-function(inputfiles,samplenames=NULL, organism=c("hsapiens","mmusculus"),outputFile="report.html", sampleRatio=1,cache=FALSE, nHVGs=1000,nPC=10,sf=10000,logFC=1,FDR=0.01){
  organism<-match.arg(organism)
  if(cache){
    cacheFile<-paste0(outputFile, ".rdata")
    if(file.exists(cacheFile)){
      load(cacheFile)
    }else{
      plotData<-prepareReportData(inputfiles=inputfiles,samplenames=samplenames, organism=organism,  sampleRatio=sampleRatio,nHVGs=nHVGs, nPC=nPC, sf=sf, logFC=logFC, FDR=FDR)
      save(plotData, file=cacheFile)
      cat("Cachefile saved")
    }
  }else{
    plotData<-prepareReportData(inputfiles=inputfiles,samplenames=samplenames, organism=organism,  sampleRatio=sampleRatio,nHVGs=nHVGs, nPC=nPC, sf=sf, logFC=logFC, FDR=FDR)
  }
  plotData$sampleRatio=sampleRatio
  plotData$nHVGs = nHVGs
  plotData$nPC=nPC
  plotData$sf=sf
  plotData$logFC=logFC
  plotData$FDR=FDR

  #reportRmd <- "E:/sqh/programs/scRNABatchQC/inst/report/scRNABatchQCreport.Rmd"
  reportRmd <- system.file("report/scRNABatchQCreport.Rmd", package="scRNABatchQC")
  
  outputFile <- getAbsolutePath(outputFile)
  output_dir = dirname(outputFile)
  output_file = basename(outputFile)
  
  cat("Output report to:", outputFile, "\n")
  rmarkdown::render(reportRmd,
                    output_dir = output_dir,
                    output_file = output_file,
                    params = list(data = plotData))
}
