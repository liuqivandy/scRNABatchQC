scRNABatchQC<-function(inputfiles,samplenames=NULL, organism=c("hsapiens","mmusculus"),outputFile="report.html", sampleRatio=1,cache=FALSE, nHVGs=1000,nPC=10,sf=10000,logFC=0.58,FDR=0.01){
  organism<-match.arg(organism)
  if(!missing(cache) & cache){
    plotData<-prepareReportData(inputfiles,samplenames, organism, outputFile, sampleRatio,nHVGs, nPC, sf, logFC, FDR)
  }else{
    plotData<-prepareReportData(inputfiles,samplenames, organism, sampleRatio, nHVGs, nPC, sf, logFC, FDR)
  }
  
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
