scRNABatchQC<-function(inputfiles,samplenames=NULL, organism=c("hsapiens","mmusculus"),outputFile="report.html", sampleRatio=1,cache=FALSE){
  if(!missing(cache) & cache){
    plotData<-prepareReportData(inputfiles,samplenames, organism, outputFile, sampleRatio)
  }else{
    plotData<-prepareReportData(inputfiles,samplenames, organism, sampleRatio)
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
