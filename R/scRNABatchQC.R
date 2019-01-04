scRNABatchQC<-function(inputfiles,samplenames=NULL, organism="hsapiens", outputFile="report.html", cache=TRUE){
  if(!missing(cache) & cache){
    plotData<-prepareReportData(inputfiles,samplenames, organism, outputFile)
  }else{
    plotData<-prepareReportData(inputfiles,samplenames, organism)
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
