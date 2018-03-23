scRNABatchQC<-function(sampleTable, organism, outputFile, cache){
  if(!missing(cache) & cache){
    plotData<-prepareReportData(sampleTable, organism, outputFile)
  }else{
    plotData<-prepareReportData(sampleTable, organism)
  }
  
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
