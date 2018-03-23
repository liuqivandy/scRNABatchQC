scRNABatchQC<-function(sampleTable, organism, outputFile, cache){
  if(!missing(cache) & cache){
    plotData<-prepareReportData(sampleTable, organism, outputFile)
  }else{
    plotData<-prepareReportData(sampleTable, organism)
  }
  
  reportRmd <- system.file("report/scRNABatchQCreport.Rmd", package="scRNABatchQC")
  
  rmarkdown::render(reportRmd,
                    output_dir = getwd(),
                    output_file = outputFile,
                    params = list(data = plotData))
}