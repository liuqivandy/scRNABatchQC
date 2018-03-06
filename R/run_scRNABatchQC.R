source("./dataFunctions.R")
source("./plotFunctions.R")
source("./prepareSCRNAData.R")
source("./prepareSCRNADataSet.R")

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

rmarkdown::render("./scRNABatchQCreport.Rmd", 
                  output_dir = "Z:/JiePing/scRNABatchQC",
                  params = list(data = sces, all = scesall))
