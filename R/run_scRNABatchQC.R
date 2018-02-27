source("R/prepareSCRNAData.R")
source("R/prepareSCRNADataSet.R")
source("R/plotFunctions.R")

sampleTable <- data.frame(Sample = c("S1", "S2", "S3"),
                        File = file.path("Z:/JiePing/scRNABatchQC", c("count1.csv", "count2.csv", "count3.csv")),
                        Transform = c(0, 0, 0))
sces <- prepareSCRNADataSet(sampleTable)

sceall <- preparePCATSNEData(sces)

rmarkdown::render("R/scRNABatchQCreport.Rmd", output_dir = "Z:/JiePing/scRNABatchQC",
                  params = list(data = sces, all = sceall))
