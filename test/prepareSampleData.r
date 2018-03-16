
rootDir<-"Z:/shengq1/20180214_scRNABatchQC/"
sampleTable<-data.frame(Sample=c("S1", "S2", "S3"),
                          File=c("qi_m1.csv", "qi_m2.csv", "s1_pan_qi_1.csv"))

genenames<-read.csv(paste0(rootDir, "mousenames.csv"), header=F)
for(sind in c(1:nrow(sampleTable))){
  sampleFile = paste0(rootDir, sampleTable$File[sind])
  sampleName = paste0(rootDir, sampleTable$Sample[sind], ".csv")
  
  dat<-read.csv(sampleFile, header=F)
  dat<-t(dat)
  rownames(dat)<-genenames$V1
  colnames(dat)<-paste0("Cell", c(1:ncol(dat)))
  write.csv(dat, file=sampleName)
}
