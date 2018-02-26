source("prepareSCRNAData.R")
source("prepareSCRNADataSet.R")
source("plotFunctions.R")

path <- "Z:/JiePing/scRNABatchQC"
counts <- list()
counts[[1]] <- as.matrix(t(read.csv(file.path(path, "qi_m2.csv"))))
counts[[2]] <- as.matrix(t(read.csv(file.path(path, "s1_pan_qi_1.csv"))))
counts[[3]] <- as.matrix(t(read.csv(file.path(path, "qi_m1.csv"))))
genename <- read.csv(file.path(path, "mousenames.csv"), as.is = T, header = F)[, 1]
rownames(counts[[1]]) <- rownames(counts[[2]]) <- rownames(counts[[3]]) <- genename

sces <- list()

for (n in 1 : length(counts)) {
  sces[[n]] <- prepareSCRNAData(counts[[n]])
}
names(sces) <- c("S1", "S2", "S3")

sceall <- preparePCATSNEData(sces)

rmarkdown::render("R/scRNABatchQCreport.Rmd", params = list(
  sces = sces,
  sceall = sceall
))
