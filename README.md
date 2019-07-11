scRNABatchQC
==========
* [Introduction](#introduction)
* [Installation](#installation)
* [Usage](#example)
* [Speed test](#speed)
<a name="introduction"/>

# Introduction

scRNABatchQC is an R package for generating a single html report to check and compare quality of multiple single cell RNA-seq datasets.

<a name="installation"/>

# Installation

Step 1:  Install pandoc (https://github.com/jgm/pandoc/releases/tag/2.2.1) (pandoc is required to convert files from markup format into html format). After installation , update your path to include the directory where pandoc’s binaries are installed. 

Step 2: Install a modified version of WebGestatR by:

	library(devtools)
	install_github("shengqh/WebGestaltR")

Step 3: Install scRNABatchQC by:

	install_github("liuqivandy/scRNABatchQC", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
  
<a name="example"/>

# Usage

After installing scRNABatchQC, use following codes to run examples:

## Example 1:

Check and compare the quality of two scRNA-seq datasets from Retinal Bipolar Neurons (Cell 2016 Aug 25;166(5):1308-1323 ). The scRNA-seq dataset should be provided by the gene(row)-cell(column) matrix and the rowname should be gene symbol. Since one dataset has more than 14,000 cells and 2,4904 genes,  **memory >= 4Gb and 64 bit system** are required. A report named "report.html" will be generated in your working directory.
	
One SingleCellExperiment object (scesMerge) containing the combined dataset and a list of SingleCellExperiment objects (sces, each object contains the preprocessed dataset and metadata for one dataset) will be returned.

```
library(scRNABatchQC)
	
start_time <- Sys.time()
result<-scRNABatchQC(inputfiles=c("https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar1.csv.gz", 
	                          "https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar5.csv.gz"),
                     organism="mmusculus")
end_time <- Sys.time()
end_time - start_time

#the number of genes and the number of cells in the combined dataset
dim(result$scesMerge)
	
#the number of genes and the number of cells in the first dataset after filtering
dim(result$sces[[1]])

```

## Example 2:

Check the quality of six seminiferous tubule (ST) datasets of murine spermatogenesis from GEO (GSE112393), each dataset has 25,000 genes and 2,600~5,500 cells
	
A report named "report.html" will be generated in your working directory

```
library(scRNABatchQC)
	
start_time <- Sys.time()
scRNABatchQC(inputfiles=c("ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3069nnn/GSM3069439/suppl/GSM3069439_ST1_DGE.txt.gz", 
                          "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3069nnn/GSM3069440/suppl/GSM3069440_ST2_DGE.txt.gz", 
                          "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3069nnn/GSM3069441/suppl/GSM3069441_ST3_DGE.txt.gz",
                          "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3069nnn/GSM3069442/suppl/GSM3069442_ST4_DGE.txt.gz",
                          "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3069nnn/GSM3069443/suppl/GSM3069443_ST5_DGE.txt.gz",
                          "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3069nnn/GSM3069444/suppl/GSM3069444_ST6_DGE.txt.gz"),
             organism="mmusculus")
end_time <- Sys.time()
end_time - start_time
```

## Example 3:
	
Check the quality of three scRNA-seq datasets from embryonic development of the thymus. 
	
```
library(scRNABatchQC)
	
start_time <- Sys.time()
scRNABatchQC(inputfiles=c("ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2883nnn/GSM2883184/suppl/GSM2883184_E12_5_wholeThy_venus_1.dge.txt.gz",
                          "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2883nnn/GSM2883185/suppl/GSM2883185_E12_5_wholeThy_venus_2.dge.txt.gz",
			  "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2883nnn/GSM2883186/suppl/GSM2883186_E12_5_wholeThy_venus_3.dge.txt.gz"),
	     organism="mmusculus")
end_time <- Sys.time()
end_time - start_time
```

<a name="speed"/>

# Speed test

|Computer|CPU|Memory|System|R version|Example1|Example2|Example3|
|---|---|---|---|---|---|---|---|
|MacBook Pro Laptop|3.1 GHz Intel Core i5|16 Gb|MacOS|3.6.0|2.5 min|3.6 min|41 sec|
|Windows Desktop|2 GHz Intel(R) Xeon(R) E5-2620|32 Gb|Windows 7|3.4.3|3 min|5 min| 1 min|
|Lenovo Laptop|1.8 GHz Intel(R) i7 8565U|16 Gb|Windows 10|3.6.0|3.4 min|5.1 min| 1.1 min|
|MacBook Pro Laptop|2.7 GHz Intel Core i5|8 Gb|MacOS|3.4.1|6.5 min|12.5 min|1.1 min|
|MacBook Pro Laptop|1.7 GHz Intel Core i5|4 Gb|MacOS|3.4.3|39 min |13 min|3 min|
|Windows Desktop|2.6 GHz Intel(R) Xeon(R) E5-2640|64 Gb|Windows 10|3.6.0|4.1 min|6 min| 1.3 min|
|ACCRE Cluster Gateway|Intel(R) Xeon(R) Gold 6154 CPU @ 3.00GHz|1000 Gb|CentOs 7|3.5.1|3.1 min|5.2 min| 28.6 sec|

