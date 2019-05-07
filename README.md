scRNABatchQC
==========
* [Introduction](#introduction)
* [Installation](#installation)
* [Usage](#example)

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

	install_github("liuqivandy/scRNABatchQC")
  
<a name="example"/>

# Usage

After installing scRNABatchQC, use following codes to run examples

	library(scRNABatchQC)
	
	#check and compare the quality of two scRNA-seq datasets from Retinal Bipolar Neurons (Cell 2016 Aug 25;166(5):1308-1323 )
	#The scRNA-seq dataset should be provided by the gene(row)-cell(column) matrix and the rowname should be gene symbol
	# Since one dataset has more than 14,000 cells and 2,4904 genes, big memory (>=16Gb) is required.
	#a report named "report.html" will be generated in your working directory
	
	#One SingleCellExperiment object (scesMerge) containing the combined dataset and a list of SingleCellExperiment objects (sces, each object contains the preprocessed dataset and metadata for one dataset) will be returned.
	
	#It took about 2.5 min on a MacBook Pro laptop with 3.1 GHz Intel Core i5 and 16 Gb memory, and about 3 min on a Windows desktop with an Intel(R) Xeon(R) CPU E5-2620 0 at 2 GHz and 32 GB memory. 	
	
	result<-scRNABatchQC(inputfiles=c("https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar1.csv.gz", "https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar5.csv.gz"),organism="mmusculus")
	
	#the number of genes and the number of cells in the combined dataset
	dim(result$scesMerge)
	
	#the number of genes and the number of cells in the first dataset after filtering
	dim(result$sces[[1]])

	#check the quality of six seminiferous tubule (ST) datasets of murine spermatogenesis from GEO (GSE112393), each dataset has 25,000 genes and 2,600~5,500 cells
	
	#a report named "report.html" will be generated in your working directory
	
	#It took about 3.6 min on a MacBook Pro laptop with 3.1 GHz Intel Core i5 and 16 Gb memory and about 5 min on  a Windows desktop with an Intel(R) Xeon(R) CPU E5-2620 0 at 2 GHz and 32 GB memory.
       scRNABatchQC(inputfiles=c("ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3069nnn/GSM3069439/suppl/GSM3069439_ST1_DGE.txt.gz", "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3069nnn/GSM3069440/suppl/GSM3069440_ST2_DGE.txt.gz", "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3069nnn/GSM3069441/suppl/GSM3069441_ST3_DGE.txt.gz","ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3069nnn/GSM3069442/suppl/GSM3069442_ST4_DGE.txt.gz","ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3069nnn/GSM3069443/suppl/GSM3069443_ST5_DGE.txt.gz","ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3069nnn/GSM3069444/suppl/GSM3069444_ST6_DGE.txt.gz"),organism="mmusculus")

	#check the quality of three scRNA-seq datasets from embryonic development of the thymus
	# It took about 41 sec on a MacBook Pro laptop with 3.1 GHz Intel Core i5 and 16 Gb memory and about 1 min on a Windows desktop with an Intel(R) Xeon(R) CPU E5-2620 0 at 2 GHz and 32 GB memory. 	   
	scRNABatchQC(inputfiles=c("ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2883nnn/GSM2883184/suppl/GSM2883184_E12_5_wholeThy_venus_1.dge.txt.gz","ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2883nnn/GSM2883185/suppl/GSM2883185_E12_5_wholeThy_venus_2.dge.txt.gz","ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2883nnn/GSM2883186/suppl/GSM2883186_E12_5_wholeThy_venus_3.dge.txt.gz"),organism="mmusculus"

