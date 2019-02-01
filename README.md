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

An example to run scRNABatchQC.
        
	
	
	
	library(scRNABatchQC)
	#This will check and compare the quality of two scRNA-seq datasets from Retinal Bipolar Neurons (Cell 2016 Aug 25;166(5):1308-1323 ) and generate a report.html in your current working directory
	scRNABatchQC(inputfiles=c("https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar1.csv.gz", "https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar5.csv.gz"),organism="mmusculus")
      #check the quality of six seminiferous tubule (ST) datasets of murine spermatogenesis from GEO (GSE112393)
      scRNABatchQC(inputfiles=c("ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3069nnn/GSM3069439/suppl/GSM3069439_ST1_DGE.txt.gz", "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3069nnn/GSM3069440/suppl/GSM3069440_ST2_DGE.txt.gz", "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3069nnn/GSM3069441/suppl/GSM3069441_ST3_DGE.txt.gz","ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3069nnn/GSM3069442/suppl/GSM3069442_ST4_DGE.txt.gz","ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3069nnn/GSM3069443/suppl/GSM3069443_ST5_DGE.txt.gz","ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3069nnn/GSM3069444/suppl/GSM3069444_ST6_DGE.txt.gz"),organism="mmusculus")
