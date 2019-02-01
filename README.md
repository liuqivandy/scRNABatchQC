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
        
	
	
	setwd("/scratch/scRNABatchQC/")     # set working directory
	library(scRNABatchQC)
	
	scRNABatchQC(inputfiles=c("https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar1.csv.gz", "https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar5.csv.gz"))
