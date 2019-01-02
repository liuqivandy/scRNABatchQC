scRNABatchQC
==========
* [Introduction](#introduction)
* [Download and installation](#download)
* [Quick Start](#example)

<a name="introduction"/>

# Introduction

scRNABatchQC is an R package for quality control of multiple single cell RNAseq data.

<a name="download"/>

# Download and installation

Step 1:  Install pandoc (https://github.com/jgm/pandoc/releases/tag/2.2.1) . After installation , update your path to include the directory where pandoc’s binaries are installed.

Step 2: Install a modified version of WebGestatR by:

	library(devtools)
	install_github("shengqh/WebGestaltR")

Step 3: Install scRNABatchQC by:

	install_github("liuqivandy/scRNABatchQC")
  
<a name="example"/>

# Quick start

Here we show the most basic steps.
        
	
	
	setwd("/scratch/scRNABatchQC/")     # set working directory
	library(scRNABatchQC)
	organism = "mmusculus"
	sampleTable <- data.frame(Sample = c("S1", "S2", "S3"), File = c("count1.csv", "count2.csv", "count3.csv"))
	scRNABatchQC(sampleTable, organism, "scRNABatchQCreport.html", cache=TRUE )
