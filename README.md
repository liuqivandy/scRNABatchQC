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

Before you install scRNABatchQC, a modified version of WebGestatR is highly recommended to be installed by:

	library(devtools)
	install_github("shengqh/WebGestaltR")

Then you can install scRNABatchQC by:

	install_github("liuqivandy/scRNABatchQC")
  
<a name="example"/>

# Quick start

Here we show the most basic steps.
        
	
	
	setwd("/scratch/scRNABatchQC/")     # set working directory
	library(scRNABatchQC)
	organism = "mmusculus"
	sampleTable <- data.frame(Sample = c("S1", "S2", "S3"), File = c("count1.csv", "count2.csv", "count3.csv"))
	scRNABatchQC(sampleTable, organism, "scRNABatchQCreport.html", cache=TRUE )
