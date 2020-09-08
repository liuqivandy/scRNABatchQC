#' QC across multiple scRNAseq datasets
#' 
#' @description Compare multiple scRNA-seq datasets simultaneously on numerous technical and biological features 

#' @param inputs string vector of file or path names,  or a list of SingleCellExperiment or Seurat v3 objects; \cr
#' inputs can be a vector of file names (or a URL starting http://, file://, etc.) of gene-by-cell count matrices, the rowname should be gene symbol; each file should be regular delimited file;  Compressed files ending .gz and .bz2 are supported.\cr
#' inputs can be a vector of path names, each of which contains the barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz provided by 10X from CellRanger >=3.0 \cr
#' inputs can also be a list of SingleCellExperiment or Seurat v3 objects
#' @param names string vector; giving names of each sample  (default: NULL); names should have the same length of inputs; if NULL, the names are S1, S2... 
#' @param nHVGs integer; the number of highly variable genes (default: 1000)
#' @param nPCs integer; the number of principal components (default: 10)
#' @param sf integer; Scale factor to normalize the single cell RNA-seq data (default: 10000)
#' @param mincounts integer; the cutoff of filtering the cell if the total number of counts in the cell less than the mincounts (default:500)
#' @param mingenes integer; the cutoff of filtering the cell if the total number of genes detected in the cell less than the mingenes (default: 200)
#' @param maxmito float; the cutoff of filtering the cell if the percentage of mtRNA reads in the cell larger than the minmito (default: 0.2) 
#' @param PCind integer; which principal component for exploring biological featues (default: 1; the first principal component will be used to find genes highly correlated with PCA 1); PCind should be less than nPC 
#' @param mtRNA string; the pattern of genenames for mitochondrial encoded RNAs ; (default: "^mt-|^MT-", the default is mtRNA genenames in human or mouse); If not human or mouse, give the gene name pattern of mtRNA
#' @param rRNA string; the pattern of genenames for ribosomal proteins; (default: "^Rp[sl][[:digit:]]|^RP[SL][[:digit:]]", the default is ribosomal protein genenames in human or mouse); If not human or mouse, give the gene name pattern of ribosomal proteins
#' @param logFC float; log fold change cutoff to select differentially expressed genes (default: 1)
#' @param FDR float; FDR cutoff to select differentially expressed genes (default: 0.01)
#' @param sampleRatio float; the ratio of cells sampled from each dataset to examine the expression similarity (default: 1)
#' @param organism  string; the organism of single cell RNAseq datasets; if supported by WebGestaltR, functional enrichment analysis will be performed (default: mmusculus);WebGestaltR supports 12 organisms, including athaliana, btaurus,celegans, cfamiliaris, drerio, sscrofa, dmelanogaster, ggallus, hsapiens, mmusculus, rnorvegicus, and scerevisiae.
#' @param outputFile  string; the name of the output file (default: report.html)
#' @param lineSize float; the line size of figures in the report (default: 1)
#' @param pointSize float; the point size of figures in the report  (default: 0.8)
#' @param chunk.size NULL or integer; default is NULL, suggesting data will be loaded into memory at one time, otherwise, the data will be loaded into memory by chunks with chunk.size
#' @param createReport logical; default is TRUE, suggesting html report file will be created
#' @return a list of SingleCellExperiment objects;
#'  \itemize{
#'  \item  {       sces: a list of SingleCellExperiment objects; each object contains technical and biological metadata for one scRNAseq dataset; see the output of  \code{\link{Process_scRNAseq} }}
#'  \item {        scesMerge: a SingleCellExperiment object containing bilogical metadata for the combined dataset and the pairwise difference across datasets; see the output of \code{\link{Combine_scRNAseq} }}
#' }
#' @examples
#' library(scRNABatchQC)  
#' output<-scRNABatchQC(inputs=c("https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar1.csv.gz","https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar5.csv.gz"))
#' 
#' # a list of SingleCellExperiment objects, one for each individual dataset
#' output$sces
#' # a SingleCellExperiment object for the combined dataset
#' output$scesMerge
#' 
#' plotDensity(output$sces, "total_counts")
#' output$sces[[1]]@metadata$hvgPathway
#' plotHVGs(output$sces)
#' output$scesMerge@metadata$diffFC$genes
#' 
#' #scRNABatchQC can run on a list of SingleCellExperiment objects
#' scRNABatchQC(inputs=output$sces)
#' 
#' #scRNABatchQC can run on a list of Seurat v3 objects
#' library(Seurat)
#' S1<-CreateSeuratObject(counts=counts(output$sces[[1]]))
#' S2<-CreateSeuratObject(counts=counts(output$sces[[2]]))
#' scRNABatchQC(inputs=list(S1,S2))
#' 
#' @import R.utils ggplot2 gplots limma data.table irlba Rtsne WebGestaltR rmdformats Matrix statmod DT RCurl SingleCellExperiment
#' @importFrom devtools session_info
#' 
#' @export
#' @seealso \code{\link{Process_scRNAseq}}, \code{\link{Combine_scRNAseq}} , \code{\link{generateReport}}
scRNABatchQC<-function(inputs,names=NULL, nHVGs=1000,nPCs=10,sf=10000,mincounts=500,mingenes=200, maxmito=0.2, PCind=1, mtRNA="^mt-|^MT-", rRNA="^Rp[sl][[:digit:]]|^RP[SL][[:digit:]]", sampleRatio=1, logFC=1,FDR=0.01, organism="mmusculus", outputFile="report.html", lineSize=1, pointSize=0.8,chunk.size=NULL, createReport=TRUE){
  isOrganismValid<-.isOrganismValid(organism)
  if(! isOrganismValid){
    organism<-NULL
  }
  sces<-Process_scRNAseq(inputs=inputs,names=names,nHVGs=nHVGs, nPCs=nPCs,sf=sf, mincounts=mincounts, mingenes=mingenes, maxmito=maxmito,PCind=PCind,mtRNA=mtRNA, rRNA=rRNA, organism=organism, chunk.size=chunk.size)
  scesMerge<-Combine_scRNAseq(sces,nHVGs=nHVGs, nPCs= nPCs, logFC=logFC,FDR=FDR,sampleRatio=sampleRatio,organism=organism)
  if(createReport){
    generateReport(sces,scesMerge, outputFile=outputFile, lineSize=lineSize, pointSize=pointSize)
  }
  return(list(sces = sces, 
              scesMerge = scesMerge ))
}

###############################
#' process scRNAseq datasets one by one to generate QC metadata
#' @description  Generate technical and biological metadata for one or multiple single-cell RNAseq datasets;each dataset is processed one by one
#' @param inputs a string vector of file or path names, or a list of SingleCellExperiment or Seurat v3 objects; \cr
#' inputs can be a string vector of file names (or a URL starting http://, file://, etc.) of gene-by-cell count matrices, the rowname should be gene symbol; each file should be regular delimited file;  Compressed files ending .gz and .bz2 are supported. \cr
#' inputs can be a string vector of path names, each of which contains the barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz provided by 10X from CellRanger >=3.0 \cr
#' inputs can also be a list of SingleCellExperiment or Seurat v3 objects
#' @param names string vector; giving the names of single-cell RNAseq datasets (default: NULL); names should have the same length of inputs; if NULL, the names are S1, S2... 
#' @param nHVGs integer; the number of highly variable genes (default: 1000)
#' @param nPCs integer; the number of principal components (default: 10)
#' @param sf integer; Scale factor to normalize the data  (default: 10000)
#' @param mincounts integer; the cutoff of filtering the cell if the total number of counts in the cell less than the mincounts (default:500)
#' @param mingenes integer; the cutoff of filtering the cell if the total number of genes detected in the cell less than the mingenes (default: 200)
#' @param maxmito float; the cutoff of filtering the cell if the percentage of mtRNA reads in the cell larger than the minmito (default: 0.2)
#' @param PCind integer; which principal component for exploring biological featues (default: 1; the first principal component will be used to find genes highly correlated with PCA 1); PCind should be less than nPC 
#' @param mtRNA string; the pattern of gene names for mitochondrial encoded RNAs ; (default: "^mt-|^MT-", the default is mtRNA gene names in human or mouse); If not human or mouse, give the gene name pattern of mtRNA
#' @param rRNA string; the pattern of gene names for ribosomal proteins; (default: "^Rp[sl][[:digit:]]|^RP[SL][[:digit:]]", the default is ribosomal protein gene names in human or mouse); If not human or mouse, give the gene name pattern of ribosomal proteins
#' @param organism string; the organism of single cell RNAseq datasets; if supported by WebGestaltR, functional enrichment analysis will be performed (defeault: mmusculus) 
#' @param chunk.size NULL or integer; default is NULL, suggesting data will be loaded into memory at one time, otherwise, the data will be loaded into memory by chunks with chunk.size

#' @return a list of SingleCellExperiment objects ;  \cr
#'  each SingleCellExperiment object containing metadata for one single cell RNAseq dataset;  \cr
#'         each SingleCellExperiment object containing several slots:
#' \itemize{
#' \item  {                        assays; ShallowSimpleListAssays object containing two sparse matrix: counts and logcounts }
#' \item  {                        rowRanges@elementMetadata; A DataFrame containining metadata for each gene, including 
#'            \itemize{  
#'                         \item {              ave.counts: the average counts  }
#'                         \item {              num.cells: the number of cells with the gene detected }
#'                         \item {              hvg: a dataframe containing mean, variance and z-score for dispersion }
#'                         \item {              genevar_by_tounts:  variance explained by the log-transformed counts }
#'                         \item {              genevar_by_features: variance explained by the log-transformed features }
#'                         \item {              genevar_by_Mt: variance explained by the log-transformed mitochondrial counts }
#'                         \item {              genevar_by_rRNA: variance explained by the log-transformed rRNA counts }
#'              }
#'         }
#'  \item {                        colData; A DataFrame containing metadata for each cell, including
#'               \itemize{
#'                           \item {                log10_total_counts: the number of counts in each cell in log10 transformed }
#'                           \item {               log10_total_features: the number of genes detected in each cell in log10 transformed }
#'                           \item {               log10_total_counts_Mt: the number of mitochondrial counts in each cell in log10 transformed }
#'                           \item {               log10_total_counts_rRNA: the number of rRNA counts in each cell in log10 transformed  }
#'                }
#'        }
#'  \item {                               metadata: A list containing other metadata, including
#'             \itemize{
#'                          \item {         rawmeta: a list of metadata for genes and cells of raw data before filtering, including
#'                                   \itemize{ 
#'                                          \item {        sf: normalization factor }
#'                                          \item {        ngenes: the number of genes }
#'                                          \item {        ncells: the number of cells }
#'                                          \item {        CellData: a dataframe containing metadata for each cell before filtering, includingtotal_counts,total_features,total_counts_Mt,total_counts_rRNA, pct_counts_rRNA (perentage of rRNA counts), pct_counts_Mt (percentage of mtRNA counts), libsize.drop(cell is filtered by library size), feature.drop (cell is filtered by the number of detected genes), mito.drop (cell is filtered by the mtRNA counts), is.drop (cell is filtered by either of library size, the number of genes or the mtRNA reads) }
#'                                          \item {        GeneData: a list containing the gene filter information (gene.keep, gene is filtered since none of the cells detect theg gene)    }
#'	                                        \item {        Cutoff: a list containing the cutoff values (count, gene, mito) for filtering cells                                                }   
#'                                         }
#'                                  }
#'                         \item {           pc1genes: a dataframe containing the genes highly correlated with the 1st (default) or the PCind principal component }
#'                         \item {           pc1Pathway: a dataframe containing the pathways enriched in pc1 (default) or the PCind genes  }
#'                         \item {           hvgPathway: a dataframe containing the pathways enriched in top n (default:1000) highly variable genes }
#'                        }                                            
#'     }
#'  }                              
#'                                

#' @examples
#' library(scRNABatchQC)
#' sces<-Process_scRNAseq(inputs=c("https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar1.csv.gz","https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar5.csv.gz"))
#' names(sces)
#' class(sces[[1]])
#' head(sces[[1]]@rowRanges@elementMetadata)
#' head(sces[[2]]@rowRanges@elementMetadata)
#' head(colData(sces[[1]]))
#' sces[[1]]@metadata$rawmeta$ngenes
#' head(sces[[1]]@metadata$rawmeta$CellData)
#' head(sces[[1]]@metadata$pc1Pathway)
#' plotDensity(sces,"total_counts")
#' plotVarianceTrend(sces)
#' plotPCPathways(sces)
#' @import R.utils ggplot2 gplots limma data.table irlba Rtsne WebGestaltR rmdformats Matrix statmod DT RCurl SingleCellExperiment
#' @importFrom devtools session_info
#' 
#' @export
#' @seealso \code{\link{Combine_scRNAseq}} , \code{\link{generateReport}}
Process_scRNAseq <- function(inputs, names=NULL, nHVGs=1000, nPCs=10,sf=10000,mincounts=500,mingenes=200, maxmito=0.2,PCind=1,mtRNA="^mt-|^MT-", rRNA="^Rp[sl][[:digit:]]|^RP[SL][[:digit:]]", organism="mmusculus",chunk.size=NULL){
  isOrganismValid<-.isOrganismValid(organism)
  if(! isOrganismValid){
    organism<-NULL
  }
  
  result <- list()
  nfiles<-length(inputs)
  
  if (is.null(names)){names=paste0("S",1:nfiles)}
  if (nfiles!=length(names)) {stop("the inputs and names should have the same length", call. = FALSE)}
  if (sum(names!=make.names(names))>0) names<-paste0("S",names)
  
  for (ind in 1:nfiles) {
    cat("Processing ", names[ind], "\n")
    if (is.list(inputs)) {
          result[[ind]] <- Process_OnescRNAseq(inputs[[ind]], nHVGs=nHVGs, nPCs=nPCs,sf=sf,mincounts=mincounts,mingenes=mingenes, maxmito=maxmito,PCind=PCind,mtRNA=mtRNA, rRNA=rRNA, organism=organism)
    } else {
          result[[ind]] <- Process_OnescRNAseq(inputs[ind],  nHVGs=nHVGs, nPCs=nPCs,sf=sf,mincounts=mincounts,mingenes=mingenes, maxmito=maxmito,PCind=PCind,mtRNA=mtRNA, rRNA=rRNA, organism=organism, chunk.size=chunk.size)
    }
  }
  
  names(result) <- names
  return(result)
}

########################################

#' combine multiple scRNAseq datasets into one and generate QC metadata
#' @description Combine and compare multiple single cell RNAseq datasets, including \cr
#'  1)combine into one dataset; \cr
#'  2)find highly variable genes, reduce dimensions using PCA and tSNE for the combined dataset;  \cr
#'  3)pairwise comparsion between any two datasets to find the top differentially expressed genes; \cr 
#'
#' @param sces a list of SingleCellExperiment objects (results from \code{\link{Process_scRNAseq}}); each object corresponds to one single cell RNAseq dataset 
#' @param nHVGs integer; the number of highly variable genes (default:1000)
#' @param nPCs integer; the number of principal components (default:10)
#' @param logFC float; log fold change cutoff to select differentially expressed genes (default: 1)
#' @param FDR float; FDR cutoff to select differentially expressed genes (default: 0.01)
#' @param sampleRatio float; the ratio of cells sampled from each dataset to examine the expression similarity(default: 1)
#' @param organism string; the organism of single cell RNAseq datasets;if supported by WebGestaltR, the functional enrichment will be performed; (defeault: mmusculus) 
#' @return a SingleCellExperiment object with several slots:
#' \itemize{
#' 	               \item {         assays; ShallowSimpleListAssays object containing one sparse matrix  logcounts (log-transformed normalized counts) }
#'                 \item {         rowRanges@elementMetadata; A Dataframe hvg containining mean, variance and z-score for each gene  }
#'                 \item {         colData; A Dataframe containing conidtion for each cell   }
#'                 \item {         metadata: A list containing other metadata, including 
#'                                 \itemize{
#'                                            \item {    reducedDims: a list containing two types of reduced dimensions: PCA and tSNE  }
#' 	                                          \item {    diffFC: differential expressed genes and pathways from pairwise comparison between any two datasets }
#'                                            \item {    other metadata including nHVGs,nPC, logFC, FDR, sampleRatio  }
#'                                         }
#'                          }
#'          }
#' @importFrom Rtsne Rtsne
#' @importFrom Matrix Matrix
#' @import SingleCellExperiment
#' @export 
#' @examples 
#' library(scRNABatchQC)
#' sces<-Process_scRNAseq(inputs=c("https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar1.csv.gz","https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar5.csv.gz"))
#' scesMerge <- Combine_scRNAseq(sces)
#' logcounts(scesMerge)[1:5,1:5]
#' head(scesMerge@rowRanges@metadata$hvg)
#' summary(scesMerge@metadata$reducedDims$PCA)
#' #visualize PCA results
#' plot(scesMerge@metadata$reducedDims$PCA$x[,1:2],pch=16,col=as.factor(scesMerge@colData$condition),xlab="PCA1",ylab="PCA2")
#' #visualize tSNE results
#' plot(scesMerge@metadata$reducedDims$tSNE,pch=16,col=as.factor(scesMerge@colData$condition),xlab="tSNE1",ylab="tSNE2")
#' scesMerge@metadata$diffFC
#' @seealso \code{\link{Process_scRNAseq}} , \code{\link{generateReport}}

Combine_scRNAseq <- function(sces, nHVGs=1000, nPCs= 10, logFC=1,FDR=0.01,sampleRatio=1,organism="mmusculus") {
  if (sampleRatio>1) stop("sampleRatio must be <=1")
  
  isOrganismValid<-.isOrganismValid(organism)
  if(! isOrganismValid){
    organism<-NULL
  }
  cat("Merging data.\n")
  pca_tsne_data <- list()
  
  if(length(sces) > 0){
    for (i in 1:length(sces)) {
      sce<-sces[[i]]
      if (!sce@metadata$valid){
        next
      }

      nsample=round(dim(sces[[i]])[2]*sampleRatio,0)
      sampleind<-sample(1:dim(sces[[i]])[2],nsample)
      mat<-logcounts(sces[[i]])[,sampleind]
      colnames(mat)<-paste0(names(sces)[i], 1:nsample)
      condition<-rep(names(sces)[i],nsample)
      
      if(length(pca_tsne_data) == 0){
        pca_tsne_data$logcounts<-mat
        pca_tsne_data$condition <-condition
      }else{
        pca_tsne_data$logcounts <- .mergeSparseMatrix(pca_tsne_data$logcounts, mat)
        pca_tsne_data$condition <-c(pca_tsne_data$condition,condition)
      }
    }
  }
  
  if(length(pca_tsne_data) > 0){
    pca_tsne_data$hvg <- .getMeanVarTrend(pca_tsne_data$logcounts)
    
    ##select the top 1000 highly variable genes for the PCA
    feature_set <-  rownames(pca_tsne_data$hvg)[order(pca_tsne_data$hvg$zval,decreasing=T)][1:nHVGs]
    
    #scevar <- apply(scesdata, 1, var)
    #feature_set <- head(order(scevar, decreasing = T), n = 500)
    
    tdata<-t(pca_tsne_data$logcounts[rownames(pca_tsne_data$logcounts)%in%feature_set, , drop = FALSE])
    nPCs<-min(nPCs, min(dim(tdata))-1)
    
    pca_tsne_data$pca <- prcomp_irlba(tdata, n= nPCs)
  
    set.seed(100)
    perplexity<-min(20, floor((nrow(pca_tsne_data$pca$x)-1)/3))
    tsne_out <- Rtsne(pca_tsne_data$pca$x, initial_dims = ncol(pca_tsne_data$pca$x), pca = FALSE, perplexity =perplexity, check_duplicates = FALSE)
    pca_tsne_data$tsne <- tsne_out$Y
    scesMerge<-SingleCellExperiment(assay=list(logcounts=pca_tsne_data$logcounts))
    scesMerge@metadata$reducedDims=list(PCA=pca_tsne_data$pca, tSNE=pca_tsne_data$tsne)
    scesMerge@colData$condition<-pca_tsne_data$condition
    scesMerge@rowRanges@metadata$hvg<-pca_tsne_data$hvg
    
    #compare conditions
    cat("Performing differential expression analysis data ...\n")
    scesMerge@metadata$diffFC <- .getDiffGenes(scesMerge, organism = organism,  logFC=logFC, FDR = FDR, geneNo = 50)
    scesMerge@metadata$logFC<- logFC
    scesMerge@metadata$FDR<-FDR
    scesMerge@metadata$sampleRatio<-sampleRatio
    scesMerge@metadata$nHVGs<-nHVGs
    scesMerge@metadata$nPCs<-nPCs
  }
  else{
    scesMerge<-NULL
  }
  return(scesMerge)
}

#################
#' generate a QC report in a html file
#' @description Generate a QC report by comparing technical and biological features across multiple scRNA-seq datasets
#' @param sces a list of SingleCellExperiment Objects; each object containing quality meta data for one single cell RNAseq dataset (results from \code{\link{Process_scRNAseq}})
#' @param scesMerge a SingleCellExperiment object generated by combining and comparing multiple single cell RNAseq datasets (results from \code{\link{Combine_scRNAseq}})
#' @param outputFile the name of the output file (default: report.html)
#' @param lineSize float; the line size of figures in the generated report (default: 1)
#' @param pointSize float; the point size of figures in the generated report (default: 0.8)
#' @examples
#' library(scRNABatchQC)
#' sces<-Process_scRNAseq(inputs=c("https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar1.csv.gz","https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar5.csv.gz"))
#' scesMerge <- Combine_scRNAseq(sces)
#' generateReport(sces,scesMerge)
#' @import R.utils ggplot2 gplots limma data.table irlba Rtsne WebGestaltR rmdformats Matrix statmod DT RCurl SingleCellExperiment knitr
#' @importFrom devtools session_info
#' 
#' @export

#' @seealso \code{\link{Process_scRNAseq}} \code{\link{Combine_scRNAseq}}
generateReport<-function(sces, scesMerge, outputFile="report.html", lineSize=1, pointSize=0.8) {
  pw <- .prepareTableSummary(sces)
  plotData <- list(sces = sces, 
                   scesMerge = scesMerge, 
                   tableSummary = pw,
                   lineSize=lineSize,
                   pointSize=pointSize)
  
  cat("Report html generated.\n")
  reportRmd <- system.file("report/scRNABatchQCreport.Rmd", package="scRNABatchQC")
  # reportRmd<-"d:/github/scRNABatchQC/inst/report/scRNABatchQCreport.Rmd"
  
  outputFile <- getAbsolutePath(outputFile)
  output_dir = dirname(outputFile)
  output_file = basename(outputFile)
  
  cat("Output report to:", outputFile, "\n")
  rmarkdown::render(reportRmd,
                    output_dir = output_dir,
                    output_file = output_file,
                    intermediates_dir=tempdir(),
                    params = list(data = plotData))
}

##############################################
#' process one scRNAseq dataset to generate QC metadata 
#' @description Generate technical and biological metadata for one single cell RNAseq dataset
#' @param input string of file or path name, a SingleCellExperiment or Seurat v3 object; \cr
#' input can be the file name (or a URL starting http://, file://, etc.) of gene-by-cell count matrix, the rowname should be gene symbol; the file should be regular delimited file;  Compressed files ending .gz and .bz2 are supported. \cr
#' input can be the path name, which contains the barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz provided by 10X from CellRanger >=3.0 \cr
#' input can also be a SingleCellExperiment or Seurat v3 object
#' @param nHVGs integer; the number of highly variable genes (default: 1000)
#' @param nPCs integer: the number of principal components (default: 10)
#' @param sf integer; Scale factor to normalize the single cell RNA-seq data (default: 10000)
#' @param mincounts integer; the cutoff of filtering the cell if the total number of counts in the cell less than the mincounts (default:500)
#' @param mingenes integer; the cutoff of filtering the cell if the total number of genes detected in the cell less than the mingenes (default: 200)
#' @param maxmito float; the cutoff of filtering the cell if the percentage of mtRNA reads in the cell larger than the minmito; (default: 0.2); 
#' @param PCind integer; which principal component for exploring biological featues (default: 1; the first principal component will be used to find genes highly correlated with PCA 1); PCind should be less than nPC 
#' @param mtRNA string; the pattern of genenames for mitochondrial encoded RNAs ; (default: "^mt-|^MT-", the default is mtRNA genenames in human or mouse); If not human or mouse, input the gene name pattern of mtRNA
#' @param rRNA string; the pattern of genenames for ribosomal proteins; (default: "^Rp[sl][[:digit:]]|^RP[SL][[:digit:]]", the default is ribosomal protein genenames in human or mouse); If not human or mouse, input the gene name pattern of ribosomal proteins
#' @param organism string; the organism of single cell RNAseq datasets; if supported by WebGestaltR, functional enrichment analysis will be performed (defeault: mmusculus) 
#' @param chunk.size NULL or integer; default is NULL, suggesting data will be loaded into memory at one time, otherwise, the data will be loaded into memory by chunks with chunk.size

#' @return a SingleCellExperiment object with several slots:
#' \itemize{
#' \item  {                            assays; ShallowSimpleListAssays object containing two sparse matrix: counts and logcounts }
#' \item  {                            rowRanges@elementMetadata; A DataFrame containining metadata for each gene, including 
#'            \itemize{  
#'                         \item {               ave.counts: the average counts  }
#'                         \item {               num.cells: the number of cells with the gene detected }
#'                         \item {               hvg: a dataframe containing mean, variance and z-score for dispersion }
#'                         \item {               genevar_by_tounts:  variance explained by the log-transformed counts }
#'                         \item {               genevar_by_features: variance explained by the log-transformed features }
#'                         \item {               genevar_by_Mt: variance explained by the log-transformed mitochondrial counts }
#'                         \item {               genevar_by_rRNA: variance explained by the log-transformed rRNA counts }
#'              }
#'         }
#'  \item {                             colData; A DataFrame containing metadata for each cell, including
#'               \itemize{
#'                           \item{                log10_total_counts: the number of counts in each cell in log10 transformed }
#'                           \item {               log10_total_features: the number of genes detected in each cell in log10 transformed }
#'                           \item {               log10_total_counts_Mt: the number of mitochondrial counts in each cell in log10 transformed }
#'                           \item {               log10_total_counts_rRNA: the number of rRNA counts in each cell in log10 transformed  }
#'                }
#' }
#'  \item {                               metadata: A list containing other metadata, including
#'                \itemize{
#'                              \item {         rawmeta: a list of metadata for genes and cells of raw data before filtering, including
#'                                    \itemize{ 
#'                                          \item {        sf: normalization factor }
#'                                          \item {        ngenes: the number of genes }
#'                                          \item {        ncells: the number of cells }
#'                                          \item {        CellData: a dataframe containing metadata for each cell before filtering, includingtotal_counts,total_features,total_counts_Mt,total_counts_rRNA, pct_counts_rRNA (perentage of rRNA counts), pct_counts_Mt (percentage of mtRNA counts), libsize.drop(cell is filtered by library size), feature.drop (cell is filtered by the number of detected genes), mito.drop (cell is filtered by the mtRNA counts), is.drop (cell is filtered by either of library size, the number of genes or the mtRNA reads) }
#'                                          \item {        GeneData: a list containing the gene filter information (gene.keep, gene is filtered since none of the cells detect theg gene)    }
#'	                                        \item {        Cutoff: a list containing the cutoff values (count, gene, mito) for filtering cells                                                }   
#'                                         }
#'                                  }
#'                         \item {           pc1genes: a dataframe containing the genes highly correlated with the 1st (default) or the PCind principal component }
#'                         \item {           pc1Pathway: a dataframe containing the pathways enriched in pc1 (default) or the PCind genes  }
#'                         \item {           hvgPathway: a dataframe containing the pathways enriched in top n (default:1000) highly variable genes }
#'                        }                                            
#'     }
#'  }                              
#'                              
#' @examples
#' library(scRNABatchQC)
#' sce<-Process_OnescRNAseq(input="https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar1.csv.gz")
#' head(sce@rowRanges@elementMetadata)
#' head(colData(sce))
#' counts(sce)[1:5,1:5]
#' logcounts(sce)[1:5,1:5]
#' sce@metadata$rawmeta$ngenes
#' head(sce@metadata$rawmeta$CellData)
#' head(sce@metadata$pc1Pathway)
#' sces=list(sce=sce)
#' plotDensity(sces)
#' @import R.utils ggplot2 gplots limma data.table irlba Rtsne WebGestaltR rmdformats Matrix statmod DT RCurl SingleCellExperiment
#' @importFrom devtools session_info
#' 
#' @export
#' @seealso \code{\link{Process_scRNAseq}}, \code{\link{Combine_scRNAseq}}
Process_OnescRNAseq <- function(input, sf=10000,mincounts=500,mingenes=200, maxmito=0.2,mtRNA="^mt-|^MT-", rRNA="^Rp[sl][[:digit:]]|^RP[SL][[:digit:]]", nHVGs=1000, nPCs=10,PCind=1, organism="mmusculus",chunk.size=NULL) {
  isOrganismValid<-.isOrganismValid(organism)
  if(! isOrganismValid){
    organism<-NULL
  }
  
  sce<-Tech_OnescRNAseq(input=input,sf=sf,mincounts=mincounts,mingenes=mingenes,maxmito=maxmito,mtRNA=mtRNA, rRNA=rRNA, chunk.size=chunk.size)
  if (ncol(sce) > 10){
    sce<-Bio_OnescRNAseq(sce,nHVGs=nHVGs,nPCs=nPCs,PCind=PCind,organism=organism)
    sce@metadata$valid=TRUE
  }else{
    sce@metadata$valid=FALSE
  }
  return(sce)
}



#####################################
#' process one scRNAseq dataset to generate technical metadata
#' @description Generate technical metadata for one single cell RNAseq dataset
#' @param input string of file or path name, a SingleCellExperiment or Seurat v3 object; \cr
#' input can be a string of the file name (or a URL starting http://, file://, etc.) of gene-by-cell count matrix, the rowname should be gene symbol; the file should be regular delimited file;  Compressed files ending .gz and .bz2 are supported. \cr 
#' input can be a string of the path name, which contains the barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz provided by 10X from CellRanger >=3.0 \cr
#' input can also be a SingleCellExperiment or a Seurat v3 object
#' @param sf integer; Scale factor to normalize the single cell RNA-seq data (default: 10000)
#' @param mincounts integer; the cutoff of filtering the cell if the total number of counts in the cell less than the mincounts (default:500)
#' @param mingenes integer; the cutoff of filtering the cell if the total number of genes detected in the cell less than the mingenes (default: 200)
#' @param maxmito  float; the cutoff of filtering the cell if the percentage of mtRNA reads in the cell larger than the minmito; (default: 0.2); 
#' @param mtRNA string; the pattern of genenames for mitochondrial encoded RNAs ; (default: "^mt-|^MT-", the default is mtRNA genenames in human or mouse); If not human or mouse, input the gene name pattern of mtRNA
#' @param rRNA string; the pattern of genenames for ribosomal proteins; (default: "^Rp[sl][[:digit:]]|^RP[SL][[:digit:]]", the default is ribosomal protein genenames in human or mouse); If not human or mouse, input the gene name pattern of ribosomal proteins
#' @param chunk.size NULL or integer; default is NULL, suggesting data will be loaded into memory at one time, otherwise, the data will be loaded into memory by chunks with chunk.size

#' @return a SingleCellExperiment object containing metadata for technical features
#' \itemize{
#'  \item  {                            assays; ShallowSimpleListAssays object containing two sparse matrix: counts and logcounts }
#'  \item  {                            rowRanges@elementMetadata; A DataFrame containining metadata for each gene, including 
#'            \itemize{  
#'                         \item {               ave.counts: the average counts  }
#'                         \item {               num.cells: the number of cells with the gene detected }
#'              }
#'         }
#'  \item {                             colData; A DataFrame containing metadata for each cell, including
#'               \itemize{
#'                           \item{                log10_total_counts: the number of counts in each cell in log10 transformed }
#'                           \item {               log10_total_features: the number of genes detected in each cell in log10 transformed }
#'                           \item {               log10_total_counts_Mt: the number of mitochondrial counts in each cell in log10 transformed }
#'                           \item {               log10_total_counts_rRNA: the number of rRNA counts in each cell in log10 transformed  }
#'                }
#'       }
#'  \item {                               metadata: A list containing other metadata, including
#'                \itemize{
#'                              \item {         rawmeta: a list of metadata for genes and cells of raw data before filtering, including
#'                                    \itemize{ 
#'                                          \item {        sf: normalization factor }
#'                                          \item {        ngenes: the number of genes }
#'                                          \item {        ncells: the number of cells }
#'                                          \item {        CellData: a dataframe containing metadata for each cell before filtering, includingtotal_counts,total_features,total_counts_Mt,total_counts_rRNA, pct_counts_rRNA (perentage of rRNA counts), pct_counts_Mt (percentage of mtRNA counts), libsize.drop(cell is filtered by library size), feature.drop (cell is filtered by the number of detected genes), mito.drop (cell is filtered by the mtRNA counts), is.drop (cell is filtered by either of library size, the number of genes or the mtRNA reads) }
#'                                          \item {        GeneData: a list containing the gene filter information (gene.keep, gene is filtered since none of the cells detect theg gene)    }
#'	                                        \item {        Cutoff: a list containing the cutoff values (count, gene, mito) for filtering cells                                                }   
#'                                         }
#'                                  }
#'                           }
#'    }
#' }
#' @import R.utils ggplot2 gplots limma data.table irlba Rtsne WebGestaltR rmdformats Matrix statmod DT RCurl SingleCellExperiment
#' @importFrom devtools session_info
#' 
#' @export
#' @examples
#' library(scRNABatchQC)
#' sce<-Tech_OnescRNAseq(input="https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar1.csv.gz")
#' plotDensity(list(sce=sce))
#' names(sce@metadata)
#' head(sce@colData$log10_total_counts)
#' @seealso \code{\link{Process_OnescRNAseq}} , \code{\link{Bio_OnescRNAseq}} 


Tech_OnescRNAseq<-function(input, sf=10000,mincounts=500,mingenes=200, maxmito=0.2,mtRNA="^mt-|^MT-", rRNA="^Rp[sl][[:digit:]]|^RP[SL][[:digit:]]", chunk.size=NULL ){
 if (is(input, "SingleCellExperiment")){
     countmat<-counts(input)
  } else if (is(input, "Seurat")) {
    if (!require("Seurat",character.only = TRUE)) { stop("Please install Seurat")}
    countmat<-GetAssayData(object = input, assay = "RNA", slot = "counts")
  } else if(dir.exists(input)){
    if(.is_10X_v3(input)) {
      countmat<-read_10X_v3(input)
    }else if(.is_10X_v2(input)) {
      countmat<-read_10X_v2(input)
    }else{
      stop(paste0("Input folder doesn't contain 10X v3/v2 files: ", input))
    }
  } else if(endsWith(input, ".h5")){
    if (!require("Seurat",character.only = TRUE)) { stop("Please install Seurat")}
    countmat<-Read10X_h5(input)
    if (is.list(countmat)){
      cat("h5 file has multiple data entries, select Gene Expression\n")
      countmat<-countmat$`Gene Expression` 
    }
  } else {
    if(is.null(chunk.size)) {
      rawdata<-fread(input,data.table=F)
      countmat<-.tosparse(rawdata[,-1])
      rownames(countmat)<-rawdata[,1]
      rm(rawdata)
      gc()
    } else { countmat<-fread_bychunk(input,chunk.size=chunk.size)}
  } 
  
  #if (!is.integer(PCind) | !is.integer(nPCs) | ! is.integer(nHVGs)) stop("nPCs, PCind and nHVGs should be integer")
  
  ##before cell and genes filtering
  rawmeta<-list()
  rawmeta$ngenes<-dim(countmat)[1]
  rawmeta$ncells<-dim(countmat)[2]
  
  is.mito <- grepl(mtRNA, rownames(countmat))
  is.rRNA<-grepl(rRNA,rownames(countmat))
  
  total_counts=Matrix::colSums(countmat)
  total_features=Matrix::colSums(countmat != 0)
  total_counts_Mt = Matrix::colSums(countmat[is.mito, ])
  total_counts_rRNA=Matrix::colSums(countmat[is.rRNA, ])  
  
  rawmeta$CellData<-data.frame(total_counts=total_counts, total_features=total_features,total_counts_Mt=total_counts_Mt,total_counts_rRNA=total_counts_rRNA,pct_counts_Mt = total_counts_Mt/total_counts,pct_counts_rRNA=total_counts_rRNA/total_counts)
  
  ##filter cells with less than 500 counts or 3 mad lower
  libsize.drop = .findOutlier(rawmeta$CellData$total_counts,  log=TRUE,type = "lower",lower.limit=mincounts)
  rawmeta$CellData<-cbind(rawmeta$CellData,libsize.drop=libsize.drop$filtered)
  #rawmeta$CellData<-cbind(rawmeta$CellData,libsize.drop = .findOutlier(rawmeta$CellData$total_counts,  log=TRUE,type = "lower"))
  ##filter cells with less than 100 genes or 3 mad lower
  feature.drop=.findOutlier(rawmeta$CellData$total_features, log=TRUE,type = "lower", lower.limit=mingenes)
  rawmeta$CellData <- cbind(rawmeta$CellData,feature.drop=feature.drop$filtered)
  
  ##filter cells with larger than 20% mtRNA or 3mad higher
  mito.drop= .findOutlier(rawmeta$CellData$pct_counts_Mt, type = "higher",upper.limit=maxmito)
  
  rawmeta$CellData <- cbind(rawmeta$CellData,mito.drop= mito.drop$filtered)
  rawmeta$CellData<- cbind(rawmeta$CellData,is.drop=(rawmeta$CellData$libsize.drop | rawmeta$CellData$feature.drop | rawmeta$CellData$mito.drop))
  
  num.cells <- Matrix::rowSums(countmat != 0)
  
  rawmeta$GeneData$gene.keep <- num.cells > 0
  rawmeta$Cutoff<-list(count=libsize.drop$cutoff,gene=feature.drop$cutoff, mito=mito.drop$cutoff)
  
  ##filtered counts 
  newcount <- countmat[rawmeta$GeneData$gene.keep, !rawmeta$CellData$is.drop]
  logcount<-newcount
  ####normalize to scale factor (sf), the default is 10000
  lib_size <- sf/rawmeta$CellData$total_counts[!rawmeta$CellData$is.drop]
  
  ###
  rowind<-logcount@i+1  
  colind<-findInterval(seq(logcount@x)-1,logcount@p[-1])+1
  logcount@x<-log2(logcount@x*lib_size[colind]+1) 
  
  filteredmeta<-list()
  filteredmeta$ngenes<-dim(newcount)[1]
  filteredmeta$ncells<-dim(newcount)[2]
  filteredmeta$mincounts<-mincounts
  filteredmeta$mingenes<-mingenes
  filteredmeta$maxmito<-maxmito

  scdata<- SingleCellExperiment(assay=list(counts=newcount,logcounts=logcount))
  scdata@colData$log10_total_counts<-log10(rawmeta$CellData$total_counts)[!rawmeta$CellData$is.drop]
  scdata@colData$log10_total_features <- log10(rawmeta$CellData$total_features)[!rawmeta$CellData$is.drop]
  scdata@colData$log10_total_counts_rRNA <- log10(rawmeta$CellData$total_counts_rRNA+1)[!rawmeta$CellData$is.drop]
  scdata@colData$log10_total_counts_Mt <- log10(rawmeta$CellData$total_counts_Mt+1)[!rawmeta$CellData$is.drop]
  
  scdata@rowRanges@elementMetadata$ave.counts <- Matrix::rowMeans(counts(scdata))
  scdata@rowRanges@elementMetadata$num.cells<- num.cells[rawmeta$GeneData$gene.keep]
  
  ##keep the orignial meta data 
  scdata@metadata$rawmeta<-rawmeta
  scdata@metadata$sf<-sf

  scdata@metadata$filteredmeta<-filteredmeta
  
  return(scdata)
}


###################################################
#' process one scRNAseq dataset to generate biological metadata
#' @description Generate metadata (HVGs, PC-related genes, pathways) for one single cell RNAseq represented by gene-cell count table

#' @param scdata a SingleCellExperiment object; (results from \code{\link{Tech_OnescRNAseq}}
#' @param nHVGs integer; the number of highly variable genes (default: 1000)
#' @param nPCs integer: the number of principal components (default: 10)
#' @param PCind integer; which principal component for exploring biological featues (default: 1; the first principal component will be used to find genes highly correlated with PCA 1); PCind should be less than nPC 
#' @param organism string; the organism of single cell RNAseq datasets; if supported by WebGestaltR, functional enrichment analysis will be performed (defeault: mmusculus) 
#' @return a new SingleCellExperiment object by adding features into the the elementMetadata and metadata slots of the input object:
#' 	                               
#'  \itemize{
#'               \item {         rowRanges@elementMetadata; 
#'                            \itemize{          
#'                                      \item {   hvg: a dataframe containing mean, variance and z-score for dispersion  }
#'                                      \item {   genevar_by_tounts:  variance explained by the log-transformed counts  }
#'                                      \item {   genevar_by_features: variance explained by the log-transformed features  } 
#'                                      \item {   genevar_by_Mt: variance explained by the log-transformed mitochondrial counts  }
#'                                      \item {   genevar_by_rRNA: variance explained by the log-transformed rRNA counts  }
#'                                   }
#'                      }
#'                \item {         metadata: A list containing other metadata, including
#'                                \itemize{
#'                                      \item {            pc1genes: a dataframe containing the genes highly correlated with the 1st (default) or the PCind principal component  }
#'                                      \item {            pc1Pathway: a dataframe containing the pathways enriched in pc1 (default) or the PCind genes  }
#'                                      \item {            hvgPathway: a dataframe containing the pathways enriched in top n (default:1000) highly variable genes  }
#'                                      }                              
#'                       }
#'  }
#' @import R.utils ggplot2 gplots limma data.table irlba Rtsne WebGestaltR rmdformats Matrix statmod DT RCurl SingleCellExperiment
#' @importFrom devtools session_info
#' 
#' @export
#' @examples
#' library(scRNABatchQC)
#' sce<-Tech_OnescRNAseq(input="https://github.com/liuqivandy/scRNABatchQC/raw/master/bioplar1.csv.gz")
#' sce<-Bio_OnescRNAseq(sce)
#' head(sce@rowRanges@elementMetadata)
#' sce@metadata$hvgPathway
#' @seealso \code{\link{Process_OnescRNAseq}} , \code{\link{Tech_OnescRNAseq}} 
Bio_OnescRNAseq<-function(scdata,nHVGs=1000, nPCs=10,PCind=1, organism="mmusculus" ){
  isOrganismValid<-.isOrganismValid(organism)
  if(! isOrganismValid){
    organism<-NULL
  }
  
  scdata@rowRanges@metadata$hvg <- .getMeanVarTrend(logcounts(scdata))
  ##explained by feature###
  scdata@rowRanges@elementMetadata$genevar_by_counts<-.getVarExplainedbyFeature(scdata,"log10_total_counts")
  scdata@rowRanges@elementMetadata$genevar_by_features<-.getVarExplainedbyFeature(scdata,"log10_total_features")
  scdata@rowRanges@elementMetadata$genevar_by_Mt<-.getVarExplainedbyFeature(scdata,"log10_total_counts_Mt")
  scdata@rowRanges@elementMetadata$genevar_by_rRNA<-.getVarExplainedbyFeature(scdata,"log10_total_counts_rRNA")
  
  ##select the top HVGs highly variable genes for the PCA, default is 1000
  hvggenes <-  rownames(scdata@rowRanges@metadata$hvg)[order(scdata@rowRanges@metadata$hvg$zval,decreasing=T)][1:nHVGs]
  
  tdata<-t(logcounts(scdata)[rownames(scdata)%in%hvggenes, , drop = FALSE])
  nPCs<-min(nPCs, min(dim(tdata))-1)
  
  pcaresult <- prcomp_irlba(tdata, n= nPCs)
  
  design <- model.matrix( ~ pcaresult$x[, PCind])
  fit <- lmFit(logcounts(scdata)[rownames(scdata)%in%hvggenes, , drop = FALSE], design)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  
  ##select the top 500 pc1 genes
  scdata@metadata$pc1genes <- topTable(fit, coef = 2, number = 500)
  ## enriched pathways in top 1000 hvgs and 500 pc1 genes
  ##query the webgestalt to see if the organism exist or not...
  if(!is.null(organism)){
    scdata@metadata$hvgPathway <- .getIndividualPathway(hvggenes,organism=organism)
    scdata@metadata$pc1Pathway <- .getIndividualPathway(rownames(scdata@metadata$pc1genes), organism=organism)
  }
  
  return(scdata)
} 
