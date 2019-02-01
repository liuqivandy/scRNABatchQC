#' prepareSCRNADataSet


#' @description  Generate statistics information for multiple single-cell RNA-seq datasets represented by gene-count matrices
#' @param inputfiles; a character vector giving the input file name (or a URL starting http://, file://, etc.) of gene-by-cell count matrices, the rowname should be gene symbol; each file should be regular delimited file;  Compressed files ending .gz and .bz2 are supported
#' @param sampleRatio float; the ratio of cells sampled from each experiment to examine the expression similarity(default: 1, use the full dataset without sampling)
#' @param nHVGs integer; the number of highly variable genes (default: 1000)
#' @param sf integer; Scale factor  (default: 10000)
 
#' @return a named list containing each scRNA data sets, 
#'                                
#'                                

#' @examples
#' library(scRNABatchQC)
#' sces<-prepareSCRNADataSet(inputfiles=c("data1.csv","data2.csv"))

#' @import R.utils ggplot2 gplots limma data.table irlba Rtsne WebGestaltR rmdformats Matrix statmod DT RCurl
#' @importFrom devtools session_info
#' 
#' @export

prepareSCRNADataSet <- function(inputfiles, samplenames=NULL,organism=c("hsapiens","mmusculus"), sampleRatio=1,nHVGs=1000, sf=10000){
    
  organism<-match.arg(organism)
 
  result <- list()
  nfiles<-length(inputfiles)
  
  if (is.null(samplenames)){samplenames=paste0("S",1:nfiles)}
  if (nfiles!=length(samplenames)) {stop("the inputfiles and samplenames donnot have the same length", call. = FALSE)}
  if (sum(samplenames!=make.names(samplenames))>0) samplenames<-paste0("S",samplenames)
 
  for (ind in 1:nfiles) {
   
    
    cat("Preparing ", samplenames[ind], "\n")
    
    result[[ind]] <- prepareSCRNAData(inputfiles[ind], organism=organism, sampleRatio=sampleRatio,nHVGs=nHVGs, sf=sf)
  }
  
  names(result) <- samplenames
  return(result)
}

#' preparePCATSNEData
#'
#' @description generate the PCA and tSNE data for the combined samples.
#'
#' @param sces a named list of scRNA data
#' @param nHVGs integer; the number of highly variable genes (default:1000)
#' @param nPC integer: the number of principal components (default:10)
#' @return a list with PCA and TSNE data
#' @importFrom Rtsne Rtsne
#' @importFrom Matrix Matrix
#' @export preparePCATSNEData
#' @examples 
#' library(scRNABatchQC)
#' sces <- prepareSCRNADataSet(inputfiles=c("data1.csv","data2.csv"))
#' sceall <- preparePCATSNEData(sces)
preparePCATSNEData <- function(sces, nHVGs=1000, nPC= 10) {

   pca_tsne_data <- list()
   
   lenind<-sapply(sces,function(x){dim(x$data)[2]})
    pca_tsne_data$condition <-rep(names(lenind),lenind)


    pca_tsne_data$logcounts<-sces[[1]]$data

     colnames(pca_tsne_data$logcounts) <- paste0(names(lenind)[1], 1:lenind[1])

    if(length(sces) > 1){

	  for (i in 2:length(sces)) {

	  	mat<-sces[[i]]$data
                colnames(mat)<-paste0(names(lenind)[i], 1:lenind[i])
	  	pca_tsne_data$logcounts <- .mergeSparseMatrix(pca_tsne_data$logcounts, mat)
	  }

     }

  

  

 

   pca_tsne_data$hvg <- .getMeanVarTrend(pca_tsne_data$logcounts)

  

  ##select the top 1000 highly variable genes for the PCA

  feature_set <-  rownames(pca_tsne_data$hvg)[order(pca_tsne_data$hvg$zval,decreasing=T)][1:nHVGs]

  #scevar <- apply(scesdata, 1, var)

  

  #feature_set <- head(order(scevar, decreasing = T), n = 500)

   

 pca_tsne_data$pca <- prcomp_irlba(t(pca_tsne_data$logcounts[rownames(pca_tsne_data$logcounts)%in%feature_set, , drop = FALSE]), n = nPC)

  
  set.seed(100)

  tsne_out <- Rtsne(pca_tsne_data$pca$x, initial_dims = ncol(pca_tsne_data$pca$x), pca = FALSE, perplexity = 20, check_duplicates = FALSE)
  
  pca_tsne_data$tsne <- tsne_out$Y


  return(pca_tsne_data)

}



