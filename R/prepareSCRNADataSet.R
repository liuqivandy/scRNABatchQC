##' prepareSCRNADataSet
##'
##' The function read multiple count table and prepare statistics information for each table
##'
##' @param sampleTable the sample table with first column as sample name, second column as file location
##' @param organism the organism for annotation
##' @return a named list of scRNA data
##' @export prepareSCRNADataSet
##' @examples 
##' #sampleTable <- data.frame(Sample = c("S1", "S2", "S3"), 
##'                            File = c("count1.csv", "count2.csv", "count3.csv"))
##' #sces <- prepareSCRNADataSet(sampleTable)
prepareSCRNADataSet <- function(inputfiles, samplenames=NULL,organism){
  result <- list()
  nfiles<-length(inputfiles)
  
  if (is.null(samplenames)){samplenames=paste0("S",1:nfiles)}
  if (nfiles!=length(samplenames)) {stop("the inputfiles and samplenames donnot have the same length", call. = FALSE)}
  if (sum(samplenames!=make.names(samplenames))>0) samplenames<-paste0("S",samplenames)
 
  for (ind in 1:nfiles) {
   
    
    cat("Preparing ", samplenames[ind], "\n")
    
    result[[ind]] <- prepareSCRNAData(inputfiles[ind], organism)
  }
  
  names(result) <- samplenames
  return(result)
}

##' preparePCATSNEData
##'
##' The function prepare statistics information from multiple scRNA dataset.
##'
##' @param sces a named list of scRNA data
##' @return a list with PCA and TSNE data
##' @importFrom Rtsne Rtsne
##' @importFrom Matrix Matrix
##' @export preparePCATSNEData
##' @examples 
##' #sces <- prepareSCRNADataSet(sampleTable)
##' #sceall <- preparePCATSNEData(sces)
preparePCATSNEData <- function(sces, ncomponents = 10, perplexity = 20) {

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

  feature_set <-  rownames(pca_tsne_data$hvg)[order(pca_tsne_data$hvg$zval,decreasing=T),][1:1000]

  #scevar <- apply(scesdata, 1, var)

  

  #feature_set <- head(order(scevar, decreasing = T), n = 500)

   

 pca_tsne_data$pca <- prcomp_irlba(t(pca_tsne_data$logcounts[rownames(pca_tsne_data$logcounts)%in%feature_set, , drop = FALSE]), n = ncomponents)

  
  set.seed(100)

  tsne_out <- Rtsne(pca_tsne_data$pca$x, initial_dims = ncol(pca_tsne_data$pca$x), pca = FALSE, perplexity = perplexity, check_duplicates = FALSE)
  
  pca_tsne_data$tsne <- tsne_out$Y


  return(pca_tsne_data)

}



