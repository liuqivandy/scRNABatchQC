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
prepareSCRNADataSet <- function(sampleTable, organism){
  result <- list()
  
  n <- 1
  for (n in 1:nrow(sampleTable)) {
    sampleName<-as.character(sampleTable[n,1])
    countFile<-as.character(sampleTable[n,2])
    cat("Preparing ", sampleName, "\n")
    counts <- as.matrix(read.csv(countFile, row.names=1, header=T))
    result[[n]] <- prepareSCRNAData(counts, organism)
  }
  names(result) <- sampleTable[, 1]
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
  
  allct <- as(sces[[1]]$data, "RsparseMatrix")
  conditions <- rep(names(sces)[1], dim(sces[[1]]$data)[2])
  colnames(allct) <- paste0(names(sces)[1], "cell", 1:dim(sces[[1]]$data)[2])

  if(length(sces) > 1){
	  for (i in 2:length(sces)) {
	  	mat <- as(sces[[i]]$data, "RsparseMatrix")
	  	colnames(mat) <- paste0(names(sces)[i], "cell", 1:dim(sces[[i]]$data)[2])

	  	allct <- .mergeSparseMatrix(allct, mat)
	  	conditions <- c(conditions, rep(names(sces)[i], dim(sces[[i]]$data)[2]))
	  }
  }
  
  lib_size <- Matrix::colSums(allct)/mean(Matrix::colSums(allct))
  
  counts_norm_lib_size <- t(apply(allct, 1, function(x) x/lib_size ))
  num.cells <- Matrix::rowSums(allct != 0)
  to.keep <- num.cells > 0
  
  scesdata <- log2(counts_norm_lib_size[to.keep, ] + 1)
  pca_tsne_data$logcounts <- scesdata
  
  scevar <- apply(scesdata, 1, var)
  
  feature_set <- head(order(scevar, decreasing = T), n = 500)
  pca_tsne_data$pca <- prcomp(t(scesdata[feature_set, , drop = FALSE]), rank. = ncomponents)
  
  pca_tsne_data$condition <- sapply(rownames(pca_tsne_data$pca$x), function(x) strsplit(x, "cell")[[1]][1])

  set.seed(100)
  tsne_out <- Rtsne(pca_tsne_data$pca$x, initial_dims = ncol(pca_tsne_data$pca$x), pca = FALSE, perplexity = perplexity, check_duplicates = FALSE)
  
  pca_tsne_data$tsne <- tsne_out$Y

  return(pca_tsne_data)
}

