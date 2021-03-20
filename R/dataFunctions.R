fread_bychunk<-function(inputfile,chunk.size=2000) {
  
  done=FALSE
  chunk = 1
  tmp = fread(inputfile,skip=(chunk-1)*chunk.size,nrow=chunk.size,data.table=F)
  data_sparse<-.tosparse(tmp[,-1])
  rownames(data_sparse)<-tmp[,1]
  if (nrow(tmp)<chunk.size) done=TRUE
  chunk=chunk+1
  while(!done)
  {
    tmp = fread(inputfile,skip=(chunk-1)*chunk.size+1,nrow=chunk.size,data.table=F)
    ind<-nrow(data_sparse)
    data_sparse<-rbind(data_sparse,.tosparse(tmp[,-1]))
    rownames(data_sparse)[(ind+1):nrow(data_sparse)]<-tmp[,1]
    chunk = chunk + 1
    if(nrow(tmp)<chunk.size) done = TRUE
    rm(tmp)
    gc()
  }
 
  return(data_sparse)
  
}

.is_10X_v2<-function(inputFolder){
  barcodeFile = paste0(inputFolder, "/barcodes.tsv")	
  featureFile = paste0(inputFolder, "/genes.tsv")	
  matrixFile = paste0(inputFolder, "/matrix.mtx")	
  
  if(!file.exists(barcodeFile)){
    return(FALSE)
  }
  if(!file.exists(featureFile)){
    return(FALSE)
  }
  if(!file.exists(matrixFile)){
    return(FALSE)
  }
  return(TRUE)
}


.is_10X_v3<-function(inputFolder){
  barcodeFile = paste0(inputFolder, "/barcodes.tsv.gz")	
  featureFile = paste0(inputFolder, "/features.tsv.gz")	
  matrixFile = paste0(inputFolder, "/matrix.mtx.gz")	
  
  if(!file.exists(barcodeFile)){
    return(FALSE)
  }
  if(!file.exists(featureFile)){
    return(FALSE)
  }
  if(!file.exists(matrixFile)){
    return(FALSE)
  }
  return(TRUE)
}

read_10X_v2<-function(inputFolder){
  barcodeFile = paste0(inputFolder, "/barcodes.tsv")	
  featureFile = paste0(inputFolder, "/genes.tsv")	
  matrixFile = paste0(inputFolder, "/matrix.mtx")	
  
  if(!file.exists(barcodeFile)){
    stop(paste0("File not found ", barcodeFile))
  }
  if(!file.exists(featureFile)){
    stop(paste0("File not found ", featureFile))
  }
  if(!file.exists(matrixFile)){
    stop(paste0("File not found ", matrixFile))
  }
  
  features = read.delim(featureFile, header=F, stringsAsFactors = F)
  barcodes = readLines(barcodeFile)
  mat = readMM(matrixFile)
  
  colnames(mat)<-barcodes
  rownames(mat)<-features$V1
  result <-aggregate.Matrix(mat, features$V2, FUN=sum, na.rm=TRUE)
  return(result)
}

read_10X_v3<-function(inputFolder){
  barcodeFile = paste0(inputFolder, "/barcodes.tsv.gz")	
  featureFile = paste0(inputFolder, "/features.tsv.gz")	
  matrixFile = paste0(inputFolder, "/matrix.mtx.gz")	
  
  if(!file.exists(barcodeFile)){
    stop(paste0("File not found ", barcodeFile))
  }
  if(!file.exists(featureFile)){
    stop(paste0("File not found ", featureFile))
  }
  if(!file.exists(matrixFile)){
    stop(paste0("File not found ", matrixFile))
  }
  
  features = read.delim(featureFile, header=F, stringsAsFactors = F)
  barcodes = readLines(barcodeFile)
  mat = readMM(matrixFile)
  
  colnames(mat)<-barcodes
  rownames(mat)<-features$V1
  result <-aggregate.Matrix(mat, features$V2, FUN=sum, na.rm=TRUE)
  return(result)
}

.detach_package <- function(pkg, character.only = FALSE)
{
  if(!character.only)
  {
    pkg <- deparse(substitute(pkg))
  }
  search_item <- paste("package", pkg, sep = ":")
  while(search_item %in% search())
  {
    detach(search_item, unload = TRUE, character.only = TRUE)
  }
  unloadNamespace(pkg)
}

.isOrganismValid<-function(organism, showWarning=TRUE){
  if(is.null(organism)){
    return(FALSE)
  }
  
  validOrganisms<-listOrganism()
  result<-organism %in% validOrganisms
  if((! result) & showWarning){
    warning(paste0("organism ", organism, " is not supported by WebGestalt, no pathway analysis will be performed. Valid organisms include: ", paste(validOrganisms,collapse=", ")))
  }
  return(result)
}

##convert the data frame to the sparse Matrix
.tosparse<-function(datatable){
  i_list <- lapply(datatable, function(x) which(x != 0))
  counts <- unlist(lapply(i_list, length), use.names = F)
  sparseMatrix(
    i = unlist(i_list, use.names = F),
    j = rep(1:ncol(datatable), counts),
    x = unlist(lapply(datatable, function(x) x[x != 0]), use.names = F),
    dims = dim(datatable),
    dimnames = dimnames(datatable))
}

#######	
.getRange<-function(values){
  rangeCount<-range(values)
  medianCount<-median(values)
  result<-paste0("[", rangeCount[1], "-", medianCount, "-", rangeCount[2], "]")
  return(result)
}

####
.getWebGestaltPathway <- function(genes, organism) {
  tryres=try(spathway<-WebGestaltR(enrichMethod = "ORA", organism = organism,
                                   enrichDatabase = "pathway_KEGG", interestGene = genes,
                                   interestGeneType = "genesymbol", referenceSet = "genome",
                                   isOutput = FALSE))
  
  if(class(tryres) == "try-error"){
    return(NULL)
  }
  
  if (is.null(spathway) | typeof(spathway) == "character") {
    return(NULL)
  } 
  
  sdata <- data.frame(Pathway = gsub(" - .*", "", spathway$description),
                      FDR = -log10(spathway$FDR), stringsAsFactors = F)
  return(sdata)
}


###
.getIndividualPathway <- function(sgenes, organism) {
  sdata <- .getWebGestaltPathway(sgenes, organism)
  return(sdata)
}

###
.getMultiplePathways <- function(sces, metaObjectName) {
  sdata<-NULL
  for (i in 1:length(sces)) {
    fNo <- which(names(sces[[i]]@metadata) == metaObjectName)
    if(length(fNo) == 0){
      next
    }
    spathway <- sces[[i]]@metadata[fNo][[1]]
    spathway$Sample <- names(sces)[i]
    sdata <- rbind(sdata, spathway)
  }
  
  if(is.null(sdata)){
    warning(paste0(metaObjectName, " is not exists in object sces"))
    return(NULL)
  }
  
  sdata$FDR[sdata$FDR == Inf] <- max(sdata$FDR[sdata$FDR != Inf]) + 1
  
  if(!("Pathway" %in% colnames(sdata))){
    sdata$Pathway<-sdata$description
  }
  
  mdata <- reshape2::dcast(sdata, Pathway ~ Sample, value.var = "FDR", fill = 0)
  
  for(sample in names(sces)){
    if(!(sample %in% colnames(mdata))){
      mdata[,sample] <- abs(rnorm(nrow(mdata), 0, 0.01))
    }
  }
  
  rownames(mdata) <- mdata$Pathway
  mdata <- as.matrix(mdata[, c(2 : ncol(mdata)), drop = F])
  
  return(mdata)
}

.getRawColData <- function(sces, feature){
  if(missing(feature)){
    stop("Need to specify feature of .getColData")
  }
  if (is.null(names(sces))) names(sces)<-1:length(sces)
  
  fNo <- which(names(sces[[1]]@metadata$rawmeta$CellData) == feature)
  if(length(fNo) == 0){
    stop(paste0("Feature ", feature, " is not exists in object sces"))
  }
  
  result<-NULL
  for (i in 1:length(sces)) {
    result<-rbind(result, data.frame(Sample=names(sces)[i], Value=sces[[i]]@metadata$rawmeta$CellData[fNo]))
  }
  colnames(result) <- c("Sample", "Value")
  return(result)
}

.getGeneData <- function(sces, feature){
  if(missing(feature)){
    stop("Need to specify feature of .getColData")
  }
  if (is.null(names(sces))) names(sces)<-1:length(sces)
  fNo <- which(names(sces[[1]]@rowRanges@elementMetadata) == feature)
  if(length(fNo) == 0){
    stop(paste0("Feature ", feature, " is not exists in object sces"))
  }
  
  result<-NULL
  for (i in 1:length(sces)) {
    result<-rbind(result, data.frame(Sample=names(sces)[i], Value=sces[[i]]@rowRanges@elementMetadata[fNo]))
  }
  colnames(result) <- c("Sample", "Value")
  return(result)
}


.getCbindRowData<-function(sces, feature){
  if(missing(feature)){
    stop("Need to specify feature of .getRowData")
  }
  
  fNo <- which(colnames(rowData(sces[[1]]$sce)) == feature)
  if(length(fNo) == 0){
    stop(paste0("Feature ", feature, " is not exists in object sces"))
  }
  
  result<-NULL
  for (i in 1:length(sces)) {
    result<-cbind(result, rowData(sces[[i]]$sce)[, fNo])
  }
  colnames(result)<-names(sces)
  return(result)
}
###############

##################################
.getDiffGenes <- function(scesall, organism, logFC=1,FDR = 0.01, geneNo = 50,chunk=1000) {
  
  
  
  if(length(unique(scesall@colData$condition)) == 1){
    return(list(genes = NULL, pathways = NULL))
  }
  design <- model.matrix( ~ 0 + as.factor(scesall@colData$condition))
  snames <- unique(scesall@colData$condition)
  colnames(design) <- snames
  cont <- c()
  compareNames <- c()
  
  for (i in 1 : (length(snames)-1)) {
    for (j in (i + 1) : length(snames)) {
      cont <- c(cont, paste0(snames[i], " - ", snames[j]))
      compareNames <- c(compareNames, paste0(snames[i], "_VS_", snames[j]))
    }
  }
  
  filtereddata<-logcounts(scesall)[scesall@rowRanges@metadata$hvg$var>0,]
  
  cat("performing differential analysis ...\n")
  ngenes <- nrow(filtereddata)
  if (ngenes > chunk) {
    by.chunk <- cut(seq_len(ngenes), ceiling(ngenes/chunk))
  } else {
    by.chunk <- factor(integer(ngenes))
  }
  coefNo <- length(cont)
  pairTables <- vector("list",coefNo)
  
  for (element in levels(by.chunk)) {
    current <- by.chunk == element
    cur.exprs <- filtereddata[current, , drop = FALSE]
    fit <- lmFit(cur.exprs, design)
    contrast.matrix <- makeContrasts(contrasts = cont, levels = design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)
    for (i in 1 : coefNo) {
      pairTables[[i]] <- rbind(pairTables[[i]],topTable(fit2, coef = i, number = dim(cur.exprs)[1], sort.by = "none"))
    }
  }
  
  names(pairTables) <- cont
  diffglist <- c()
  for (i in 1:coefNo) {
    diffvs <- pairTables[[i]][abs(pairTables[[i]]$logFC) > logFC & pairTables[[i]]$adj.P.Val < FDR, ]
    diffgenes <- rownames(diffvs)[order(abs(diffvs$logFC), decreasing = TRUE)][1:min(geneNo, dim(diffvs)[1])]
    diffglist <- unique(c(diffglist, diffgenes))
  }
  
  diffglist <- na.omit(diffglist)
  
  mDiffFC <- NULL
  mDiffPathway <- NULL
  if(length(diffglist) > 0){
    diffFC <- NULL
    for (i in 1:coefNo) {
      matchid <- rownames(pairTables[[i]]) %in% diffglist
      diffFC <- rbind(diffFC, data.frame(Comparison = cont[i], 
                                         Gene = rownames(pairTables[[i]])[matchid], 
                                         LogFold = pairTables[[i]]$logFC[matchid]))
    }
    mDiffFC <- reshape2::dcast(diffFC, Gene ~ Comparison, value.var = "LogFold", fill = 0)
    rownames(mDiffFC) <- mDiffFC$Gene
    
    for (con in cont) {
      if (!(con %in% colnames(mDiffFC))) {
        mDiffFC[, con] <- rnorm(nrow(mDiffFC), 0, 0.01)
      }
    }
    mDiffFC <- mDiffFC[, -1, drop = F] 
    
    if (!is.null(organism)) {
      diffPathList <- NULL
      for (i in 1:coefNo) {
        cat("pathway analysis of", i, ":", cont[i], "\n")
        
        diffvs <- pairTables[[i]][abs(pairTables[[i]]$logFC) > logFC & pairTables[[i]]$adj.P.Val < FDR, ]
        if (nrow(diffvs) > 1) {
          alldiffgenes <- rownames(diffvs)
          pathList <- .getWebGestaltPathway(alldiffgenes, organism)
          if (!is.null(pathList)) {
            pathList$Comparison <- cont[[i]]
            diffPathList <- rbind(diffPathList, pathList)
          }
        }
      }
      
      if (!is.null(diffPathList)) {
        infDiffIndex <- diffPathList$FDR == Inf
        if (sum(infDiffIndex) > 0) {
          maxFdr <- max(diffPathList$FDR[diffPathList$FDR != Inf])
          diffPathList$FDR[infDiffIndex] <- maxFdr + 1
        }
        mDiffPathway <- reshape2::dcast(diffPathList, Pathway ~ Comparison, value.var = "FDR", fill = 0)
        for (con in cont){
          if (!(con %in% colnames(mDiffPathway))) {
            mDiffPathway[, con] <- abs(rnorm(nrow(mDiffPathway), 0, 0.01))
          }
        }
        rownames(mDiffPathway) <- mDiffPathway$Pathway
        mDiffPathway <- mDiffPathway[, -1, drop = F]
      }
    }
  }
  
  r <- list(genes = mDiffFC, pathways = mDiffPathway)
  return(r)
}

### Biological features similarity
### select the top 50 genes (adjustable) to plot the heatmap and pathway analysis
### .getBiologicalSimilarity(sces, objectName="hvg",valueName="zval")
### .getBiologicalSimilarity(sces, objectName="pc1genes",valueName="logFC")
.getBiologicalSimilarity <- function(sces, objectName="hvg",valueName="zval",topn=50,defaultValue=0) {
  #filterIndex  <- which(colnames(sobj) == filterName)
  genelist <- c()
  for (i in 1:length(sces)) {
    
    if (objectName=="pc1genes") sobj=sces[[i]]@metadata$pc1genes
    #sobj <- sobj[sobj[, filterIndex] < 0.01, ]
    if (objectName=="hvg") {
      sobj=sces[[i]]@rowRanges@metadata$hvg
      sobj<-sobj[order(sobj[valueName],decreasing=T),]}
  }
  sgene <- rownames(head(sobj,topn))
  genelist <- c(genelist, sgene)
  genelist <- unique(genelist)
  
  sdata <- NULL
  for (i in 1:length(sces)) {
    
    
    if (objectName=="pc1genes") sobj=sces[[i]]@metadata$pc1genes
    #sobj <- sobj[sobj[, filterIndex] < 0.01, ]
    if (objectName=="hvg")   sobj=sces[[i]]@rowRanges@metadata$hvg
    
    matchid <- rownames(sobj) %in% genelist
    filtered <- sobj[matchid, ]
    Value=filtered[valueName]
    if (objectName=="pc1genes") Value=abs(Value )
    sdata <- rbind(sdata, data.frame(Sample = names(sces)[i], 
                                     Feature = rownames(filtered), 
                                     Value = Value))
  }
  
  mdata <- reshape2::dcast(sdata, Feature ~ Sample, value.var = valueName, fill = defaultValue)
  
  mdatax <- as.matrix(mdata[, c(2:ncol(mdata))])
  rownames(mdatax) <- mdata$Feature
  return(mdatax)
}

.mergeSparseMatrix <- function(mat1, mat2) {
  mat1_diff_mat2_rownames <- setdiff(rownames(mat2), rownames(mat1))
  mat2_diff_mat1_rownames <- setdiff(rownames(mat1), rownames(mat2))
  
  allrown <- union(rownames(mat1), rownames(mat2))
  
  suppmat1 <- Matrix(nrow = length(mat1_diff_mat2_rownames), ncol = ncol(mat1), 0)
  rownames(suppmat1) <- mat1_diff_mat2_rownames
  colnames(suppmat1) <- colnames(mat1)
  
  suppmat2 <- Matrix(nrow = length(mat2_diff_mat1_rownames), ncol = ncol(mat2), 0)
  rownames(suppmat2) <- mat2_diff_mat1_rownames
  colnames(suppmat2) <- colnames(mat2)
  
  mat1_more <- rbind(mat1, suppmat1)
  mat2_more <- rbind(mat2, suppmat2)
  
  mat <- cbind(mat1_more[allrown, ], mat2_more[allrown, ])
  
  return(mat)
}

.findOutlier <- function (dat, nmads = 3, log=FALSE, type = c("lower", "higher"), upper.limit=NA, lower.limit=NA) {
  
  if ((!is.na(lower.limit) & lower.limit<0) | (!is.na(upper.limit) & upper.limit<0)) stop("the lower or upper cutoff should be NA or >0")
  if(log){
    dat<-log10(dat); 
    if (!is.na(upper.limit)) upper.limit=log10(upper.limit); 
    if (!is.na(lower.limit)) lower.limit=log10(lower.limit)
  }
  
  med <- median(dat, na.rm = TRUE)
  mad <- mad(dat, center = med, na.rm = TRUE)
  
  diff.val <- nmads * mad
  upper.limit <- min(upper.limit, med + diff.val,na.rm=TRUE)
  lower.limit <- max(lower.limit, med - diff.val,na.rm=TRUE)
  
  type <- match.arg(type)
  if (type == "lower") {
    upper.limit <- Inf
  } else if (type == "higher") {
    lower.limit <- -Inf
  }
  result<-dat < lower.limit | upper.limit < dat
  cutoff<-ifelse (type=="lower",lower.limit,upper.limit)
  if (log) {cutoff=10^cutoff}
  return(list(filtered=result,cutoff=cutoff))
}

.getVarExplainedbyFeature <- function(sce, feature, chunk = 1000) {
  exprs_mat <- logcounts(sce)
  
  feature.data <- sce@colData[[feature]]
  expf<-mean(feature.data)
  sdf<-sd(feature.data)
  
  ngenes <- nrow(exprs_mat)
  ncells<-ncol(exprs_mat)
  
  if (ngenes > chunk) {
    by.chunk <- cut(seq_len(ngenes), ceiling(ngenes/chunk))
  } else {
    by.chunk <- factor(integer(ngenes))
  }
  
  R_squared <- numeric(ngenes)
  
  for (element in levels(by.chunk)) {
    current <- by.chunk == element
    cur.exprs <- exprs_mat[current, , drop = FALSE]
    Exy= cur.exprs%*%feature.data/ncells-sce@rowRanges@metadata$hvg$mean[current]*expf
    sdxy=sdf*sqrt(sce@rowRanges@metadata$hvg$var[current])
    R_squared[current] <- (Exy/sdxy) ^ 2 
  }
  
  Pct_Var_Explained <- 100 * R_squared
  return(Pct_Var_Explained)
  
} 

##find the highly variable genes and the mean-variance trend
##return the hvginfo, including the mean, variance and zval of each gene
## return the mean-variance trend
.getMeanVarTrend<-function(data,chunk=1000){
  ngenes <- nrow(data)
  ncol<-ncol(data)
  
  if (ngenes > chunk) {
    by.chunk <- cut(seq_len(ngenes), ceiling(ngenes/chunk))
  } else {
    by.chunk <- factor(integer(ngenes))
  }
  
  gene_mean <- gene_var<- numeric(ngenes)
  for (element in levels(by.chunk)) {
    current <- by.chunk == element
    cur.exprs <- data[current, , drop = FALSE]
    gene_mean[current]<-tmpMean <- Matrix::rowMeans(cur.exprs)
    gene_var[current]<- Matrix::rowSums((cur.exprs-tmpMean)^2)/(ncol-1)
  }
  
  data_x_bin<-cut(x=gene_mean,breaks=150)
  mean_x<-tapply(X=gene_mean,INDEX=data_x_bin,FUN=mean)
  mean_y <- tapply(X = gene_var, INDEX = data_x_bin, FUN=mean)
  sd_y<- tapply(X=gene_var,INDEX=data_x_bin,FUN=sd)
  gene.dispersion.scaled <- (gene_var - mean_y[as.numeric(x = data_x_bin)])/sd_y[as.numeric(x = data_x_bin)]
  gene.dispersion.scaled[is.na(x = gene.dispersion.scaled)] <- 0
  hvginfo<-data.frame(mean=gene_mean,var=gene_var,zval=gene.dispersion.scaled)
  rownames(hvginfo) <- rownames (data)
  #data_x_bin_plot<-cut(x=mean_x,breaks=20)
  # mean_x_plot<-tapply(X=mean_x,INDEX=data_x_bin_plot, FUN=mean)
  # mean_y_plot<-tapply(X=mean_y,INDEX=data_x_bin_plot,FUN=mean)
  #hvginfo<-hvginfo[order(hvginfo$zval,decreasing=T),]
  
  return(hvginfo=hvginfo)
}

.prepareTableSummary <- function(sces) {
  pw <- matrix(nrow = length(sces), ncol = 17)
  
  for (i in 1:length(sces)) {
    pw[i, 1] <- names(sces)[i]
    pw[i, 2] <- sum(sces[[i]]@metadata$rawmeta$CellData$total_counts) # Count
    pw[i, 3] <- sces[[i]]@metadata$rawmeta$ncells #Cell
    pw[i, 4] <- sces[[i]]@metadata$rawmeta$ngenes #Gene
    pw[i, 5] <- paste0("[", summary(sces[[i]]@metadata$rawmeta$CellData$total_counts)[1],
                       "-", summary(sces[[i]]@metadata$rawmeta$CellData$total_counts)[3],
                       "-", summary(sces[[i]]@metadata$rawmeta$CellData$total_counts)[6], "]") # R-Count
    
    pw[i, 6] <- paste0("[", summary(sces[[i]]@metadata$rawmeta$CellData$total_features)[1],
                       "-", summary(sces[[i]]@metadata$rawmeta$CellData$total_features)[3],
                       "-", summary(sces[[i]]@metadata$rawmeta$CellData$total_features)[6], "]") # R-Gene
    
    pw[i, 7] <- paste0(format(as.numeric(as.character(max(sces[[i]]@metadata$rawmeta$CellData$pct_counts_Mt,na.rm=T)*100)), digits = 2, nsmall = 1), "%") #mtRNA
    
    pw[i, 8] <- paste0(format(as.numeric(as.character(max(sces[[i]]@metadata$rawmeta$CellData$pct_counts_rRNA,na.rm=T)*100)), digits = 2, nsmall = 1), "%") # rRNA
    pw[i, 9] <- sum(sces[[i]]@metadata$rawmeta$CellData$libsize.drop,na.rm=T) # F-Count
    pw[i, 10] <- as.integer(sces[[i]]@metadata$rawmeta$Cutoff$count) #cutoff by counts
    pw[i, 11] <- sum(sces[[i]]@metadata$rawmeta$CellData$feature.drop,na.rm=T) # F-Gene
    pw[i, 12] <- as.integer(sces[[i]]@metadata$rawmeta$Cutoff$gene)   #cutoff by gene
    pw[i, 13] <- sum(sces[[i]]@metadata$rawmeta$CellData$mito.drop,na.rm=T) # F-mtRNA
    pw[i, 14] <- paste0(round(sces[[i]]@metadata$rawmeta$Cutoff$mito,2)*100,"%")   #cutoff by mtRNA percentage
    pw[i, 15] <- sum(sces[[i]]@metadata$rawmeta$CellData$is.drop,na.rm=T)  #F-total
    pw[i, 16] <- sces[[i]]@metadata$rawmeta$ncells - as.integer(pw[i, 9]) #valid cell
    pw[i, 17] <- sces[[i]]@metadata$valid  #Valid
  }
  
  colnames(pw) <- c("EID", "Count", "Cell", "Gene", "R-Count", "R-Gene", 
                    "mtRNA", "rRNA", "F-Count", "C-Count","F-Gene","C-Gene", "F-mt","C-mt", "F", "V-Count", "Valid")
  
  return(as.data.frame(pw))
}
