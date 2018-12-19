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

.getRange<-function(values){
  rangeCount<-range(values)
  medianCount<-median(values)
  result<-paste0("[", rangeCount[1], "-", medianCount, "-", rangeCount[2], "]")
  return(result)
}

.getWebGestaltPathway <- function(genes, organism) {
  spathway <- WebGestaltR(enrichMethod = "ORA", organism = organism,
                          enrichDatabase = "pathway_KEGG", interestGene = genes,
                          interestGeneType = "genesymbol", referenceSet = "genome",
                          is.output = FALSE)
  if (is.null(spathway) | typeof(spathway) == "character") {
    return(NULL)
  } else {
    sdata <- data.frame(Pathway = gsub(" - .*", "", spathway$description),
                        FDR = -log10(spathway$FDR), stringsAsFactors = F)
    return(sdata)
  }
}

.getIndividualPathway <- function(sobj, organism) {

  sgenes <- rownames(sobj)
  sdata <- .getWebGestaltPathway(sgenes, organism)
  return(sdata)
}

.getMultiplePathway <- function(sces, metaObjectName) {
  sdata<-NULL
  for (i in 1:length(sces)) {
    fNo <- which(names(sces[[i]]) == metaObjectName)
    if(length(fNo) == 0){
      next
    }
    spathway <- sces[[i]][fNo][[1]]
    spathway$Sample <- names(sces)[i]
    sdata <- rbind(sdata, spathway)
  }
  
  if(is.null(sdata)){
    stop(paste0(metaObjectName, " is not exists in object sces"))
  }
  
  sdata$FDR[sdata$FDR == Inf] <- max(sdata$FDR[sdata$FDR != Inf]) + 1
  
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

.getColData <- function(sces, feature){
  if(missing(feature)){
    stop("Need to specify feature of .getColData")
  }
  
  fNo <- which(names(sces[[1]]) == feature)
  if(length(fNo) == 0){
    stop(paste0("Feature ", feature, " is not exists in object sces"))
  }
  
  result<-NULL
  for (i in 1:length(sces)) {
    result<-rbind(result, data.frame(Sample=names(sces)[i], Value=sces[[i]][fNo]))
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

.getDiffGenes <- function(scesall, organism, FDR = 0.01, geneNo = 50) {
  if(length(unique(scesall$condition)) == 1){
    return(list(genes = NULL, pathways = NULL))
  }
  
  design <- model.matrix( ~ 0 + as.factor(scesall$condition))
  snames <- unique(scesall$condition)
  colnames(design) <- snames
  
  cont <- c()
  compareNames <- c()
  
  for (i in 1 : (length(snames)-1)) {
    for (j in (i + 1) : length(snames)) {
      cont <- c(cont, paste0(snames[i], " - ", snames[j]))
      compareNames <- c(compareNames, paste0(snames[i], "_VS_", snames[j]))
    }
  }
  
  cat("performing differential analysis ...\n")
  
  fit <- lmFit(scesall$logcounts, design)
  contrast.matrix <- makeContrasts(contrasts = cont, levels = design)
  
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)
  
  coefNo <- length(cont)
  pairTables <- list()
  
  for (i in 1 : coefNo) {
    pairTables[[i]] <- topTable(fit2, coef = i, num = dim(scesall$logcounts)[1], sort.by = "none")
  }
  names(pairTables) <- cont
  
  diffglist <- c()
  for (i in 1:coefNo) {
    diffvs <- pairTables[[i]][abs(pairTables[[i]]$logFC) > 1 & pairTables[[i]]$adj.P.Val < FDR, ]
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
    mDiffFC <- dcast(diffFC, Gene ~ Comparison, value.var = "LogFold", fill = 0)
    rownames(mDiffFC) <- mDiffFC$Gene
    
    for (con in cont) {
      if (!(con %in% colnames(mDiffFC))) {
        mDiffFC[, con] <- rnorm(nrow(mDiffFC), 0, 0.01)
      }
    }
    mDiffFC <- mDiffFC[, -1, drop = F] 
    
    if (!missing(organism)) {
      diffPathList <- NULL
      for (i in 1:coefNo) {
        cat("pathway analysis of", i, ":", cont[i], "\n")
        
        diffvs <- pairTables[[i]][abs(pairTables[[i]]$logFC) > 1 & pairTables[[i]]$adj.P.Val < FDR, ]
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
        mDiffPathway <- dcast(diffPathList, Pathway ~ Comparison, value.var = "FDR", fill = 0)
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
### select the top 50 genes (adjustable) with FDR<0.01
### .getBiologicalSimilarity(sces, objectName="hvg", filterName="FDR", valueName="bio")
### .getBiologicalSimilarity(sces, objectName="pc1genes", filterName="adj.P.Val", valueName="logFC")
.getBiologicalSimilarity <- function(sces, objectName, filterName, valueName, defaultValue = 0) {
  objIndex <- which(names(sces[[1]]) == objectName)
  sobj <- sces[[1]][objIndex][[1]]
  filterIndex  <- which(colnames(sobj) == filterName)
  valueIndex  <- which(colnames(sobj) == valueName)
  
  genelist <- c()
  for (i in 1:length(sces)) {
    sobj <- sces[[i]][objIndex][[1]]
    sobj <- sobj[sobj[, filterIndex] < 0.01, ]
    sgene <- rownames(sobj)[order(abs(sobj[, valueName]), decreasing = TRUE)][1:min(50, dim(sobj)[1])]
    genelist <- c(genelist, sgene)
  }
  genelist <- unique(genelist)
  
  sdata <- NULL
  for (i in 1:length(sces)) {
    sobj <- sces[[i]][objIndex][[1]]
    matchid <- rownames(sobj) %in% genelist
    filtered <- sobj[matchid, ]
    sdata <- rbind(sdata, data.frame(Sample = names(sces)[i], 
                                     Feature = rownames(filtered), 
                                     Value = filtered[, valueIndex]))
  }
  
  mdata <- dcast(sdata, Feature ~ Sample, value.var = "Value", fill = defaultValue)
  rownames(mdata) <- mdata$Feature
  mdata <- as.matrix(mdata[, c(2:ncol(mdata))])
  
  return(mdata)
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

.findOutlier <- function (dat, nmads = 3, type = c("lower", "higher"), upper.limit=NA, lower.limit=NA) {
  
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
  
  return(dat < lower.limit | upper.limit < dat)
}

.getVarExplainedbyFeature <- function(sce, feature, chunk = 1000) {
  
  exprs_mat <- sce$data
 
  
  x <- sce[feature][[1]]

  
  ngenes <- nrow(sce$data)
  
  if (ngenes > chunk) {
    by.chunk <- cut(seq_len(ngenes), ceiling(ngenes/chunk))
  } else {
    by.chunk <- factor(integer(ngenes))
  }
  
  R_squared <- numeric(ngenes)
  
  for (element in levels(by.chunk)) {
    current <- by.chunk == element
    cur.exprs <- exprs_mat[current, , drop = FALSE]
   
    R_squared[current] <- cor(t(cur.exprs),x,use="pairwise.complete.obs") ^ 2 
  }
  
  
  
  Pct_Var_Explained <- 100 * R_squared
  return(Pct_Var_Explained)
}



##find the topn highly variable genes and the mean-variance trend
##return the hvginfo, including the mean, variance and zval of top 500 gene
## return the mean-variance trend
.getMeanVarTrend<-function(sce,chunk=1000,topn=500){
    
	ngenes <- nrow(sce$data)
  
  if (ngenes > chunk) {
    by.chunk <- cut(seq_len(ngenes), ceiling(ngenes/chunk))
  } else {
    by.chunk <- factor(integer(ngenes))
  }
  
  gene_mean <- gene_var<- numeric(ngenes)
  
  
  for (element in levels(by.chunk)) {
    current <- by.chunk == element
    cur.exprs <- sce$data[current, , drop = FALSE]
   
    gene_mean[current] <- apply(cur.exprs,1,mean)
	gene_var[current]<-apply(cur.exprs,1,var)
	
  }
 
  data_x_bin<-cut(x=gene_mean,breaks=150)
  mean_x<-tapply(X=gene_mean,INDEX=data_x_bin,FUN=mean)
  mean_y <- tapply(X = gene_var, INDEX = data_x_bin, FUN=mean)
  sd_y<- tapply(X=gene_var,INDEX=data_x_bin,FUN=sd)
  gene.dispersion.scaled <- (gene_var - mean_y[as.numeric(x = data_x_bin)])/sd_y[as.numeric(x = data_x_bin)]
  gene.dispersion.scaled[is.na(x = gene.dispersion.scaled)] <- 0
  
  
  hvginfo<-data.frame(mean=gene_mean,var=gene_var,zval=gene.dispersion.scaled)
  rownames(hvginfo) <- rownames (sce$data)
  hvginfo<-hvginfo[order(hvginfo$zval,decreasing=T),]
  
  data_x_bin_plot<-cut(x=mean_x,breaks=20)
  mean_x_plot<-tapply(X=mean_x,INDEX=data_x_bin_plot, FUN=mean)
  mean_y_plot<-tapply(X=mean_y,INDEX=data_x_bin_plot,FUN=mean)
  
  return(list(hvginfo=hvginfo[1:topn,],trend=data.frame(x=mean_x_plot,y=mean_y_plot)))
 }



.prepareTableSummary <- function(sces) {
  pw <- matrix(nrow = length(sces), ncol = 13)
  
  for (i in 1:length(sces)) {
    pw[i, 1] <- names(sces)[i]
    pw[i, 2] <- sum(sces[[i]]$rawdata) # Count
    pw[i, 3] <- dim(sces[[i]]$rawdata)[2] #Cell
    pw[i, 4] <- dim(sces[[i]]$rawdata)[1] #Gene
    pw[i, 5] <- paste0("[", summary(Matrix::colSums(sces[[i]]$rawdata))[1],
                       "-", summary(Matrix::colSums(sces[[i]]$rawdata))[3],
                       "-", summary(Matrix::colSums(sces[[i]]$rawdata))[6], "]") # R-Count
    
    pw[i, 6] <- paste0("[", summary(Matrix::colSums(sces[[i]]$rawdata != 0))[1],
                       "-", summary(Matrix::colSums(sces[[i]]$rawdata != 0))[3],
                       "-", summary(Matrix::colSums(sces[[i]]$rawdata != 0))[6], "]") # R-Gene
    
    pw[i, 7] <- paste0(format(as.numeric(as.character(max(sces[[i]]$pct_counts_Mt))), digits = 2, nsmall = 1), "%") #mtRNA
    
    pw[i, 8] <- paste0(format(as.numeric(as.character(max(sces[[i]]$pct_counts_rRNA))), digits = 2, nsmall = 1), "%") # rRNA
    pw[i, 9] <- sum(sces[[i]]$libsize.drop) # F-Count
    pw[i, 10] <- sum(sces[[i]]$feature.drop) # F-Gene
    pw[i, 11] <- 0 # F-rRNA
    pw[i, 12] <- sum(sces[[i]]$mito.drop) # F-mtRNA
    pw[i, 13] <- sum(as.numeric(pw[i, 9:12]))
  }
  
  colnames(pw) <- c("EID", "Count", "Cell", "Gene", "R-Count", "R-Gene", 
                    "mtRNA", "rRNA", "F-Count", "F-Gene", "F-rRNA", "F-mtRNA", "F")
  
  return(as.data.frame(pw))
}
