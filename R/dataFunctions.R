
.getWebGestaltPathway <- function(genes, organism) {
  spathway<-WebGestaltR(enrichMethod="ORA",organism=organism,
                        enrichDatabase="pathway_KEGG",interestGene=genes,
                        interestGeneType="genesymbol",referenceSet="genome",
                        is.output=FALSE)
  sdata<-data.frame(Pathway=gsub(" - .*", "", spathway$description),
                    FDR=-log10(spathway$FDR),
                    stringsAsFactors = F)
  return(sdata)
}

.getIndividualPathway <- function(sobj, filterName, organism) {
  filterIndex  <- which(colnames(sobj) == filterName)
  
  sobj<-sobj[sobj[, filterIndex] < 0.01, ]
  sgenes <- rownames(sobj)
  sdata<-.getWebGestaltPathway(sgenes, organism)
  return(sdata)
}

.getMultiplePathway <- function(sces, metaObjectName) {
  fNo <- which(names(metadata(sces[[1]]$sce)) == metaObjectName)
  if(fNo == 0){
    stop(paste0(metaObjectName, " is not exists in object sces"))
  }
  
  sdata<-NULL
  for (i in 1:length(sces)) {
    sce<-sces[[i]]$sce
    spathway = metadata(sce)[[fNo]]
    spathway$Sample = names(sces)[i]
    sdata<-rbind(sdata, spathway)
  }
  
  sdata$FDR[sdata$FDR==Inf]<-max(sdata$FDR[sdata$FDR!=Inf]) + 1
  
  mdata<-reshape2::dcast(sdata, Pathway~Sample, value.var="FDR", fill=0)
  rownames(mdata)<-mdata$Pathway
  mdata<-as.matrix(mdata[,c(2:ncol(mdata))])

  return(mdata)
}

.getColData<-function(sces, feature){
  if(missing(feature)){
    stop("Need to specify feature of .getColData")
  }
  
  fNo <- which(colnames(colData(sces[[1]]$sce)) == feature)
  if(fNo == 0){
    stop(paste0("Feature ", feature, " is not exists in object sces"))
  }
  
  result<-NULL
  for (i in 1:length(sces)) {
    result<-rbind(result, data.frame(Sample=names(sces)[i], Value=colData(sces[[i]]$sce)[, fNo]))
  }
  return(result)
}

.getCbindRowData<-function(sces, feature){
  if(missing(feature)){
    stop("Need to specify feature of .getRowData")
  }
  
  fNo <- which(colnames(rowData(sces[[1]]$sce)) == feature)
  if(fNo == 0){
    stop(paste0("Feature ", feature, " is not exists in object sces"))
  }
  
  result<-NULL
  for (i in 1:length(sces)) {
    result<-cbind(result, rowData(sces[[i]]$sce)[, fNo])
  }
  colnames(result)<-names(sces)
  return(result)
}

.getDiffGenes <- function(sceall, organism, FDR = 0.01, geneNo = 50) {
  design <- model.matrix( ~ 0 + as.factor(sceall$condition))
  snames <- unique(sceall$condition)
  colnames(design) <- snames
  
  cont <- c()
  compareNames <- c()
  
  for (i in 1 : (length(snames)-1)) {
    for (j in (i + 1) : length(snames)) {
      cont <- c(cont, paste0(snames[i], " - ", snames[j]))
      compareNames <- c(compareNames, paste0(snames[i], "_VS_", snames[j]))
    }
  }
  
  fit <- lmFit(logcounts(sceall), design)
  contrast.matrix <- makeContrasts(contrasts = cont, levels = design)
  
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)
  
  coefNo <- length(cont)
  pairTables <- list()
  
  for (i in 1 : coefNo) {
    pairTables[[i]] <- topTable(fit2, coef = i, num = dim(sceall)[1], sort.by = "none")
  }
  names(pairTables) <- cont
  
  diffPathList<-NULL
  diffglist <- c()
  for (i in 1:coefNo) {
    diffvs <- pairTables[[i]][abs(pairTables[[i]]$logFC) > 1 & pairTables[[i]]$adj.P.Val < FDR, ]
    diffgenes <- rownames(diffvs)[order(abs(diffvs$logFC), decreasing = TRUE)][1:min(geneNo, dim(diffvs)[1])]
    diffglist <- unique(c(diffglist, diffgenes))
    
    if(!missing(organism)){
      pathList<-.getWebGestaltPathway(rownames(diffvs), organism)
      pathList$Comparison <- cont[[i]]
      diffPathList<-rbind(diffPathList, pathList)
    }
  }
  
  diffFC<-NULL
  for (i in 1:coefNo) {
    matchid <- rownames(pairTables[[i]]) %in% diffglist
    diffFC <- rbind(diffFC, data.frame(Comparison=cont[i], Gene=rownames(pairTables[[i]])[matchid], AbsLogFold=abs(pairTables[[i]]$logFC[matchid])))
  }
  mDiffFC<-dcast(diffFC, Gene ~ Comparison, value.var="AbsLogFold", fill=0)
  rownames(mDiffFC)<-mDiffFC$Gene
  mDiffFC<-mDiffFC[,-1] 
  
  mDiffPathway<-dcast(diffPathList, Pathway ~ Comparison, value.var="FDR", fill=0)
  rownames(mDiffFC)<-mDiffFC$Gene
  mDiffFC<-mDiffFC[,-1] 
  
  
  return(mDiffFC)
}
