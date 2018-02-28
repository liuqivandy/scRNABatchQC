
getHighVarGenes <- function(sces, FDR = 0.01, geneNo = 50) {
  hvgenes <- list()
  
  hvglist <- c()
  for (i in 1:length(sces)) {
    hvout <- sces[[i]]$hvg[sces[[i]]$hvg$FDR < FDR, ]
    hvgenes[[i]] <- rownames(hvout)[order(hvout$bio, decreasing = TRUE)][1:min(geneNo, dim(hvout)[1])]
    hvglist <- unique(c(hvglist, hvgenes[[i]]))
  }
  names(hvgenes) <- names(sces)
  
  matchid <- rownames(sces[[1]]$sce) %in% hvglist
  
  hvgbio <- sces[[1]]$hvg$bio[matchid]
  for (i in 2:length(sces)) {
    hvgbio <- cbind(hvgbio, sces[[i]]$hvg$bio[matchid])
  }
  
  rownames(hvgbio) <- rownames(sces[[1]]$sce)[matchid]
  colnames(hvgbio) <- names(sces)
  return(hvgbio)
}

getTopPC1genes <- function(sces, FDR = 0.01, geneNo = 50) {
  pcgenes <- list()
  
  pcglist <- c()
  for (i in 1:length(sces)) {
    pcout <- sces[[i]]$pc1genes[sces[[i]]$pc1genes$adj.P.Val < FDR, ]
    pcgenes[[i]] <- rownames(pcout)[order(abs(pcout$logFC), decreasing = TRUE)][1:min(geneNo, dim(pcout)[1])]
    pcglist <- unique(c(pcglist, pcgenes[[i]]))
  }
  names(pcgenes) <- names(sces)
  
  matchid <- rownames(sces[[1]]$sce) %in% pcglist
  
  pcgFC <- sces[[1]]$pc1genes$logFC[matchid]
  for (i in 2:length(sces)) {
    pcgFC <- cbind(pcgFC, sces[[i]]$pc1genes$logFC[matchid])
  }
  
  rownames(pcgFC) <- rownames(sces[[1]]$sce)[matchid]
  colnames(pcgFC) <- names(sces)
  return(pcgFC)
}

getDiffGenes <- function(sceall, FDR = 0.01, geneNo = 50) {
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
  
  diffgenes <- list()
  diffglist <- c()
  
  for (i in 1:coefNo) {
    diffvs <- pairTables[[i]][abs(pairTables[[i]]$logFC) > 1 & pairTables[[i]]$adj.P.Val < FDR, ]
    diffgenes[[i]] <- rownames(diffvs)[order(abs(diffvs$logFC), decreasing = TRUE)][1:min(geneNo, dim(diffvs)[1])]
    diffglist <- unique(c(diffglist, diffgenes[[i]]))
  }
  names(diffgenes) <- compareNames
  
  matchid <- rownames(pairTables[[1]]) %in% diffglist
  
  diffFC <- pairTables[[1]]$logFC[matchid]
  for (i in 2:coefNo) {
    diffFC <- cbind(diffFC, pairTables[[i]]$logFC[matchid])
  }
  rownames(diffFC) <- rownames(pairTables[[1]])[matchid]
  colnames(diffFC) <- compareNames
  return(diffFC)
}

