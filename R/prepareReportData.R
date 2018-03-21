prepareReportData<-function(sampleTable, organism, cachePrefix){
  if(!missing(cachePrefix)){
    scesFile = paste0(cachePrefix, "_sces.rdata")
    if(file.exists(scesFile)){
      load(scesFile)
    }
    else{
      sces <- prepareSCRNADataSet(sampleTable, organism)
      save(sces, file=scesFile)
      
    }
  }else{
    sces <- prepareSCRNADataSet(sampleTable, organism)
  }
  
  if(!missing(cachePrefix)){
    scesAllFile = paste0(cachePrefix, "_scesall.rdata")
    if(file.exists(scesAllFile)){
      load(scesAllFile)
    }
    else{
      scesall <- preparePCATSNEData(sces)
      save(scesall, file=scesAllFile)
      
    }
  }else{
    scesall <- preparePCATSNEData(sces)
  }

  if(!missing(cachePrefix)){
    diffFCFile = paste0(cachePrefix, "_diffFC.rdata")
    if(file.exists(diffFCFile)){
      load(diffFCFile)
    }
    else{
      diffFC <- .getDiffGenes(scesall, organism = organism,  FDR = 0.01, geneNo = 50)
      save(diffFC, file=diffFCFile)
    }
  }else{
    diffFC <- .getDiffGenes(scesall, organism = organism,  FDR = 0.01, geneNo = 50)
  }
  
  hvgPathways<-.getMultiplePathway(sces, metaObjectName="hvgPathway")
  pc1Pathways<-.getMultiplePathway(sces, metaObjectName="pc1Pathway")
  plotData<-list(sces=sces, scesall=scesall, diffFC=diffFC, hvgPathways=hvgPathways, pc1Pathways=pc1Pathways)
  return(plotData)
}
