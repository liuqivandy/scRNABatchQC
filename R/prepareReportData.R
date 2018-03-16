prepareReportData<-function(sampleTable, organism){
  sces <- prepareSCRNADataSet(sampleTable, organism)
  scesall <- preparePCATSNEData(sces)
  diffFC <- .getDiffGenes(scesall, organism = organism,  FDR = 0.01, geneNo = 50)
  hvgPathways<-.getMultiplePathway(sces, metaObjectName="hvgPathway")
  pc1Pathways<-.getMultiplePathway(sces, metaObjectName="pc1Pathway")
  plotData<-list(sces=sces, scesall=scesall, diffFC=diffFC, hvgPathways=hvgPathways, pc1Pathways=pc1Pathways)
  return(plotData)
}
