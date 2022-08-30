#' Circuit scoring 
#' @param dataRow gene expression matrix with gene in rows
#' @param clusterRef the cluster indices of all models
#' @param cenMedRef cluster centers
#' @param cutOffM cluster radii
#' @param gene_list gene clustering output 
#' @param inNodes a list of input nodes
#' @param topol_cgr the current CG circuit topology
#' @param modelsCGr the number of RACIPE models to be simulated (10000)
#' @return scoresOuts: 1: the total score, 2: the number of noisy models
#' @export
#' @import sRACIPE
#' @import SummarizedExperiment
#' @importFrom stats median dist
inNWsimilarityRefCGr<-function(dataRow, clusterRef, cenMedRef, cutOffM, gene_list, inNodes, topol_cgr, modelsCGr = 10000)
{
#  require(sRACIPE)
  #topol_cgr<-newt_top
  #dataRow doesn't include expression(noise) of input genes 
  #topol_cgr<-lastTop
  #inNodes<-0
  numb_genes<-dim(dataRow)[1]#=number of non-input species(genes)
  modelsRef<-dim(dataRow)[2]
  numbCGRn<-length(gene_list)#it includes input and output nodes
  numbClusterRef<-length(unique(clusterRef))#it doesn't include input nodes
  scoresOut<-rep(0,4)
  
  df_lastTop<-convAdjTop(topol_cgr, numbCGRn, 1:numbCGRn)
  
  #SIMULATE CGR-CIRCUIT
  #set.seed(42)
  racCGR<-suppressMessages(sracipeSimulate(circuit = df_lastTop, numModels = modelsCGr, plots = FALSE, 
                          integrateStepSize = 0.02, simulationTime = 100))
  
  
  dataRowGenePert<-assay(racCGR,1)#This is our  data 
  dataColGenePert<-t(as.matrix(dataRowGenePert))
  logP<-log2(dataColGenePert)
  scLogP<-scale(logP)
  dataRowGeneCGR<-t(scLogP)#This is our coarse-grained data 
  
  if (setequal(0,inNodes)==FALSE) #If I have input species
  {
    cgrNamer<-rownames(dataRowGeneCGR)
    indInNodes<-NULL
    for(t in 1:length(inNodes))
    {
      ndr<-which(cgrNamer==as.character(inNodes[t]))[1]
      indInNodes<-c(indInNodes,ndr)
    }
    dataRowGeneCGR<-dataRowGeneCGR[-indInNodes,]
  }
  
  
  dataColGeneCGR<-t(as.matrix(dataRowGeneCGR))
  
  ###EXTEND THE COARSE-GRAINED DATA TO NUMBER OF GENES DATA
  ###CALCULATE THE CENTROIDS AND MEDOIDS OF EACH CLUSTER OF REF DATA
  ###AND CUTOFF DISTANCES WRT ITSELF AND CGR DATA 
  
  cgrNames<-rownames(dataRowGeneCGR)
  newRownames<-NULL
  dataSimBig<-NULL #matrix(0, nrow=numb_genes, ncol=modelsCGr)
  rfill<-0
  for(ila in 1:length(cgrNames))
  {
    tla<-as.numeric(cgrNames[ila])
    mla<-length(gene_list[[tla]])
    for(lla in (rfill+1):(rfill+mla))
    {
      dataSimBig<-rbind(dataSimBig,dataRowGeneCGR[ila,])#dataSimBig[lla,]<-dataRowGeneCGR[ila,]
      newRownames<-c(newRownames, gene_list[[tla]][lla-rfill])
    }
    rfill<-rfill+mla  
  }
  dataSimBig<-as.matrix(dataSimBig)
  rownames(dataSimBig)<-newRownames
  
  
  #dataSimBig is our simulated data
  
  
  colnames(dataSimBig) <- seq_len(ncol(dataSimBig))
  
  permGenes<-rep(0, rfill)
  for(ila in 1:rfill)
  {
    permGenes[ila]<-which(newRownames==rownames(dataRow)[ila])[1]
  }
  
  shuffledDataSim<-dataSimBig[permGenes,]
  
  ###shuffledDataSim will be compared to "experimental"/reference data
  
  
  clusterCutSimM<-matrix(0, nrow=modelsCGr, ncol=3)
  
  
  for (tla in 1:modelsCGr)
  {
    
    dmedMod<-dist(rbind(shuffledDataSim[,tla],cenMedRef[(numbClusterRef+1):(2*numbClusterRef),]), method="euclidean")
    dIntm<-dmedMod[1:numbClusterRef]
    ratDm<-dIntm/cutOffM[2,]
    indClm<-which(ratDm==min(ratDm))[1]
    if(ratDm[indClm]<1)
    {
      clusterCutSimM[tla,1]<-indClm
      
    }
    
    clusterCutSimM[tla,2]<-indClm
    
    if(ratDm[indClm]>1)
    {
      clusterCutSimM[tla,1]<-0
    }
    
    clusterCutSimM[tla,3]<-ratDm[indClm]
    
  }
  
  indZe<-which(clusterCutSimM[,1]==0)
  scoresOut[1]<-length(indZe)
  
  sizeCl<-rep(0,numbClusterRef)
  sizeClref<-rep(0,numbClusterRef)
  avgCl<-rep(0,numbClusterRef)
  chScore<-FALSE
  penMisCl<-FALSE

  for(ila in 1:numbClusterRef)
  {
    indCAv<-which(clusterCutSimM[,1]==ila)
    sizeCl[ila]<-length(indCAv)
    sizeClref[ila]<-length(which(clusterRef==ila))
    if(sizeCl[ila]>0)
    {
     avgCl[ila]<-median(clusterCutSimM[indCAv,3])
    }
    if(sizeCl[ila]==0)
    {
     avgCl[ila]<-10
    }

    if(sizeCl[ila]<(0.5*sizeClref[ila]))
    {
      chScore<-TRUE
    }
    
  }
  
  scoresOut[2]<-(modelsCGr/numbClusterRef)*sum(avgCl)
  if(scoresOut[1]==0)
  {
    scoresOut[3]<-0
  }
  if(scoresOut[1]>0)
  {
    scoresOut[3]<-sum(clusterCutSimM[indZe,3])
  }
  scoresOut[4]<-scoresOut[2]+scoresOut[3]
  if(chScore==TRUE)
  {
    scoresOut[4]<-scoresOut[4]+50000
  }
  
  return(c(scoresOut[4],scoresOut[1]))
} 

