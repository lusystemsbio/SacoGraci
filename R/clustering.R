#' RACIPE simulations
#' @param dftop: circuit topology
#' @param nModels: number of RACIPE models generated
#' @param integrateStepSize: step size for the ODE integration
#' @param simulationTime: simulation time of ODE for each RACIPE model
#' @return logscData: simulated data after log scaling and standardization
#' @export
#' @import sRACIPE
#' @import SummarizedExperiment
gen_RACIPE <-function(dftop, nModels, integrateStepSize = 0.02, simulationTime = 200)
{
#  require(sRACIPE)
  rac2agCore<-sRACIPE::sracipeSimulate(circuit = dftop, numModels = nModels, plots = FALSE, 
                              integrateStepSize, simulationTime)
  
  
  dataRowGene<-assay(rac2agCore,1)
  dataColGene<-t(as.matrix(dataRowGene))
  
  dataM<-as.matrix(dataColGene)
  dataLog<-log2(dataM)
  logscData<-scale(dataLog)
  return(logscData)
}

#' Model Clustering By hierarchical clustering analysis (HCA)
#' @param logscData: simulated data after log scaling and standardization
#' @param numbClust: number of model clusters
#' @return res: model clustering output: (1: cluster sizes; 2: clustering rearranged data; 3: cluster indices)
#' @export
modClustHCA<-function(logscData,numbClust)
{

  hr<-hclust(as.dist(1-cor(t(logscData), method="pearson")), method="ward.D2")
  
  myCl<-cutree(hr,k=numbClust)
  
  dataRearr<-NULL
  clSize<-rep(0,numbClust)
  for (i in 1:numbClust)
  {
    inI<-which(myCl==i)
    clSize[i]<-length(inI)
    dataRearr<-rbind(dataRearr,logscData[inI,])
  }
  

  cluster_ind = integer()
  for(i in 1: length(clSize)){
    cluster_ind = c(cluster_ind, rep(i, clSize[i]))
  }
  
  res<-list(clSize = clSize, dataRearr = dataRearr, clInd = cluster_ind)
  
  return(res)
}

#' Generating HCA Heatmap 
#' @param logscData: simulated data after log scaling and standardization
#' @export
gen_heatmap_hca <- function(logscData){
  dist.pear = function(x) as.dist(1-cor(t(x)))
  hclust.ward = function(x) hclust(x, method="ward.D2")
  
  ht3emt<-heatmap(logscData, distfun=dist.pear, hclustfun=hclust.ward,scale="none")
}

#' Generating PCA Scatterplot
#' @param logscData: simulated data after log scaling and standardization
#' @export
gen_pca_plot <- function(logscData){
  pca_results <- prcomp(logscData, center = TRUE, scale = TRUE)
  var_explained <- pca_results$sdev^2/sum(pca_results$sdev^2)
  plot(pca_results$x[,1], pca_results$x[,2], type = "p",
       xlab = paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       ylab = paste0("PC2: ",round(var_explained[2]*100,1),"%"))
}

#' MODEL CLUSTERING BY K-MEANS
#' @param data: data matrix for k-means clustering (either logscData or projected data)
#' @param numbClust: number of model clusters
#' @param clustCenters: coordinates of the cluster center
#' @return res: model clustering output
#' @export
modClustKmeans<-function(data, numbClust, clustCenters)
{
  
  clSize<-rep(0,numbClust)
  permModels<-rep(0,dim(pcData)[2])
  kCent<-kmeans(data, centers=clustCenters,iter.max=10000,nstart=1)
  for(t in 1:numbClust)
  {
    inI<-which(kCent$cluster==t)
    clSize[t]<-length(inI)
    permModels<-c(permModels,inI)
  }
  
  cluster_ind = integer()
  for(i in 1: length(clSize)){
    cluster_ind = c(cluster_ind, rep(i, clSize[i]))
  }
  res<-list(clSize = clSize, permM = permModels, clInd = cluster_ind)
  
  return(res)
}

#' GENE CLUSTERING BY INDIVIDUAL MODELS (HCA)
#' @param logscData: simulated data after log scaling and standardization
#' @param numbGeneClust: number of gene clusters
#' @return gene_list: gene clustering output
#' @export
geneClustInd<-function(logscData,numbGeneClust)
{
  hr<-hclust(as.dist(1-cor(logscData, method="pearson")), method="ward.D2")
  myCl<-cutree(hr,k=numbGeneClust)
  geneNames<-colnames(logscData)
  
  gene_list<-list()
  for (i in 1:numbGeneClust)
  {
    gene_list[[i]]<-geneNames[which(myCl==i)]
  }
  
  return(gene_list)
}

#' GENE CLUSTERING BY MEDIAN VALUES
#' @param clustData: model clustering data 
#' @param clSize: size of each model clusters (an output of the model clustering results)
#' @param numbGeneClust: number of gene clusters
#' @return gene_list: gene clustering output
#' @export
geneClustMedian<-function(clustData,clSize, numbGeneClust)
{
  
  cumsCL<-cumsum(clSize)
  medM<-apply(clustData[1:cumsCL[1],],2,median)
  
  for(i in 2:length(clSize))
  {
    mCl<-apply(clustData[(cumsCL[i-1]+1):(cumsCL[i]),],2,median)
    medM<-rbind(medM,mCl)
  }
  hr<-hclust(as.dist(1-cor(medM, method="pearson")), method="ward.D2")
  myCl<-cutree(hr,k=numbGeneClust)
  geneNames<-colnames(clustData)
  
  gene_list<-list()
  for (i in 1:numbGeneClust)
  {
    gene_list[[i]]<-geneNames[which(myCl==i)]
  }
  
  return(gene_list)
}

#' Reordering the gene clusters in the gene expression matrix
#' @param logscData: current gene expression data matrix
#' @param gene_list: gene clustering output
#' @param geneGroupOrder: desired order of gene groups
#' @return a list of reordered data (logscData) and an updated gene list (gene_list)
#' @export
reordering<-function(logscData, gene_list, geneGroupOrder = NULL) {
  shufCol<-NULL
  numbGeneClust = length(gene_list)
  numbGene = ncol(logscData)
  
  if(is.null(geneGroupOrder)) geneGroupOrder = 1:numbGeneClust
  
  for(i in 1:numbGeneClust)
  {
    shufCol<-c(shufCol,gene_list[[geneGroupOrder[i]]])
  }
  permCol<-rep(0,numbGene)
  for(t in 1:numbGene)
  {
    permCol[t]<-which(colnames(logscData)==shufCol[t])[1]
  }
  
  gene_list_new<-list()
  for(i in 1:numbGeneClust)
  {
    gene_list_new[[i]]<-gene_list[[geneGroupOrder[i]]]
  }
  
  return(list(data = logscData[,permCol], gene_list = gene_list_new))
}

#' Identify the center and radius of each model cluster
#' @param data: gene expression matrix
#' @param clusterRef: cluster indices of all models 
#' @param percThr: Threshold of permutation test 
#' @return : a list of (center, variance, radius(centroid & medoid)) for each cluster
#' @export
centMedVarCutDistPerc<-function(data, clusterRef, percThr=0.01)
{
  dataRow = t(as.matrix(data))
  
  numb_genes<-dim(dataRow)[1]
  modelsRef<-dim(dataRow)[2]
  numbClusterRef<-length(unique(clusterRef))
  
  centMedRef<-matrix(0, nrow=2*numbClusterRef, ncol=numb_genes)
  
  for (i in 1:numbClusterRef)
  {
    centMedRef[i,]<-apply(dataRow[,which(clusterRef==i)], 1, mean)
    centMedRef[i+numbClusterRef,]<-apply(dataRow[,which(clusterRef==i)],1, median)
  }
  
  sampMRef<-NULL
  for (i in 1:numb_genes)
  {
    sCGr<-sample(dataRow[i,], size=modelsRef, replace=TRUE)
    sampMRef<-cbind(sampMRef, sCGr)
  } 
  
  cutOffDistM<-matrix(0, nrow=2, ncol=numbClusterRef)
  for (j in 1:numbClusterRef)
  {
    distCRef<-dist(rbind(centMedRef[j,],sampMRef), method="euclidean")
    sortDCR<-sort(distCRef[1:modelsRef])
    pvalC05<-sortDCR[as.integer(percThr*modelsRef)]
    cutOffDistM[1,j]<-pvalC05
    
    distMRef<-dist(rbind(centMedRef[j+numbClusterRef,],sampMRef), method="euclidean")
    sortDMR<-sort(distMRef[1:modelsRef])
    pvalM05<-sortDMR[as.integer(percThr*modelsRef)]
    cutOffDistM[2,j]<-pvalM05
  }
  
  varCl<-rep(0, numbClusterRef)
  
  for (i in 1:numbClusterRef)
  {
    clsize<-length(which(clusterRef==i))
    dcref<-dist(rbind(centMedRef[i,], t(dataRow[,which(clusterRef==i)])), method="euclidean")
    dintr<-dcref[1:clsize]
    varCl[i]<-1/clsize*sum(dintr*dintr)
  }
  colnames(centMedRef)<-rownames(dataRow)
  
  return(list(center = centMedRef, variance = varCl, radius = cutOffDistM))
  
}
