###EX5 DATA (EMT network)

library(SacoGraci)
library(sRACIPE)
set.seed(42)

# Load circuit topology
top_emt = read.csv("data/ex3_net.csv", header=T)

# RACIPE simulation
racEMT = gen_RACIPE(top_emt, 500)

# Clustering
gen_heatmap_hca(racEMT)

###Eliminate input nodes/genes
badN<-c("miR9","miR30c","miR205")
badI<-NULL
for(i in 1:length(badN))
{
  indb<-which(colnames(racEMT)==badN[i])[1]
  badI<-c(badI,indb) 
}
racEMT=racEMT[,-badI]
pca=prcomp(racEMT, scale=FALSE)

plot(pca$x[,1], pca$x[,2], xlab="PC1", ylab="PC2")


pcEMT=cbind(pca$x[,1],pca$x[,2])

centM=c(-4,0)
centM=rbind(centM,c(0.5,1.5))
centM=rbind(centM,c(1.5,-1.5))
centM=rbind(centM,c(4,0))

resK<-modClustKmeans_cc(pcEMT,numbClust=4, clustCenters=centM)
data_reordered<-racEMT[resK$permM,]
myclGenes<-geneClustMedian(data_reordered,clSize=resK$clSize,numbGeneClust=4)

shG = c(3,2,1,4)# my designed order
processed_results = recordering(data_reordered, shG, myclGenes)

# permutation to compute the radii of clusters
myd01_statCl<-centMedVarCutDistPerc(processed_results$data, resK$clInd, 0.01)

# generate initial starting CG circuits
inTopsM <-gaInitial_gen(top_emt, processed_results$gene_list, 90)

# end of data processing, ready for circuit optimization (see Tutorial)

