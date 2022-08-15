###EX5 DATA (EMT network)

library(SacoGraci)
set.seed(42)

# Load circuit topology
top_ex3 = read.csv("ex3_net.csv", header=T)

# RACIPE simulation
### use this for testing:  racEx3 = gen_RACIPE(top_ex3, 100)
racEx3 = gen_RACIPE(dftop = top_ex3, nModels = 10000)

# Clustering
gen_heatmap_hca(racEx3)

###Eliminate input nodes/genes
badNi = c("miR9","miR30c","miR205")
badI = NULL
for(i in 1:length(badN)){
  indb = which(colnames(racEx3)==badN[i])[1]
  badI = c(badI,indb) 
}
racEx3_new = racEx3[,-badI]
pca = prcomp(racEx3_new, scale=FALSE)

plot(pca$x[,1], pca$x[,2], xlab="PC1", ylab="PC2")

pcEMT=cbind(pca$x[,1],pca$x[,2])

centM=c(-4,0)
centM=rbind(centM,c(0.5,1.5))
centM=rbind(centM,c(1.5,-1.5))
centM=rbind(centM,c(4,0))

resK<-modClustKmeans(pcEMT,numbClust=4, clustCenters=centM)
data_reordered<-racEx3_new[resK$permM,]
myclGenes<-geneClustMedian(data_reordered,clSize=resK$clSize,numbGeneClust=4)

shG = c(3,2,1,4)# my designed order
processed_results = reordering(data_reordered, myclGenes, shG)
 
# permutation to compute the radii of clusters
myd01_statCl<-centMedVarCutDistPerc(processed_results$data, resK$clInd, 0.01)

# generate initial starting CG circuits
inTopsM <-gaInitial_gen(top_ex3, processed_results$gene_list, 90)

# end of data processing, ready for circuit optimization (see Tutorial)
