###EX5 DATA (GSD network)

library(SacoGraci)
library(sRACIPE)
set.seed(43)
# Load circuit topology
top_ex5 = read.csv("data/ex5_net.csv", header=T)

# RACIPE simulation
racEx5 = gen_RACIPE(top_ex5, 500)

# Clustering
gen_heatmap_hca(racEx5)

my2clag = modClustHCA(racEx5,numbClust = 6)
data_reordered=my2clag$dataRearr

myclGenes<-geneClustMedian(data_reordered,clSize=my2clag$clSize,numbGeneClust=5)

shG = c(4,5,3,1,2)# my designed order
processed_results = recordering(data_reordered, shG, myclGenes)

# permutation to compute the radii of clusters
myd01_statCl<-centMedVarCutDistPerc(processed_results$data, my2clag$clInd, 0.01)

# generate initial starting CG circuits
inTopsM <-gaInitial_gen(top_ex5, processed_results$gene_list, 90)

# end of data processing, ready for circuit optimization (see Tutorial)
