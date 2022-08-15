###EX4 DATA (SCLC network)

library(SacoGraci)
set.seed(43)

top_ex4 = read.csv("ex4_net.csv", header=T)
### use this for testing: racEx4 = gen_RACIPE(top_ex4, 100)
racEx4 = gen_RACIPE(dftop = top_ex4, nModels = 10000)

gen_heatmap_hca(logscData = racEx4)  # heatmap visualization
gen_pca_plot(logscData = racEx4)  # pca visualization

output_model_grouping = modClustHCA(logscData = racEx4,numbClust = 3)  # model clustering with HCA, we specify 3 clusters
data_reordered = output_model_grouping$dataRearr  # reordered gene expression data by the clustering outcome
clusterRef = output_model_grouping$clInd  # cluster indices of all models

output_gene_grouping = geneClustInd(logscData = data_reordered,numbGeneClust = 4)  # gene clustering with HCA, we specify 4 cluster

shG = c(2,3,4,1)  # an example of the designed order to arrange genes. (optional; if not provided, the default sequence will be used) 
output_processed = reordering(logscData = data_reordered, gene_list = output_gene_grouping, geneGroupOrder = shG) # reordered gene expression data by the clustering outcome

data_processed = output_processed$data
gene_list = output_processed$gene_list

stat_clusters = centMedVarCutDistPerc(data = data_processed, clusterRef = clusterRef, percThr = 0.01)

inTopsM = gaInitial_gen(circuit_top = top_ex4, gene_list = gene_list, numbNewTop = 90)

# example code for circuit optimization
circuit1 = opt_MH(network_top = top_ex4, data = data_processed, clusterRef = clusterRef, 
                  cenMedRef = stat_clusters$center, cutOffM = stat_clusters$radius, 
                  gene_list = gene_list, init_top = inTopsM[1,], 
                  output = "Results1", nRepeat= 1, nIter = 10, modelsCGr = 100, 
                  tempM = 60)

