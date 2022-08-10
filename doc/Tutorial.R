## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE, fig.width = 7, fig.height = 6, fig.align = "center")

## -----------------------------------------------------------------------------
library(SacoGraci)
set.seed(43)
top_ex = read.csv("Ex_topology.csv", header=T)

## -----------------------------------------------------------------------------
data_exp = readRDS("exp_data_RACIPE_ex.RDS")

## -----------------------------------------------------------------------------
gen_heatmap_hca(logscData = data_exp)  # heatmap visualization
gen_pca_plot(logscData = data_exp)  # pca visualization

output_model_grouping = modClustHCA(logscData = data_exp,numbClust = 3)  # model clustering with HCA, we specify 3 clusters
data_reordered = output_model_grouping$dataRearr  # reordered gene expression data by the clustering outcome
clusterRef = output_model_grouping$clInd  # cluster indices of all models

output_gene_grouping = geneClustInd(logscData = data_reordered,numbGeneClust = 4)  # gene clustering with HCA, we specify 4 cluster

shG = c(2,3,4,1)  # an example of the designed order to arrange genes. (optional; if not provided, the default sequence will be used) 
output_processed = reordering(logscData = data_reordered, gene_list = output_gene_grouping, geneGroupOrder = shG) # reordered gene expression data by the clustering outcome

data_processed = output_processed$data
gene_list = output_processed$gene_list

## -----------------------------------------------------------------------------
stat_clusters <-centMedVarCutDistPerc(data = data_processed, clusterRef = clusterRef, percThr = 0.01)

## -----------------------------------------------------------------------------
inTopsM <-gaInitial_gen(circuit_top = top_ex, gene_list = gene_list, numbNewTop = 90)

## ---- message=FALSE-----------------------------------------------------------
circuit1 = opt_MH(network_top = top_ex, data = data_processed, clusterRef = clusterRef, 
                  cenMedRef = stat_clusters$center, cutOffM = stat_clusters$radius, 
                  gene_list = gene_list, init_top = inTopsM[1,], 
                  output = "Results1", nRepeat= 1, nIter = 10, modelsCGr = 100, 
                  tempM = 60)

## -----------------------------------------------------------------------------
sessionInfo()

