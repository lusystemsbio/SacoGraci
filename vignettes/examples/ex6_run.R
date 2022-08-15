###EX6 DATA (OVAL420 network)

library(SacoGraci)
set.seed(42)

# Load circuit topology
top_ex6 = read.csv("ex6_net.csv", header=T)

# RACIPE simulation
### use this for testing: racEx6 = gen_RACIPE(top_ex6, 100)
racEx6 = gen_RACIPE(dftop = top_ex6, nModels = 10000)

# Clustering
gen_heatmap_hca(racEx6)

my2clag = modClustHCA(racEx6,numbClust = 3)
data_reordered = my2clag$dataRearr
my2clag410Genes = geneClustInd(data_reordered,numbGeneClust = 4)

shG = c(2,3,1,4)  # my designed order
processed_results = reordering(data_reordered, my2clag410Genes, shG)

# permutation to compute the radii of clusters
myd01_statCl = centMedVarCutDistPerc(processed_results$data, my2clag$clInd, 0.01)

# generate initial starting CG circuits
inTopsM = gaInitial_gen(top_ex6, processed_results$gene_list, 10)

# end of data processing, ready for circuit optimization (see Tutorial). In this case, only MH is needed
