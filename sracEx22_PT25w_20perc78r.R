library(sRACIPE)
library(igraph)
library(doParallel)

convTopAdj<-function(data_t, nGenes, colN=1:nGenes)
{
 adjM<-rep(0, nGenes*nGenes)
 nEd<-dim(data_t)[1]
 data_tn<-matrix(0,nrow=nEd, ncol=3)
 if(is.list(data_t[,1]))
 {
  
  for(i in 1:3)
   {
    data_tn[,i]<-unlist(data_t[,i])
   }
 
 }
 else
 {

   for(i in 1:3)
   {
    data_tn[,i]<-data_t[,i]
   }

 }

 for (edg in 1:nEd)
  {#edg<-1
    vsou<-which(colN==data_tn[edg,1])[1]
    vtgt<-which(colN==data_tn[edg,2])[1]
    vtype<-data_tn[edg,3]
    adjM[(vsou-1)*nGenes+vtgt]<-vtype
  }

 return(adjM)
}



convAdjTop<-function(adjMatL, numbG, colN=1:numbG)
{
#colN=vector of gene names=1:numbG by default
adjMatLin<-unlist(adjMatL)
nnEntr<-which(adjMatLin!=0)
numb_conn_cgr<-length(nnEntr)

cgr_top<-NULL

for (i in 1:numb_conn_cgr)
{
  posN<-nnEntr[i]
  colposin<-posN%%(numbG)
  ###Stopped HERE
  if(colposin==0)
   {
    colpos<-numbG
    ropos<-(posN%/%(numbG))
   }

  if(colposin!=0)
   {
    colpos<-colposin
    ropos<-(posN%/%(numbG))+1
   }

  cgr_top<-rbind(cgr_top, c(colN[ropos], colN[colpos], adjMatLin[posN]))  

}

colnames(cgr_top)<-c("Source", "Target", "Type")

data_cgr_top<-as.data.frame(cgr_top)

return(data_cgr_top)

}





find_in_list<-function(xli, n, elem)
{
indli<-0
for ( i in 1:n)
{
if (is.element(elem, xli[[i]]))
   {indli<-i}
}
return(indli)
}


bottom_up_gene_grouping<-function(corM, n_cgnodes)
{
#returns the grouping of the genes in n_cgnodes groups 
#given their correlation matrix corM
#corM is supposed to have colnames, rownames non-empty

copycorGenes<-corM
copycorGenes[copycorGenes==1]<-0

#rn<-colnames(copycorGenes)
coln<-colnames(copycorGenes)

numb_genes<-dim(copycorGenes)[1]

list_bup<-list()
for (i in 1:numb_genes)
{
list_bup[[i]]<-coln[i]
}


while(length(list_bup)>n_cgnodes)
{

   #Find maximum correlation in the corr matrix and its position

   maxcor<-max(copycorGenes)
   maxpos<-which(copycorGenes==maxcor)[1]
   rowpos1<-maxpos%%numb_genes
   if(rowpos1==0)
   {
    rowpos<-numb_genes
    colpos<-(maxpos%/%numb_genes)
   }

   if(rowpos1!=0)
   {
    rowpos<-rowpos1
    colpos<-(maxpos%/%numb_genes)+1
   }
   g1<-coln[rowpos]
   g2<-coln[colpos]

   #Merge the two gene groups in one of them and delete the other

   posg1<-find_in_list(list_bup, length(list_bup), g1)
   posg2<-find_in_list(list_bup, length(list_bup), g2)
   if(posg1 !=posg2)
   {list_bup[[posg1]]<-c(list_bup[[posg1]], list_bup[[posg2]])
    list_bup[[posg1]]<-unique(list_bup[[posg1]])
    list_bup[[posg2]]<-NULL
   }
   copycorGenes[rowpos,colpos]<--5
   copycorGenes[colpos,rowpos]<--5


}

return(list_bup)

}




###Builds the adjacency matrix with probabilities

adjMatCGrProb<-function(data_top, gene_list)
{

#Find the connections in the course-grained network
n_cgnodes<-length(gene_list)

adj_cgr<-list()
m<-1
matTop<-as.matrix(data_top)

for( i in 1:n_cgnodes)
{

 tmem1<-is.element(matTop[,1], gene_list[[i]])
 out_ind<-which(tmem1)
 numb_out<-length(gene_list[[i]]) 
 for (j in 1:n_cgnodes)
 {
     tmem2<-is.element(matTop[,2], gene_list[[j]]) 
     in_ind<-which(tmem2)
     numb_in<-length(gene_list[[j]])

     edgij<-intersect(out_ind, in_ind)
     if (length(edgij)>0)
     { 
        
        act_len<-length(which(matTop[edgij,3]==1))
        inh_len<-length(which(matTop[edgij,3]==2))
        adj_cgr[[m]]<-c(act_len/(numb_in*numb_out), inh_len/(numb_in*numb_out), 
        1-act_len/(numb_in*numb_out)-inh_len/(numb_in*numb_out))
     }  
     if (length(edgij)==0)
     {
       adj_cgr[[m]]<-c(0,0,1)
     }
     m<-m+1
  }

}
return(adj_cgr)

}





initialize_topology<-function(adj_matrProb, numb_cgnodes)
{

#returns the adj. matrix of the most dense topology as a vector(byrow)
adj_Matr<-rep(0, numb_cgnodes*numb_cgnodes)
for (i in 1:numb_cgnodes)
{
  for (j in 1:numb_cgnodes)
  {
    probS<-adj_matrProb[[(i-1)*numb_cgnodes+j]]
    prob1<-probS[1]
    prob2<-probS[2]
    if ((prob1+prob2)>0)
    {
      if (prob1>=prob2)
      {
         adj_Matr[(i-1)*numb_cgnodes+j]<-1
      }
      else
      {
         adj_Matr[(i-1)*numb_cgnodes+j]<-2
      }
    }
    else
     {
         adj_Matr[(i-1)*numb_cgnodes+j]<-0
     }    
  }
}

df_adjMat<-convAdjTop(adj_Matr, numb_cgnodes, 1:numb_cgnodes)
grAdj<-graph_from_data_frame(df_adjMat, directed = FALSE, vertices = NULL)
vconn<-vertex_connectivity(grAdj, checks = TRUE)
if(vconn>0)
{
   return(adj_Matr)
}
else 
{ 
   return(rep(0, numb_cgnodes*numb_cgnodes))
}

}




find_topology<-function(top_list, samp_top, minN, maxN)
{
flagF<-FALSE
i<-maxN
while ( (i>=minN)&& (flagF==FALSE))
{
 if (all(samp_top==top_list[[i]]))
 {
   flagF<-TRUE
 }
 else
 {
   i<-i-1
 }
}

if(i<minN)
{
   return(0)
}
else
{
  return(i)
}

}





##################
####NEW FUNCTIONS FOR DEALING WITH INPUT/OUTPUT NODES/GENES

#adjMatCGrProb<-function(data_top, gene_list)

find_inout<-function(adjMatrPr, numb_cgnodes)
{
#adjMatrPr=adjacency matrix with probabilities of CG-circuits 
#with number of nodes=numb_cgnodes

#numb_cgnodes<-length(gene_list)
noedgeV<-c(0,0,1)

input_nodes<-NULL
for (i in 1:numb_cgnodes)
{
  t<-0
  for (j in 1:numb_cgnodes)
  {
    if (all(noedgeV==adjMatrPr[[(j-1)*numb_cgnodes+i]]))
    {
       t<-t+1
    }
  }
  
  if(t==numb_cgnodes)
  {
     input_nodes<-c(input_nodes,i)
  } 
} 


output_nodes<-NULL
for (i in 1:numb_cgnodes)
{
  t<-0
  for (j in 1:numb_cgnodes)
  {
    if (all(noedgeV==adjMatrPr[[(i-1)*numb_cgnodes+j]]))
    {
       t<-t+1
    }
  }
  
  if(t==numb_cgnodes)
  {
     output_nodes<-c(output_nodes,i)
  } 
} 

if(length(input_nodes)>0)
{
 inNodes<-input_nodes
}
if(length(input_nodes)==0)
{
 inNodes<-0
}


if(length(output_nodes)>0)
{
 outNodes<-output_nodes
}
if(length(output_nodes)==0)
{
 outNodes<-0
}

resNod<-list()
resNod[[1]]<-inNodes
resNod[[2]]<-outNodes

return(resNod)

}




inNrepsample_topology<-function(adj_matrPr, numb_cgnodes, inNodes, old_top)
{

#resInOut<-find_inout(adj_matrPr, numb_cgnodes)
#inNodes<-resInOut[[1]]
#outNodes<-resInOut[[2]]

allEd<-NULL
for (t in 1:length(adj_matrPr))
{
  allEd<-c(allEd, adj_matrPr[[t]])   
}
   
indEd0<-which(allEd>0)
indEd1<-which(allEd<1)
indEd<-intersect(indEd0, indEd1)
numb_posEdges<-length(indEd)

probChEd<-allEd[indEd]
connGraph<-FALSE
backTr<-1

while((connGraph==FALSE)&&(backTr==1))
{

new_top<-old_top

ntrials<-1

####DO RANDOM SAMPLING
while((connGraph==FALSE)&& (ntrials<=10))
{
  combEd<-sample(indEd, size=1, prob=probChEd)
  s<-combEd%%3
  if( s!=0)
   {
     chEdge<-(combEd%/%3)+1
   }
  if(s==0)
   {
     chEdge<-(combEd%/%3)
   }
  

  if (s==old_top[chEdge])#force a change
   {
    if (s==0)
     {
       indS<-3
     }
    if(s!=0)
     {
       indS<-s
     }
    
    probS<-adj_matrPr[[chEdge]]
  
    if(sum(probS[-indS])>0)
     {
       newV<-c(1,2,0)
       newV<-newV[-indS]
       s<-sample(newV, size=1, prob=probS[-indS])
     }
        
   }
  
  
  if (s!=old_top[chEdge])
  {
    new_top[chEdge]<-s
    adj_Matr<-matrix(new_top, nrow=numb_cgnodes, ncol=numb_cgnodes, byrow=TRUE)    
    
    agrInN<-FALSE
    sumInc<-apply(adj_Matr,2,sum)
   
    indInN<-which(sumInc==0)
    if(length(indInN)==0)
    {
      indInN<-0
    }
    
    agrInN<-setequal(inNodes, indInN)
    
    
    agrInOut<-agrInN

    #agrOutN<-FALSE
    #sumOut<-apply(adj_Matr,1,sum)
    
    #indOutN<-which(sumOut==0)
    #if(length(indOutN)==0)
    #{
    #  indOutN<-0
    #}
    
    #agrOutN<-setequal(outNodes, indOutN)
   
    #agrInOut<-(agrInN && agrOutN)


    df_adjMat<-convAdjTop(new_top, numb_cgnodes, 1:numb_cgnodes)#Mat(adj_Matr)
    grAdj<-graph_from_data_frame(df_adjMat, directed = FALSE, vertices = NULL)
      
    moveCond<-(length(V(grAdj)$name)<numb_cgnodes)|| (agrInOut==FALSE)
     
    if(moveCond==TRUE)
           {
             ntrials<-ntrials+1
             new_top<-old_top
           }

    if (moveCond==FALSE)
           {
             vconn<-vertex_connectivity(grAdj, checks = TRUE)
             
             if(vconn>0) 
                 {
                   connGraph<-TRUE       
                 }
             if (vconn==0)
                 {
                   ntrials<-ntrials+1
                   new_top<-old_top
                 }                                 
               
           } #if(moveCond==FALSE) 
   
          
   }#if(s!=old..)

   
}#while(ntr<=10 && connGr=F


####DO EXHAUSTIVE SEARCH FOR AN EDGE CHANGE
edgeTr<-1

  while ((connGraph==FALSE)&& (edgeTr<=numb_posEdges))
   {
     combEd<-indEd[edgeTr]
     s<-combEd%%3
     if( s!=0)
       {
        chEdge<-(combEd%/%3)+1
       }
     if(s==0)
       {
        chEdge<-(combEd%/%3)
       }

     
     if (s==old_top[chEdge])#force a change
       {
        if (s==0)
          {
           indS<-3
          }
        if(s!=0)
          {
           indS<-s
          }
    
        probS<-adj_matrPr[[chEdge]]
  
        if(sum(probS[-indS])>0)
          {
           newV<-c(1,2,0)
           newV<-newV[-indS]
           s<-sample(newV, size=1, prob=probS[-indS])
          }
       }

      if (s!=old_top[chEdge])
       {
        new_top[chEdge]<-s
        adj_Matr<-matrix(new_top, nrow=numb_cgnodes, ncol=numb_cgnodes, byrow=TRUE)    
     
        agrInN<-FALSE
        sumInc<-apply(adj_Matr,2,sum)
   
        indInN<-which(sumInc==0)
        if(length(indInN)==0)
        {
          indInN<-0
        }
        
        agrInN<-setequal(inNodes, indInN)
        
        agrInOut<-agrInN

        #agrOutN<-FALSE
        #sumOut<-apply(adj_Matr,1,sum)
    
        #indOutN<-which(sumOut==0)
        #if(length(indOutN)==0)
        #{
        # indOutN<-0
        #}
        
        #agrOutN<-setequal(outNodes, indOutN)

        #agrInOut<-(agrInN && agrOutN)


        df_adjMat<-convAdjTop(new_top, numb_cgnodes, 1:numb_cgnodes)#Mat(adj_Matr)
        grAdj<-graph_from_data_frame(df_adjMat, directed = FALSE, vertices = NULL)
      
        moveCond<-(length(V(grAdj)$name)<numb_cgnodes)|| (agrInOut==FALSE)
     
        if (moveCond==TRUE)
           {
             edgeTr<-edgeTr+1
             new_top<-old_top
           }
        if (moveCond==FALSE)
           {
             vconn<-vertex_connectivity(grAdj, checks = TRUE)

             if(vconn>0) 
                 {
                   connGraph<-TRUE       
                 }
             if (vconn==0)
                 {
                   edgeTr<-edgeTr+1
                   new_top<-old_top                    
                 }
               
           } #if(moveCond==FALSE) 
           
       
       
       }#if(s!=old..)
     
   }#while(connGraph==F && edgeTr<=, i.e. end of exhaustive search


if(edgeTr>numb_posEdges)
  {
    backtrInd<-backtrInd-1
  }

}###end backtracking while loop

if(connGraph==TRUE)
  {
     return(new_top)
  }

if(backtr<1)
  {
     return (rep(0, numb_cgnodes*numb_cgnodes))
  }
    


}#end function




inNWsimilarityRefCGr<-function(dataRow, clusterRef, cenMedRef, cutOffM, gene_list, inNodes, topol_cgr)
{
#topol_cgr<-newt_top
#dataRow doesn't include expression(noise) of input genes 
numb_genes<-dim(dataRow)[1]#=number of non-input species(genes)
modelsRef<-dim(dataRow)[2]
modelsCGr<-10000
numbCGRn<-length(gene_list)#it includes input and output nodes
numbClusterRef<-length(unique(clusterRef))#it doesn't include input nodes
scoresOut<-rep(0,4)


df_lastTop<-convAdjTop(topol_cgr, numbCGRn, 1:numbCGRn)

		
 #SIMULATE CGR-CIRCUIT
 #set.seed(42)
 racCGR<-sracipeSimulate(circuit = df_lastTop, numModels = modelsCGr, plots = FALSE, 
                        integrateStepSize = 0.03, simulationTime = 120)


 racNCGr<-sracipeNormalize(racCGR)

 dataRowGeneCGR<-assay(racNCGr,1)#This is our coarse-grained data
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
 
#write.table(t(dim(dataColGeneCGR)), file=fErr, append=TRUE, sep="\t", row.names=F, col.names=F)

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
#scoresOut[1]<-sum(clusterCutSimM[,3])
scoresOut[1]<-length(indZe)

sizeCl<-rep(0,numbClusterRef)
sizeClref<-rep(0,numbClusterRef)
avgCl<-rep(0,numbClusterRef)
chScore<-FALSE
for(ila in 1:numbClusterRef)
{
indCAv<-which(clusterCutSimM[,1]==ila)
sizeCl[ila]<-length(indCAv)
sizeClref[ila]<-length(which(clusterRef==ila))
avgCl[ila]<-median(clusterCutSimM[indCAv,3])
if(sizeCl[ila]<(0.5*sizeClref[ila]))
{
 chScore<-TRUE
}

}


scoresOut[2]<-(modelsCGr/numbClusterRef)*sum(avgCl)
scoresOut[3]<-sum(clusterCutSimM[indZe,3])
scoresOut[4]<-scoresOut[2]+scoresOut[3]
if(chScore==TRUE)
{
 scoresOut[4]<-scoresOut[4]+50000
}

return(c(scoresOut[4],scoresOut[1]))

} 



difLen<-function(samptop, reftop)
{
 return(length(which(samptop!=reftop)))
}


newdelEdges_top<-function(samp_top, numb_genes, numb_del, pos_del, list_topl)
{
#deletes numb_del edges from a samp_top that has numb_genes
#samp_top is given as a linearized adj matrix of size numbGenes*numbGenes
#pos_del are the edges(as index in linearized adj matrix) that can be deleted 
#while trying to maintain connectivity
#make sure the obtained topology is not in list_topl 


#n<-length(list_topl)
npossDel<-length(pos_del)
allDel<-NULL
probV<-rep(1/npossDel, npossDel)
backTr<-1

while ((length(allDel)<numb_del)&& (backTr>0))
{
 
 vconn<-0
 zeroTot<-1
 indFind<-0
 vcont<-((vconn==0)|| (zeroTot>0))||(indFind!=0)
 
 indExit<-length(which(probV!=0))
 countDel<-1
 cprobV<-probV
 while(vcont &&(countDel<=indExit))
 {
    delOne<-sample(1:npossDel, size=1, prob=cprobV)
    cprobV[delOne]<-0
    tallDel<-c(allDel, delOne)

    tnew_adj<-samp_top
    tnew_adj[pos_del[tallDel]]<-0
    tnew_top<-convAdjTop(tnew_adj, numb_genes) 
    grAdj<-graph_from_data_frame(tnew_top, directed = FALSE, vertices = NULL)
    vconn<-vertex_connectivity(grAdj, checks = TRUE)
    
    adj_Matr<-matrix(tnew_adj, nrow=numb_genes, ncol=numb_genes, byrow=TRUE)    
     
    sumInc<-apply(adj_Matr,2,sum)
    sumOut<-apply(adj_Matr,1,sum)
    zeroS<-length(which(sumInc==0))
    zeroOut<-length(which(sumOut==0))
    zeroTot<-zeroS

    if( length(allDel)==(numb_del-1))
    {
    indFind<-find_topology(list_topl, tnew_adj, 1,length(list_topl))
    }

    vcont<-((vconn==0)|| (zeroTot>0))||(indFind!=0)
 
    countDel<-countDel+1
 }     

 if ((countDel>indExit)&& (vcont=TRUE))
  {
   
   laDel<-length(allDel)
   if (laDel>0)
    { 
       #probV[allDel[laDel]]<-1/npossDel
       allDel<-allDel[-laDel]
       backTr<-length(allDel)
    }
   if (laDel==0)
    {
       backTr<-0
    }
  }
 else
  {
     allDel<-c(allDel, delOne)
     probV[delOne]<-0
     backTr<-length(allDel)
  }

}

if (backTr>0)
{
fin_top<-samp_top
fin_top[pos_del[allDel]]<-0
#adj_fin<-convTopAdj(fin_top, numb_genes)
return(fin_top)
}

if (backTr==0)
{
 return(rep(0, numb_genes*numb_genes))
}

}






signChanges_top<-function(samp_top, numbGenes, numbChanges, pos_changes, list_top)
{
#tries to do a number of numbChanges sign changes to edges of samp_top
#samp_top is given as a linearized adj matrix of size numbGenes*numbGenes
#pos_changes are the edges(as index in linearized adj matrix) for which 
#we can flip the sign
#need to make sure the obtained top is not in list_top
#assumes numbChanges>0
#tries to do it 3 times, it goes 1 level deep

#n<-length(list_top)
npossCh<-length(pos_changes)

m<-1
vcont<-TRUE

while(vcont && (m<=3))
{

probCV<-rep(1/npossCh, npossCh)
backTr<-numbChanges

cTop<-samp_top
chOne<-sample(1:npossCh, size=numbChanges)

cTop[pos_changes[chOne]]<-3-cTop[pos_changes[chOne]]
adj_Matr<-matrix(cTop, nrow=numbGenes, ncol=numbGenes, byrow=TRUE)  

sumInc<-apply(adj_Matr,2,sum)
sumOut<-apply(adj_Matr,1,sum)
zeroS<-length(which(sumInc==0))
zeroOut<-length(which(sumOut==0))
zeroTot<-zeroS

indF<-find_topology(list_top, cTop, 1,length(list_top))

vcont<-(zeroTot>0)||(indF!=0)
probCV[chOne]<-0

while(vcont && (backTr==numbChanges))
{
  #for (t in backTr:numbChanges)
  #{
   # cTop[pos_changes[chOne[t]]]<-3-cTop[pos_changes[chOne[t]]]
  #}
     cTop[pos_changes[chOne[numbChanges]]]<-3-cTop[pos_changes[chOne[numbChanges]]]

  #if (backTr==numbChanges)
  #{ 
    cprobCV<-probCV
    while (vcont && (sum(cprobCV)>0))
    {
      ccTop<-cTop
      singleCh<-sample(1:npossCh, size=1, prob=cprobCV)
      cprobCV[singleCh]<-0
      ccTop[pos_changes[singleCh]]<-3-ccTop[pos_changes[singleCh]]
      
      adj_Matr<-matrix(ccTop, nrow=numbGenes, ncol=numbGenes, byrow=TRUE)  
      sumInc<-apply(adj_Matr,2,sum)
      sumOut<-apply(adj_Matr,1,sum)
      zeroS<-length(which(sumInc==0))
      zeroOut<-length(which(sumOut==0))
      zeroTot<-zeroS

      indF<-find_topology(list_top, ccTop, 1,length(list_top))

      vcont<-(zeroTot>0)||(indF!=0)
         
     }

     if (vcont==TRUE)
       {
        backTr<-backTr-1
       }
     else
       {
        cTop<-ccTop
       }
   
  
   
}#end while(vcont && bactkTr==numbChanges)

m<-m+1

}#end while(vcont && m<=3)

if (vcont==TRUE)
{
  return(rep(0,numbGenes*numbGenes))
}
else
{
  return (cTop)
}

}#end function




######THIS FUNCTION TOGETHER WITH THE PREVIOUS THREE ARE TAKEN FROM GAHYBRID FILE

gaInitial_gen<-function(pert_top, gene_list, nGenes, numbNewTop)
{
#pert_top<-allMCs
#gene_list<-gene_6cg
#nGenes<-24
#numbNewTop<-50

  dfPert<-convAdjTop(pert_top,nGenes)
  adj_matrProb<-adjMatCGrProb(dfPert, gene_list)
  numb_cgnodes<-length(gene_list)

  #Find the densest topology
  #You'll delete edges from it and change edge signs

  inTop<-initialize_topology(adj_matrProb, numb_cgnodes)
  #init_topF<-convAdjTop(inTop,numb_cgnodes)

  #Find the edges that allow sign changes and those that allow deletion

  bothSigns<-NULL
  canDel<-NULL
  for (i in 1:(numb_cgnodes*numb_cgnodes))
  {
    probE<-adj_matrProb[[i]]
    if((probE[1]>0) && (probE[2]>0))
      {
        bothSigns<-c(bothSigns,i)
      }
    if(probE[3]>0)
      {
        canDel<-c(canDel,i)
      }
   }#end for  

   edgeInd<-which(inTop!=0)
   edgeCDel<-intersect(canDel, edgeInd)
   maxDel<-min(length(edgeCDel), length(edgeInd)-numb_cgnodes)
   #eDel<-is.element(edgeInd, canDel)   
   pos_Del<-edgeCDel

   #First we'll delete edges, then we'll make sign changes
   newT<-1
   list_topCG<-list()
   list_topCG[[1]]<-inTop

   #First produce 2*numbNewTop(<300) different topologies
   #numbNewTop will be 100 max
   iter<-0
   while((newT<2*numbNewTop)&& (iter<300))
   { 
     pTop<-inTop
     if (maxDel>0)
     {
      raUn<-runif(1)
      numbDel<-as.integer(raUn*maxDel)#this is number of deletions on inTop
      
      if (numbDel>0)
       {
          adjNew<-newdelEdges_top(inTop,numb_cgnodes, numbDel,pos_Del, list_topCG)
          pTop<-adjNew
       }

      }#end if(maxDel>0)

      #Now do sign changes, make sure you don't get duplicates
     
      if (sum(pTop)>0)
      {
         edgepTop<-which(pTop!=0)
         pos_Sign<-intersect(edgepTop, bothSigns)   
         #pos_Sign<-which(epSign) #these are the edges of pTop that admit sign changes
         maxSign<-length(pos_Sign)
         if (maxSign>0)
         {
          
           #Try to do sign changes 3 times, stop when you get new topology
           u<-1
           badCond<-TRUE
           while(badCond && (u<=3))
           {
             raUn<-runif(1)
             numbSign<-as.integer(raUn*maxSign)#this is number of sign changes on inTop
             if (numbSign>0)
             {
                ppTop<-signChanges_top(pTop, numb_cgnodes, numbSign, pos_Sign, list_topCG)
                if(sum(ppTop)!=0)
                {
                  pTop<-ppTop
                  badCond<-FALSE
                }
                
              }
              else
              {
                badCond<-FALSE
              }

              u<-u+1
            }#end while

           }#end if (maxSign>0)

       }#end if (sum(pTop)>0)

       if (sum(pTop)>0)
        {
           newT<-newT+1
           list_topCG[[newT]]<-pTop
        }

       iter<-iter+1

   }#end while(newT<2*numbNewTop...

   distM<-matrix(0, nrow=newT, ncol=newT)
   for(i in 1:newT)
    {
      for(j in 1:newT)
      {
        distM[i,j]<-difLen(list_topCG[[i]], list_topCG[[j]])
      }
    }

   sumDist<-apply(distM, 1, sum)
   sortDist<-sort(sumDist, decreasing=TRUE, index.return=TRUE)
   numbGet<-min(numbNewTop, newT)
   resM<-matrix(0, nrow=numbGet, ncol=numb_cgnodes*numb_cgnodes)
   for(lto in 1:numbGet)
    {
      resM[lto,]<-list_topCG[[sortDist$ix[lto]]]
    }

  return(resM)

}






###################HERE STARTS THE MAIN PROGRAM FOR SYNTHETIC DATA
###ASSUME WE HAVE ONLY THE TOPOLOGY FILE
###AND WE KNOW NUMBcLUSTrEF AND NUMBcgrNODES

expdata_top<-read.csv("Repress3by3ext.csv", header=T)
dataCol1<-read.csv("dataExt22_clusteredCorr.csv", header=F)#work with already clustered data
coln<-as.character(dataCol1[1,])
dataColGeneexp<-as.matrix(dataCol1[-c(1),])
colnames(dataColGeneexp)<-coln

dataRowGeneexp<-t(dataColGeneexp)

rownames(dataRowGeneexp)<-coln


#dataRowGene is our dataReference-"exptl" data

cl7Size<-c(2344, 934, 348, 2307, 944, 1094, 2029)
cumClSize<-cumsum(cl7Size)
clusterCutRef<-rep(0, 10000)
clusterCutRef[1:cumClSize[1]]<-1
clusterCutRef[(1+cumClSize[1]):(cumClSize[2])]<-2
clusterCutRef[(cumClSize[2]+1):(cumClSize[3])]<-3
clusterCutRef[(cumClSize[3]+1):(cumClSize[4])]<-4
clusterCutRef[(cumClSize[4]+1):(cumClSize[5])]<-5
clusterCutRef[(cumClSize[5]+1):(cumClSize[6])]<-6
clusterCutRef[(cumClSize[6]+1):(cumClSize[7])]<-7

###DO GENE GROUPING IN numbCGRnodes GROUPS
m1<-apply(dataColGeneexp[1:cumClSize[1],],2,median)
m2<-apply(dataColGeneexp[(1+cumClSize[1]):(cumClSize[2]),],2,median)
m3<-apply(dataColGeneexp[(1+cumClSize[2]):(cumClSize[3]),],2,median)
m4<-apply(dataColGeneexp[(1+cumClSize[3]):(cumClSize[4]),],2,median)
m5<-apply(dataColGeneexp[(1+cumClSize[4]):(cumClSize[5]),],2,median)
m6<-apply(dataColGeneexp[(1+cumClSize[5]):(cumClSize[6]),],2,median)
m7<-apply(dataColGeneexp[(1+cumClSize[6]):(cumClSize[7]),],2,median)

mMat41<-rbind(rbind(rbind(m1,m2),m3),m4)
mMat42<-rbind(rbind(m5,m6),m7)
mMat<-rbind(mMat41,mMat42)
cMed<-as.matrix(mMat)
colnames(cMed)<-coln
corMed<-cor(cMed,method="pearson")
gene_3cg<-bottom_up_gene_grouping(corMed,3)

gene_grouping<-gene_3cg


#numbCGRnodes<-3

cenDf<-read.csv("Myd20_CenMedRef_Ex22Corr.csv", header=T)
myd_cenMedRef<-as.matrix(cenDf) 
cutOff<-read.csv("Myd20_CutOffM_Ex22Corr.csv", header=T)
myd_cutOffM<-as.matrix(cutOff)

colnames(myd_cenMedRef)<-coln

fileOutTop<-"Ex22corr_PT25w_20perc78_tops.txt"
fileOutErrAcc<-"Ex22corr_PT25w_20perc78_ScoreAccRate.txt"
fileOutStuff<-"Ex22corr_PT25w_20perc78_Stuff.txt"
fileAllSamp<-"Ex22corr_PT25w_20perc78_tops_allSampled.txt"
fileAllInitTop<-"Ex22corr_PT25w_20perc78_tops_initial.txt"
fileSampUnsc<-"Ex22corr_PT25w_20perc78_tops_SampUnsc.txt"

numb_Sim<-2
ordSim<-9


allPert<-read.csv("PerturbedTop_top10_all_Ex22.csv",header=F)
pertM<-as.matrix(allPert[(numb_Sim*(ordSim-1)+1):(numb_Sim*ordSim),1:81])

selSeq<-seq(from=2*numb_Sim*(ordSim-1)+1, to=2*numb_Sim*ordSim)

allInit<-read.csv("All_initial_tops_Ex22.csv",header=F)
inTopsM<-as.matrix(allInit[selSeq,])


adjProblist<-list()
for(i in 1:numb_Sim)
{
dfPert<-convAdjTop(pertM[i,], dim(dataRowGeneexp)[1])
adjProblist[[i]]<-adjMatCGrProb(dfPert, gene_grouping)
}


numbThr<-50

logAlpha<-log(0.4)
tempArray<-c(1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.5,2.0,2.5,3.0,3.5,4.0,6,9,11,13,20,28,40,55,70,90,120)
tempArrayV<-rep(tempArray,numb_Sim)

tempArrayL<-list() 


nChainsV<-rep(length(tempArray),numb_Sim) #20=initial number of chains
nChains<-sum(nChainsV)
numbCGRn<-length(gene_grouping)

posM<-NULL

list_init<-list()
numIn<-NULL
for(t in 1:numb_Sim)
{
inTopr<-gaInitial_gen(pertM[t,], gene_grouping, dim(dataRowGeneexp)[1], 70)
inToptry<-rbind(inTopsM[(2*t-1):(2*t),],inTopr)
inTopw<-inToptry[!duplicated(inToptry),]
inTop<-as.matrix(inTopw)
numNew<-dim(inTop)[1]
numIn<-c(numIn, numNew)

list_init[[t]]<-as.matrix(inToptry)


write.table(inTop, file = fileAllInitTop, append = TRUE, sep = "\t", row.names = F, col.names = F)

posM<-rbind(posM,inToptry[1:nChainsV[t],])

tempArrayL[[t]]<-tempArray
}

write.table(t(numIn),file=fileOutStuff,append=TRUE, sep="\t", row.names=F, col.names=F)
  


cl<-makeCluster(numbThr)#takes 52 minutes to run
registerDoParallel(cl)

clusterExport(cl, varlist=c('dataRowGeneexp', 'clusterCutRef','myd_cenMedRef',  'myd_cutOffM', 'gene_grouping'))


#####CALCULATE SCORES INITIAL TOP'S THEN WRITE THEM OUT
####BEFORE THAT WRITE THEM IN UNSCORED FILE

unscSamp<-cbind(rep(0,nChains), posM)
write.table(unscSamp, file=fileSampUnsc, append=TRUE, sep="\t", row.names=F, col.names=F)


list_sampled<-list()
scSampV<-NULL
thrownOut<-NULL

redIn<-posM[!duplicated(posM),]
realIn<-dim(redIn)[1]

for(i in 1:realIn)
{
list_sampled[[i]]<-redIn[i,]
}

calc_posM<-NULL 
for(t in 1:realIn)
{
partM<-matrix(rep(redIn[t,],5),nrow=5,ncol=dim(redIn)[2],byrow=T)
calc_posM<-rbind(calc_posM, partM)
}

nSimul<-5*realIn
nrep<-ceiling(nSimul/numbThr)


red_lastSc<-NULL
if(nrep>1)
{
   for(s in 1:(nrep-1))
   {
    resPar<-foreach(j=1:numbThr, .combine=c, .packages=c('sRACIPE','igraph')) %dopar% inNWsimilarityRefCGr(dataRowGeneexp, 
         clusterCutRef, myd_cenMedRef, myd_cutOffM, gene_grouping, 0, calc_posM[(s-1)*numbThr+j,]) 
     
    vecRes<-rep(0,numbThr*2)
    for(m in 1:(numbThr*2))
    {
        vecRes[m]<-resPar[[m]]
    }
    matRes<-matrix(vecRes,nrow=numbThr,ncol=2,byrow=TRUE)
    red_lastSc<-rbind(red_lastSc,matRes)
   }
}
  
if(nSimul>((nrep-1)*numbThr+1))
{
  resPar2<-foreach(j=((nrep-1)*numbThr+1):nSimul, .combine=c, .packages=c('sRACIPE','igraph')) %dopar% inNWsimilarityRefCGr(dataRowGeneexp, 
         clusterCutRef, myd_cenMedRef, myd_cutOffM, gene_grouping, 0, calc_posM[j,]) 
    
   vecRes<-rep(0,(nSimul-(nrep-1)*numbThr)*2)
   for(m in 1:((nSimul-(nrep-1)*numbThr)*2))
   {
     vecRes[m]<-resPar2[[m]]
   }
   matRes<-matrix(vecRes,nrow=nSimul-(nrep-1)*numbThr,ncol=2,byrow=TRUE)
   red_lastSc<-rbind(red_lastSc, matRes)
}

for(t in 1:realIn)
{
  pMa<-red_lastSc[(5*t-4):(5*t),]
  scSampV<-c(scSampV, mean(pMa[,1]))
  thrownOut<-c(thrownOut, pMa[5,2])
}

last_scM<-NULL
for(i in 1:(dim(posM)[1]))
{
indF<-find_topology(list_sampled, posM[i,], 1, realIn)
last_scM<-rbind(last_scM, c(scSampV[indF],thrownOut[indF]))
}

selSeq<-c(1)
for(i in 1:(numb_Sim-1))
{
nextNu<-nChainsV[i]+selSeq[i]
selSeq<-c(selSeq,nextNu)
}

xoutIn<-cbind(posM[selSeq,], last_scM[selSeq,])
matInd<-matrix(rep(0,2*numb_Sim),nrow=numb_Sim, ncol=2,byrow=T)
matInd[,2]<-1:numb_Sim
xout<-cbind(matInd,xoutIn)
write.table(xout, file=fileOutTop, append=TRUE, sep="\t", row.names=F, col.names=F)

xouts<-cbind(rep(0,realIn),cbind(cbind(redIn,scSampV),thrownOut))
write.table(xouts, file=fileAllSamp, append=TRUE, sep="\t", row.names=F, col.names=F)

errAc<-rep(0, 2*nChainsV[1]+1)
for(i in 1:nChainsV[1])
{
 errAc[2*i]<-last_scM[i,1]
}
write.table(t(errAc), file=fileOutErrAcc, append=TRUE, sep="\t", row.names=F, col.names=F)

if(numb_Sim>1)
{

for(t in 2:numb_Sim)
{
  errAc<-rep(0, 2*nChainsV[t]+1)
  for(i in 1:nChainsV[t])
  {
    errAc[2*i]<-last_scM[sum(nChainsV[1:(t-1)])+i,1]
  }
  write.table(t(errAc), file=fileOutErrAcc, append=TRUE, sep="\t", row.names=F, col.names=F)
}

}


######START THE ITERATIONS

for(itAdd in 1:3)
{

numbSwapsA<-rep(0,nChains+2*numb_Sim)
accRate<-rep(0,nChains)
nUp<-list()
nDown<-list()
upDown<-list()
permArray<-list()
swEven<-list()

logSwapA<-list()
diffA<-list()
for(i in 1:numb_Sim)
{
permArray[[i]]<-1:nChainsV[i]
swEven[[i]]<-seq(from=nChainsV[i]-1, to=1, by=-1)
nUp[[i]]<-rep(0,nChainsV[i])
nDown[[i]]<-rep(0,nChainsV[i])
upDown[[i]]<-rep(0,nChainsV[i])
nUp[[i]][1]<-1
nDown[[i]][nChainsV[i]]<-1
upDown[[i]][1]<-1
upDown[[i]][nChainsV[i]]<--1
logSwapA[[i]]<-rep(0,nChainsV[i]-1)
diffA[[i]]<-rep(0,nChainsV[i])
}

endSeq<-rep(0,numb_Sim)
endSeq[1]<-nChainsV[1]
for(t in 2:numb_Sim)
{
  endSeq[t]<-endSeq[t-1]+nChainsV[t]
}


selSeq<-c(1)
for(i in 1:(numb_Sim-1))
{
  nextNu<-nChainsV[i]+selSeq[i]
  selSeq<-c(selSeq,nextNu)
}


indSw<-c(0)
for(u in 2:numb_Sim)
{
  indSw<-c(indSw,sum(nChainsV[1:(u-1)]))
}

numbIt<-50*itAdd

###START THE ITERATIONS

for(itn in 1:numbIt)
{
 #update the topologies of each chain
 newTop<-NULL
 indNew<-NULL
 nSamp<-length(list_sampled)

 new_scM<-last_scM
 newPosM<-posM


 for(i in 1:nChains)
 {
    inAdj<-1
    for(t in 2:numb_Sim)
    {
      if(i>endSeq[t-1])
      {
        inAdj<-t
      }
    }
    adjProbM<-adjProblist[[inAdj]]
    new_top<-inNrepsample_topology(adjProbM, numbCGRn, 0, posM[i,])
    indF<-find_topology(list_sampled, new_top, 1, nSamp)
    if(sum(new_top)>0)
    {
     if(indF==0)
     {
       newTop<-rbind(newTop,new_top)
       indNew<-c(indNew,i)           
     }
     
     newPosM[i,]<-new_top

    }#end if (sum(new_top)>0)
   
    if(sum(new_top)==0)
    {
      newPosM[i,]<-posM[i,]
      new_scM[i,]<-last_scM[i,]
    }
  
 }
 
 nUpd<-length(indNew)

 
 if(nUpd>0)
 {
  #newTop<-as.matrix(newTop)
  write.table(t(c(itn,indNew)),file=fileOutStuff,append=TRUE, sep="\t", row.names=F, col.names=F)
  
  unscSamp<-cbind(rep(itn,nUpd), newTop)
  write.table(unscSamp, file=fileSampUnsc, append=TRUE, sep="\t", row.names=F, col.names=F)

  redM<-newTop[!duplicated(newTop),]
  realSim<-as.integer(length(redM)/(numbCGRn*numbCGRn))
  if(realSim==1)
  {
    redM<-matrix(redM,nrow=1,ncol=numbCGRn*numbCGRn,byrow=T)
  }
  
  calc_posM<-NULL 
  for(t in 1:realSim)
  {
    partM<-matrix(rep(redM[t,],5),nrow=5,ncol=dim(redM)[2],byrow=T)
    calc_posM<-rbind(calc_posM, partM)
  }

  nSimul<-5*realSim
  nrep<-ceiling(nSimul/numbThr)


  red_newSc<-NULL
  if(nrep>1)
    {
     for(s in 1:(nrep-1))
     {
      resPar<-foreach(j=1:numbThr, .combine=c, .packages=c('sRACIPE','igraph')) %dopar% inNWsimilarityRefCGr(dataRowGeneexp, 
         clusterCutRef, myd_cenMedRef, myd_cutOffM, gene_grouping, 0, calc_posM[(s-1)*numbThr+j,]) 
     
      vecRes<-rep(0,numbThr*2)
      for(m in 1:(numbThr*2))
      {
         vecRes[m]<-resPar[[m]]
      }
      matRes<-matrix(vecRes,nrow=numbThr,ncol=2,byrow=TRUE)
      red_newSc<-rbind(red_newSc,matRes)
     }
    }
  
   if(nSimul>((nrep-1)*numbThr+1))
    {
     resPar2<-foreach(j=((nrep-1)*numbThr+1):nSimul, .combine=c, .packages=c('sRACIPE','igraph')) %dopar% inNWsimilarityRefCGr(dataRowGeneexp, 
         clusterCutRef, myd_cenMedRef, myd_cutOffM, gene_grouping, 0, calc_posM[j,]) 
    
     vecRes<-rep(0,(nSimul-(nrep-1)*numbThr)*2)
     for(m in 1:((nSimul-(nrep-1)*numbThr)*2))
     {
       vecRes[m]<-resPar2[[m]]
     }
     matRes<-matrix(vecRes,nrow=nSimul-(nrep-1)*numbThr,ncol=2,byrow=TRUE)
     red_newSc<-rbind(red_newSc, matRes)
    }

    for(t in 1:realSim)
    {
     pMa<-red_newSc[(5*t-4):(5*t),]
     scSampV<-c(scSampV, mean(pMa[,1]))
     thrownOut<-c(thrownOut, pMa[5,2])
    }

    xouts<-cbind(rep(itn,realSim),cbind(cbind(redM,scSampV[(nSamp+1):(nSamp+realSim)]),thrownOut[(nSamp+1):(nSamp+realSim)]))
    write.table(xouts, file=fileAllSamp, append=TRUE, sep="\t", row.names=F, col.names=F)

    for(t in 1:realSim)
    {
      nSamp<-nSamp+1
      list_sampled[[nSamp]]<-redM[t,]
    }

    new_scInM<-NULL

    for(i in 1:nUpd)
    {
     indF<-find_topology(list_sampled, newTop[i,], 1, nSamp)
     new_scInM<-rbind(new_scInM, c(scSampV[indF],thrownOut[indF]))
    }

    new_scM[indNew,]<-new_scInM  
 
  }#end if(nUpd>0)

  nnew<-setdiff(1:nChains,indNew)
  if(length(nnew)>0)
  {
   for(t in nnew)
   {
    indF<-find_topology(list_sampled, newPosM[t,], 1, nSamp)
    new_scM[t,1]<-scSampV[indF]
    new_scM[t,2]<-0
   }
  }

 
 #SEE IF WE ACCEPT NEW SCORES AND THEN DO SWAPS
 for( i in 1:nChains)
 {
   if(last_scM[i,1]==new_scM[i,1])
   {
      accRate[i]<-accRate[i]-1
   }
   rnU<-runif(1)
   if(rnU<exp((last_scM[i,1]-new_scM[i,1])/tempArrayV[i]))
      {
        last_scM[i,]<-new_scM[i,]
        posM[i,]<-newPosM[i,]
        accRate[i]<-accRate[i]+1
      }
    
 } 
 
 #WRITE SCORES BEFORE THE SWAPS

 errAc<-rep(itn, 2*nChainsV[1]+1)
 for(i in 1:nChainsV[1])
 {
   errAc[2*i]<-last_scM[i,1]
   errAcc[2*i+1]<-accRate[i]
 }
 write.table(t(errAc), file=fileOutErrAcc, append=TRUE, sep="\t", row.names=F, col.names=F)

 if(numb_Sim>1)
 {
   for(t in 2:numb_Sim)
   {
    errAc<-rep(itn, 2*nChainsV[t]+1)
    for(i in 1:nChainsV[t])
    {
      errAc[2*i]<-last_scM[sum(nChainsV[1:(t-1)])+i,1]
      errAcc[2*i+1]<-accRate[sum(nChainsV[1:(t-1)])+i]
    }
    write.table(t(errAc), file=fileOutErrAcc, append=TRUE, sep="\t", row.names=F, col.names=F)
   }
 }

 
 #NOW DO THE SWAPS

 decR<-1
 rDos<-runif(1)
 if(rDos<decR)
 {
  #DO All CHAINS (10,9),(9,8),...,(2,1)
  for(s in 1:numb_Sim)
  {
   for(t in swEven[[s]])
   {
     rnUn<-runif(1)
     swRate<-exp((last_scM[indSw[s]+t+1,1]-last_scM[indSw[s]+t,1])*(1/tempArrayV[indSw[s]+t+1]-1/tempArrayV[indSw[s]+t]))
     if(rnUn<swRate)#if accept swap
     {
        swErr<-last_scM[indSw[s]+t+1,]
        last_scM[indSw[s]+t+1,]<-last_scM[indSw[s]+t,]
        last_scM[indSw[s]+t,]<-swErr
        swPos<-posM[indSw[s]+t+1,]
        posM[indSw[s]+t+1,]<-posM[indSw[s]+t,]
        posM[indSw[s]+t,]<-swPos
     
     
        if((t==1)&&((permArray[[s]][2]==1)&&(upDown[[s]][1]==-1)))
         {
           numbSwapsA[nChains+2*s]<-numbSwapsA[nChains+2*s]+1
         }    
        if(t==1)
         {
           upDown[[s]][permArray[[s]][2]]<-1
         }
        if((t==(nChainsV[s]-1))&&((permArray[[s]][t]==1)&&(upDown[[s]][1]==1)))
         {
           numbSwapsA[nChains+2*s-1]<-numbSwapsA[nChains+2*s-1]+1
         }     
        if(t==(nChainsV[s]-1))
         {
           upDown[[s]][permArray[[s]][t]]<--1
         }
        if(upDown[[s]][permArray[[s]][t]]==1)
         {
           nUp[[s]][t+1]<-nUp[[s]][t+1]+1
         }
        if(upDown[[s]][permArray[[s]][t+1]]==-1)
	   {
           nDown[[s]][t]<-nDown[[s]][t]+1
         }

        swPerm<-permArray[[s]][t]
        permArray[[s]][t]<-permArray[[s]][t+1]
        permArray[[s]][t+1]<-swPerm

        numbSwapsA[indSw[s]+t]<-numbSwapsA[indSw[s]+t]+1
     }

     if(swRate<1)
      {
        logSwapA[[s]][t]<-logSwapA[[s]][t]+log(swRate)
      }
   }#end for(t in ...)

   numbSwapsA[endSeq[s]]<-numbSwapsA[endSeq[s]]+1
  
  }#end for(s in ...)
 }#end if(rDos<...)


#### NOW WRITE THE POSITIONS SCORES ARRAYS
  
  xoutIn<-cbind(posM[selSeq,], last_scM[selSeq,])
  matInd<-matrix(rep(itn,2*numb_Sim),nrow=numb_Sim, ncol=2,byrow=T)
  matInd[,2]<-1:numb_Sim
  xout<-cbind(matInd,xoutIn)
  write.table(xout, file=fileOutTop, append=TRUE, sep="\t", row.names=F, col.names=F)

  errAc<-rep(itn, 2*nChainsV[1]+1)
  for(i in 1:nChainsV[1])
  {
   errAc[2*i]<-last_scM[i,1]
   errAcc[2*i+1]<-accRate[i]
  }
  write.table(t(errAc), file=fileOutErrAcc, append=TRUE, sep="\t", row.names=F, col.names=F)

  if(numb_Sim>1)
  {
   for(t in 2:numb_Sim)
   {
    errAc<-rep(itn, 2*nChainsV[t]+1)
    for(i in 1:nChainsV[t])
    {
      errAc[2*i]<-last_scM[sum(nChainsV[1:(t-1)])+i,1]
      errAcc[2*i+1]<-accRate[sum(nChainsV[1:(t-1)])+i]
    }
    write.table(t(errAc), file=fileOutErrAcc, append=TRUE, sep="\t", row.names=F, col.names=F)
   }
  }

  
  
  if(itn==numbIt)
   {
     write.table(t(numbSwapsA), file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)
     for(s in 1:numb_Sim)
     {
       write.table(t(permArray[[s]]), file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)
       write.table(t(nUp[[s]]), file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)
       write.table(t(nDown[[s]]), file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)
       write.table(t(upDown[[s]]), file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)
     }
   }


}#end for (itn in 1:...)

for(s in 1:numb_Sim)
{
  for (t in 1:(nChainsV[s]-1))
  {
    logSwapA[[s]][t]<-logSwapA[[s]][t]/numbSwapsA[endSeq[s]]
    if((nUp[[s]][t]+nDown[[s]][t])>0)
     {
       diffA[[s]][t]<-nUp[[s]][t]/(nUp[[s]][t]+nDown[[s]][t])
     }
  }
}

write.table("Log Swap rates and diffusion function are", file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)
for(s in 1:numb_Sim)
{
  write.table(t(logSwapA[[s]]), file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)
  write.table(t(diffA[[s]]), file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)
}

    

###ADD THE NEW TEMPERATURES
rLogs<-rep(0, nChains-1)
sumRV<-rep(0,numb_Sim)


write.table("RLogs are", file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)

for(s in 1:numb_Sim)
{
  for(i in 1:(nChainsV[s]-1))
  {
   if(numbSwapsA[indSw[s]+i]>0)
   {
     rLogs[indSw[s]+i]<-log(numbSwapsA[indSw[s]+i]/numbSwapsA[endSeq[s]])/logAlpha
   }
   else
   {
     rLogs[indSw[s]+i]<-logSwapA[[s]][i]/logAlpha
   }
  }
  write.table(t(rLogs[(indSw[s]+1):(indSw[s]+nChainsV[s])]), file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)

}


for(s in 1:numb_Sim)
{
  newTempA<-NULL
  for(i in 1:(nChainsV[s]-1)) 
  {    
     lrat<-as.integer(sqrt(rLogs[indSw[s]+i]))
     sumRV[s]<-sumRV[s]+lrat

     for(m in 1:(lrat+1))
     { 
       newT<-tempArrayV[indSw[s]+i]+(tempArrayV[indSw[s]+i+1]-tempArrayV[indSw[s]+i])*(m-1)/(lrat+1)
       newTempA<-c(newTempA, newT)
     }
  }
  newTempA<-c(newTempA, tempArrayV[indSw[s]+nChainsV[s]])
  tempArrayL[[s]]<-newTempA
}


write.table("SumRV is", file = fileOutStuff, append = TRUE, sep = "\t", row.names = F, col.names = F)
write.table(t(sumRV), file = fileOutStuff, append = TRUE, sep = "\t", row.names = F, col.names = F)
write.table("New temperature arrays are", file = fileOutStuff, append = TRUE, sep = "\t", row.names = F, col.names = F)
for(s in 1:numb_Sim)
{
  write.table(t(tempArrayL[[s]]), file = fileOutStuff, append = TRUE, sep = "\t", row.names = F, col.names = F)
}



if(sum(sumRV)>0)
{
  newTops<-NULL
  tempposM<-posM[1:nChainsV[1],]

  for(s in 1:numb_Sim)
  {
    if (sumRV[s]>0)
    { 
      newTops<-rbind(newTops,list_init[[s]][(nChainsV[s]+1):(nChainsV[s]+sumRV[s]),])
  
      tempposM<-rbind(tempposM,list_init[[s]][(nChainsV[s]+1):(nChainsV[s]+sumRV[s]),])
    }
  
    if(s<numb_Sim)
    {
      tempposM<-rbind(tempposM,posM[(endSeq[s]+1):(endSeq[s+1]),])
    }

  }

  posM<-tempposM

  unscSamp<-cbind(rep(0,as.integer(sum(sumRV))), newTops)
  write.table(unscSamp, file=fileSampUnsc, append=TRUE, sep="\t", row.names=F, col.names=F)


  ##CALCULATE SCORES FOR THE ADDED TOPOLOGIES
  
  inredM<-newTops[!duplicated(newTops),]
  #inrealSim<-dim(inredM)[1]
  inrealSim<-as.integer(length(inredM)/(numbCGRn*numbCGRn))
  if(inrealSim==1)
  {
      inredM<-matrix(inredM,nrow=1,ncol=numbCGRn*numbCGRn,byrow=T)
  }

  nTops<-NULL
  for(t in 1:inrealSim)
  {
    indF<-find_topology(list_sampled, inredM[t,], 1, nSamp)
    if(indF==0)
     {
       nTops<-rbind(nTops,inredM[t,])           
     }
  }

  realSim<-as.integer(length(nTops)/(numbCGRn*numbCGRn))
  if(realSim==1)
  {
      nTops<-matrix(nTops,nrow=1,ncol=numbCGRn*numbCGRn,byrow=T)
  }

  if(realSim>0)
  {
    calc_posM<-NULL 
    for(t in 1:realSim)
    {
      partM<-matrix(rep(nTops[t,],5),nrow=5,ncol=dim(nTops)[2],byrow=T)
      calc_posM<-rbind(calc_posM, partM)
    }

    nSimul<-5*realSim
    nrep<-ceiling(nSimul/numbThr)


    red_newSc<-NULL
    if(nrep>1)
    {
      for(s in 1:(nrep-1))
      {
       resPar<-foreach(j=1:numbThr, .combine=c, .packages=c('sRACIPE','igraph')) %dopar% inNWsimilarityRefCGr(dataRowGeneexp, 
         clusterCutRef, myd_cenMedRef, myd_cutOffM, gene_grouping, 0, calc_posM[(s-1)*numbThr+j,]) 
     
       vecRes<-rep(0,numbThr*2)
       for(m in 1:(numbThr*2))
       {
         vecRes[m]<-resPar[[m]]
       }
       matRes<-matrix(vecRes,nrow=numbThr,ncol=2,byrow=TRUE)
       red_newSc<-rbind(red_newSc,matRes)
      }
    }
  
    if(nSimul>((nrep-1)*numbThr+1))
    {
      resPar2<-foreach(j=((nrep-1)*numbThr+1):nSimul, .combine=c, .packages=c('sRACIPE','igraph')) %dopar% inNWsimilarityRefCGr(dataRowGeneexp, 
         clusterCutRef, myd_cenMedRef, myd_cutOffM, gene_grouping, 0, calc_posM[j,]) 
    
      vecRes<-rep(0,(nSimul-(nrep-1)*numbThr)*2)
      for(m in 1:((nSimul-(nrep-1)*numbThr)*2))
      {
        vecRes[m]<-resPar2[[m]]
      }
      matRes<-matrix(vecRes,nrow=nSimul-(nrep-1)*numbThr,ncol=2,byrow=TRUE)
      red_newSc<-rbind(red_newSc, matRes)
    }

    for(t in 1:realSim)
    {
     pMa<-red_newSc[(5*t-4):(5*t),]
     scSampV<-c(scSampV, mean(pMa[,1]))
     thrownOut<-c(thrownOut, pMa[5,2])
    }

    xouts<-cbind(rep(0,realSim),cbind(cbind(nTops,scSampV[(nSamp+1):(nSamp+realSim)]),thrownOut[(nSamp+1):(nSamp+realSim)]))
    write.table(xouts, file=fileAllSamp, append=TRUE, sep="\t", row.names=F, col.names=F)

    for(t in 1:realSim)
    {
      nSamp<-nSamp+1
      list_sampled[[nSamp]]<-nTops[t,]
    }

   }#enf if(realSim>0)

   last_scM<-NULL
   ntop<-dim(posM)[1]
   for(t in 1:ntop)
   {
     indF<-find_topology(list_sampled, posM[t,], 1, nSamp)
     last_scM<-rbind(last_scM, c(scSampV[indF],thrownOut[indF]))
   }

   
   tempArrayV<-NULL
   for(s in 1:numb_Sim)
   {
      tempArrayV<-c(tempArrayV,tempArrayL[[s]])
      nChainsV[s]<-nChainsV[s]+sumRV[s]
   }
   nChains<-sum(nChainsV)


}#end if(sum(sumRV)>0)


write.table("New number of temperatures are", file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)
write.table(t(nChainsV), file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)
   

}#end for (itAdd in 1:3)


###DO LAST RUN
#resR<-runChains(dataRowGeneexp, clusterCutRef, myd_cenMedRef, myd_varRef, myd_cutOffM, gene_6cg, adjProbM, list_sampled, 
#scSampV, posM, last_scM, tempArray, numbThr, nChains, 400, fileOutTop, fileOutErrAcc,fileOutStuff,fileAllSamp, 1)


numbSwapsA<-rep(0,nChains+2*numb_Sim)
accRate<-rep(0,nChains)
nUp<-list()
nDown<-list()
upDown<-list()
permArray<-list()
swEven<-list()

logSwapA<-list()
diffA<-list()
for(i in 1:numb_Sim)
{
permArray[[i]]<-1:nChainsV[i]
swEven[[i]]<-seq(from=nChainsV[i]-1, to=1, by=-1)
nUp[[i]]<-rep(0,nChainsV[i])
nDown[[i]]<-rep(0,nChainsV[i])
upDown[[i]]<-rep(0,nChainsV[i])
nUp[[i]][1]<-1
nDown[[i]][nChainsV[i]]<-1
upDown[[i]][1]<-1
upDown[[i]][nChainsV[i]]<--1
logSwapA[[i]]<-rep(0,nChainsV[i]-1)
diffA[[i]]<-rep(0,nChainsV[i])
}

endSeq<-rep(0,numb_Sim)
endSeq[1]<-nChainsV[1]
for(t in 2:numb_Sim)
{
  endSeq[t]<-endSeq[t-1]+nChainsV[t]
}


selSeq<-c(1)
for(i in 1:(numb_Sim-1))
{
  nextNu<-nChainsV[i]+selSeq[i]
  selSeq<-c(selSeq,nextNu)
}


indSw<-c(0)
for(u in 2:numb_Sim)
{
  indSw<-c(indSw,sum(nChainsV[1:(u-1)]))
}

numbIt<-1100

###START THE ITERATIONS

for(itn in 1:numbIt)
{
 #update the topologies of each chain
 newTop<-NULL
 indNew<-NULL
 nSamp<-length(list_sampled)

 new_scM<-last_scM
 newPosM<-posM


 for(i in 1:nChains)
 {
    inAdj<-1
    for(t in 2:numb_Sim)
    {
      if(i>endSeq[t-1])
      {
        inAdj<-t
      }
    }
    adjProbM<-adjProblist[[inAdj]]
    new_top<-inNrepsample_topology(adjProbM, numbCGRn, 0, posM[i,])
    indF<-find_topology(list_sampled, new_top, 1, nSamp)
    if(sum(new_top)>0)
    {
     if(indF==0)
     {
       newTop<-rbind(newTop,new_top)
       indNew<-c(indNew,i)           
     }
     
     newPosM[i,]<-new_top

    }#end if (sum(new_top)>0)
   
    if(sum(new_top)==0)
    {
      newPosM[i,]<-posM[i,]
      new_scM[i,]<-last_scM[i,]
    }
  
 }
 
 nUpd<-length(indNew)

 
 if(nUpd>0)
 {
  #newTop<-as.matrix(newTop)
  write.table(t(c(itn,indNew)),file=fileOutStuff,append=TRUE, sep="\t", row.names=F, col.names=F)
  
  unscSamp<-cbind(rep(itn,nUpd), newTop)
  write.table(unscSamp, file=fileSampUnsc, append=TRUE, sep="\t", row.names=F, col.names=F)

  redM<-newTop[!duplicated(newTop),]
  realSim<-as.integer(length(redM)/(numbCGRn*numbCGRn))
  if(realSim==1)
  {
    redM<-matrix(redM,nrow=1,ncol=numbCGRn*numbCGRn,byrow=T)
  }

   
  calc_posM<-NULL 
  for(t in 1:realSim)
  {
    partM<-matrix(rep(redM[t,],5),nrow=5,ncol=dim(redM)[2],byrow=T)
    calc_posM<-rbind(calc_posM, partM)
  }

  nSimul<-5*realSim
  nrep<-ceiling(nSimul/numbThr)


  red_newSc<-NULL
  if(nrep>1)
    {
     for(s in 1:(nrep-1))
     {
      resPar<-foreach(j=1:numbThr, .combine=c, .packages=c('sRACIPE','igraph')) %dopar% inNWsimilarityRefCGr(dataRowGeneexp, 
         clusterCutRef, myd_cenMedRef, myd_cutOffM, gene_grouping, 0, calc_posM[(s-1)*numbThr+j,]) 
     
      vecRes<-rep(0,numbThr*2)
      for(m in 1:(numbThr*2))
      {
         vecRes[m]<-resPar[[m]]
      }
      matRes<-matrix(vecRes,nrow=numbThr,ncol=2,byrow=TRUE)
      red_newSc<-rbind(red_newSc,matRes)
     }
    }
  
   if(nSimul>((nrep-1)*numbThr+1))
    {
     resPar2<-foreach(j=((nrep-1)*numbThr+1):nSimul, .combine=c, .packages=c('sRACIPE','igraph')) %dopar% inNWsimilarityRefCGr(dataRowGeneexp, 
         clusterCutRef, myd_cenMedRef, myd_cutOffM, gene_grouping, 0, calc_posM[j,]) 
    
     vecRes<-rep(0,(nSimul-(nrep-1)*numbThr)*2)
     for(m in 1:((nSimul-(nrep-1)*numbThr)*2))
     {
       vecRes[m]<-resPar2[[m]]
     }
     matRes<-matrix(vecRes,nrow=nSimul-(nrep-1)*numbThr,ncol=2,byrow=TRUE)
     red_newSc<-rbind(red_newSc, matRes)
    }

    for(t in 1:realSim)
    {
     pMa<-red_newSc[(5*t-4):(5*t),]
     scSampV<-c(scSampV, mean(pMa[,1]))
     thrownOut<-c(thrownOut, pMa[5,2])
    }

    xouts<-cbind(rep(itn,realSim),cbind(cbind(redM,scSampV[(nSamp+1):(nSamp+realSim)]),thrownOut[(nSamp+1):(nSamp+realSim)]))
    write.table(xouts, file=fileAllSamp, append=TRUE, sep="\t", row.names=F, col.names=F)

    for(t in 1:realSim)
    {
      nSamp<-nSamp+1
      list_sampled[[nSamp]]<-redM[t,]
    }

    new_scInM<-NULL

    for(i in 1:nUpd)
    {
     indF<-find_topology(list_sampled, newTop[i,], 1, nSamp)
     new_scInM<-rbind(new_scInM, c(scSampV[indF],thrownOut[indF]))
    }

    new_scM[indNew,]<-new_scInM  
 
  }#end if(nUpd>0)

  nnew<-setdiff(1:nChains,indNew)
  if(length(nnew)>0)
  {
   for(t in nnew)
   {
    indF<-find_topology(list_sampled, newPosM[t,], 1, nSamp)
    new_scM[t,1]<-scSampV[indF]
    new_scM[t,2]<-0
   }
  }

 
 #SEE IF WE ACCEPT NEW SCORES AND THEN DO SWAPS
 for( i in 1:nChains)
 {
   if(last_scM[i,1]==new_scM[i,1])
   {
      accRate[i]<-accRate[i]-1
   }
   rnU<-runif(1)
   if(rnU<exp((last_scM[i,1]-new_scM[i,1])/tempArrayV[i]))
      {
        last_scM[i,]<-new_scM[i,]
        posM[i,]<-newPosM[i,]
        accRate[i]<-accRate[i]+1
      }
    
 } 
 
 #WRITE SCORES BEFORE THE SWAPS

 errAc<-rep(itn, 2*nChainsV[1]+1)
 for(i in 1:nChainsV[1])
 {
   errAc[2*i]<-last_scM[i,1]
   errAcc[2*i+1]<-accRate[i]
 }
 write.table(t(errAc), file=fileOutErrAcc, append=TRUE, sep="\t", row.names=F, col.names=F)

 if(numb_Sim>1)
 {
   for(t in 2:numb_Sim)
   {
    errAc<-rep(itn, 2*nChainsV[t]+1)
    for(i in 1:nChainsV[t])
    {
      errAc[2*i]<-last_scM[sum(nChainsV[1:(t-1)])+i,1]
      errAcc[2*i+1]<-accRate[sum(nChainsV[1:(t-1)])+i]
    }
    write.table(t(errAc), file=fileOutErrAcc, append=TRUE, sep="\t", row.names=F, col.names=F)
   }
 }

 
 #NOW DO THE SWAPS

 decR<-1
 rDos<-runif(1)
 if(rDos<decR)
 {
  #DO All CHAINS (10,9),(9,8),...,(2,1)
  for(s in 1:numb_Sim)
  {
   for(t in swEven[[s]])
   {
     rnUn<-runif(1)
     swRate<-exp((last_scM[indSw[s]+t+1,1]-last_scM[indSw[s]+t,1])*(1/tempArrayV[indSw[s]+t+1]-1/tempArrayV[indSw[s]+t]))
     if(rnUn<swRate)#if accept swap
     {
        swErr<-last_scM[indSw[s]+t+1,]
        last_scM[indSw[s]+t+1,]<-last_scM[indSw[s]+t,]
        last_scM[indSw[s]+t,]<-swErr
        swPos<-posM[indSw[s]+t+1,]
        posM[indSw[s]+t+1,]<-posM[indSw[s]+t,]
        posM[indSw[s]+t,]<-swPos
     
     
        if((t==1)&&((permArray[[s]][2]==1)&&(upDown[[s]][1]==-1)))
         {
           numbSwapsA[nChains+2*s]<-numbSwapsA[nChains+2*s]+1
         }    
        if(t==1)
         {
           upDown[[s]][permArray[[s]][2]]<-1
         }
        if((t==(nChainsV[s]-1))&&((permArray[[s]][t]==1)&&(upDown[[s]][1]==1)))
         {
           numbSwapsA[nChains+2*s-1]<-numbSwapsA[nChains+2*s-1]+1
         }     
        if(t==(nChainsV[s]-1))
         {
           upDown[[s]][permArray[[s]][t]]<--1
         }
        if(upDown[[s]][permArray[[s]][t]]==1)
         {
           nUp[[s]][t+1]<-nUp[[s]][t+1]+1
         }
        if(upDown[[s]][permArray[[s]][t+1]]==-1)
	   {
           nDown[[s]][t]<-nDown[[s]][t]+1
         }

        swPerm<-permArray[[s]][t]
        permArray[[s]][t]<-permArray[[s]][t+1]
        permArray[[s]][t+1]<-swPerm

        numbSwapsA[indSw[s]+t]<-numbSwapsA[indSw[s]+t]+1
     }

     if(swRate<1)
      {
        logSwapA[[s]][t]<-logSwapA[[s]][t]+log(swRate)
      }
   }#end for(t in ...)

   numbSwapsA[endSeq[s]]<-numbSwapsA[endSeq[s]]+1
  
  }#end for(s in ...)
 }#end if(rDos<...)


#### NOW WRITE THE POSITIONS SCORES ARRAYS
  
  xoutIn<-cbind(posM[selSeq,], last_scM[selSeq,])
  matInd<-matrix(rep(itn,2*numb_Sim),nrow=numb_Sim, ncol=2,byrow=T)
  matInd[,2]<-1:numb_Sim
  xout<-cbind(matInd,xoutIn)
  write.table(xout, file=fileOutTop, append=TRUE, sep="\t", row.names=F, col.names=F)

  errAc<-rep(itn, 2*nChainsV[1]+1)
  for(i in 1:nChainsV[1])
  {
   errAc[2*i]<-last_scM[i,1]
   errAcc[2*i+1]<-accRate[i]
  }
  write.table(t(errAc), file=fileOutErrAcc, append=TRUE, sep="\t", row.names=F, col.names=F)

  if(numb_Sim>1)
  {
   for(t in 2:numb_Sim)
   {
    errAc<-rep(itn, 2*nChainsV[t]+1)
    for(i in 1:nChainsV[t])
    {
      errAc[2*i]<-last_scM[sum(nChainsV[1:(t-1)])+i,1]
      errAcc[2*i+1]<-accRate[sum(nChainsV[1:(t-1)])+i]
    }
    write.table(t(errAc), file=fileOutErrAcc, append=TRUE, sep="\t", row.names=F, col.names=F)
   }
  }

  
  
  if(itn==numbIt)
   {
     write.table(t(numbSwapsA), file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)
     for(s in 1:numb_Sim)
     {
       write.table(t(permArray[[s]]), file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)
       write.table(t(nUp[[s]]), file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)
       write.table(t(nDown[[s]]), file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)
       write.table(t(upDown[[s]]), file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)
     }
   }


}#end for (itn in 1:...)

for(s in 1:numb_Sim)
{
  for (t in 1:(nChainsV[s]-1))
  {
    logSwapA[[s]][t]<-logSwapA[[s]][t]/numbSwapsA[endSeq[s]]
    if((nUp[[s]][t]+nDown[[s]][t])>0)
     {
       diffA[[s]][t]<-nUp[[s]][t]/(nUp[[s]][t]+nDown[[s]][t])
     }
  }
}

write.table("Log Swap rates and diffusion function are", file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)
for(s in 1:numb_Sim)
{
  write.table(t(logSwapA[[s]]), file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)
  write.table(t(diffA[[s]]), file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)
}


   
stopCluster(cl)
##############################################################




