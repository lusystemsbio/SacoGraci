#' @title Identify input and output nodes from a circuit topology
#' @description Input nodes: nodes with only outward edges; Output nodes: nodes with only inward edges
#' @param adjMatrPr: Adjacency matrix with probabilities of CG-circuits derived from the full network
#' @param num_cgnodes: number of CG nodes (an output of gene clustering)
#' @return resNod: a list of input nodes (1) and output nodes (2)
#' @export
find_inout<-function(adjMatrPr, numb_cgnodes)
{
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

#' Mismatches between the sampled and reference circuit topologies
#' @param sam_top: Adjacency matrix of the sampled circuit topology
#' @param ref_top: Adjacency matrix of the reference circuit topology
#' @return : number of mismatching edges
#' @export
difLen<-function(samp_top, ref_top)
{
  return(length(which(samp_top!=ref_top)))
}

#' Deleting some edges from a circuit topology while trying to maintain a connected circuit
#' @param samp_top: the current circuit topology
#' @param numb_genes: number of genes of the circuit
#' @param numb_del: number of edges to be deleted
#' @param pos_del: indices of edges that can be deleted
#' @param inNodes: a list of input nodes, output from the function find_inout
#' @param list_top: list of sampled circuit topologies, we make sure they are not sampled again
#' @return adj_fin: the modified circuit topology
#' @export
#' @import igraph
delEdges_topNew<-function(samp_top, numb_genes, numb_del, pos_del, inNodes, list_top)
{
#  require(igraph)
  n<-length(list_top)
  npossDel<-length(pos_del)
  allDel<-NULL
  probV<-rep(1/npossDel, npossDel)
  cdata_top<-samp_top
  backTr<-sum(probV)
  
  while ((length(allDel)<numb_del)&& (backTr>0))
  {
    
    vconn<-1
    agrInN<-TRUE
    indFind<-0
    vcont<-((vconn>0)&& (agrInN==TRUE))&&(indFind==0)
    
    indExit<-length(which(probV!=0))
    countDel<-length(allDel)
    
    while(vcont &&(countDel<numb_del))
    {
      delOne<-sample(1:npossDel, size=1, prob=probV)
      probV[delOne]<-0
      allDel<-c(allDel, delOne)
      tnew_top<-samp_top[-pos_del[allDel],]
      grAdj<-graph_from_data_frame(tnew_top, directed = FALSE, vertices = NULL)
      vconn<-vertex_connectivity(grAdj, checks = TRUE)
      tnew_adj<-convTopAdj(tnew_top, numb_genes,1:numb_genes) 
      adj_Matr<-matrix(tnew_adj, nrow=numb_genes, ncol=numb_genes, byrow=TRUE)    
      
      sumInc<-apply(adj_Matr,2,sum)
      indInN<-which(sumInc==0)
      if(length(indInN)==0)
      {
        indInN<-0
      }
      
      agrInN<-setequal(inNodes, indInN)
      
      if( length(allDel)==(numb_del))
      {
        indFind<-find_topology(list_top, tnew_adj, 1, n)
      }
      
      vcont<-((vconn>0)&& (agrInN==TRUE))&&(indFind==0)
      
      countDel<-countDel+1
    }     
    
    if (vcont==FALSE)
    {
      
      laDel<-length(allDel)
      if (laDel>0)
      { 
        allDel<-allDel[-laDel]       
      }
      
    }
    
    backTr<-sum(probV)
  }
  
  if (length(allDel)==numb_del)
  {
    fin_top<-samp_top[-pos_del[allDel],]
    adj_fin<-convTopAdj(fin_top, numb_genes, 1:numb_genes)
    return(adj_fin)
  }else{
    return(rep(0, numb_genes*numb_genes))
  }
  
}

#' @title Making sign changes to some edges of a circuit topology 
#' @description this process will be tried 3 times if no new circuit topology is sampled
#' @param samp_top: the current circuit topology
#' @param numb_genes: number of genes of the circuit
#' @param numb_changes: number of edges with sign changes (>0)
#' @param pos_changes: indices of edges that can have sign changes
#' @param inNodes: a list of input nodes, output from the function find_inout
#' @param list_top: list of sampled circuit topologies, we make sure they are not sampled again
#' @return cTop: the modified circuit topology
#' @export
signChanges_top<-function(samp_top, numb_genes, numb_changes, pos_changes, inNodes, list_top)
{
  npossCh<-length(pos_changes)
  
  m<-1
  vcont<-TRUE
  
  while(vcont && (m<=3))
  {
    
    probCV<-rep(1/npossCh, npossCh)
    backTr<-numb_changes
    
    cTop<-samp_top
    chOne<-sample(1:npossCh, size=numb_changes)
    
    cTop[pos_changes[chOne]]<-3-cTop[pos_changes[chOne]]
    adj_Matr<-matrix(cTop, nrow=numb_genes, ncol=numb_genes, byrow=TRUE)  
    
    sumInc<-apply(adj_Matr,2,sum)
    indInN<-which(sumInc==0)
    if(length(indInN)==0)
    {
      indInN<-0
    }
    
    agrInN<-setequal(inNodes, indInN)
    
    indF<-find_topology(list_top, cTop, 1,length(list_top))
    
    vcont<-(agrInN==FALSE)||(indF!=0)
    probCV[chOne]<-0
    
    while(vcont && (backTr==numb_changes))
    {
      cTop[pos_changes[chOne[numb_changes]]]<-3-cTop[pos_changes[chOne[numb_changes]]]
      
      cprobCV<-probCV
      while (vcont && (sum(cprobCV)>0))
      {
        ccTop<-cTop
        singleCh<-sample(1:npossCh, size=1, prob=cprobCV)
        cprobCV[singleCh]<-0
        ccTop[pos_changes[singleCh]]<-3-ccTop[pos_changes[singleCh]]
        
        adj_Matr<-matrix(ccTop, nrow=numb_genes, ncol=numb_genes, byrow=TRUE)  
        sumInc<-apply(adj_Matr,2,sum)
        indInN<-which(sumInc==0)
        if(length(indInN)==0)
        {
          indInN<-0
        }
        
        agrInN<-setequal(inNodes, indInN)
        
        
        indF<-find_topology(list_top, ccTop, 1,length(list_top))
        
        vcont<-(agrInN==FALSE)||(indF!=0)
        
      }
      
      if (vcont==TRUE)
      {
        backTr<-backTr-1
      }
      else
      {
        cTop<-ccTop
      }
      
    }#end while(vcont && bactkTr==numb_changes)
    
    m<-m+1
    
  }#end while(vcont && m<=3)
  
  if (vcont==TRUE)
  {
    return(rep(0,numb_genes*numb_genes))
  }
  else
  {
    return (cTop)
  }
  
}#end function

#' Generating a series of randomly generated initial CG circuit topologies
#' @param circuit_top: the original circuit topology
#' @param gene_list: gene clustering output 
#' @param numbNewTop: number of new circuit topologies to be generated
#' @return resM: a list of randomly generated initial CG circuit topologies
#' @export
gaInitial_gen<-function(circuit_top, gene_list, numbNewTop)
{
  adjMProb<-adjMatCGrProb(circuit_top, gene_list)
  numb_cgnodes<-length(gene_list)
  
  inoN<-find_inout(adjMProb, numb_cgnodes)
  inNo<-inoN[[1]]
  
  #Find the densest topology
  #You'll delete edges from it and change edge signs
  
  inTop<-initialize_topology(adjMProb, numb_cgnodes)
  
  
  #Find the edges that allow sign changes and those that allow deletion
  
  bothSigns<-NULL
  canDel<-NULL
  for (i in 1:(numb_cgnodes*numb_cgnodes))
  {
    probE<-adjMProb[[i]]
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
  pos_del<-edgeCDel
  
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
      numb_del<-as.integer(raUn*maxDel)#this is number of deletions on inTop
      
      if (numb_del>0)
      { 
        cinTop<-convAdjTop(inTop,numb_cgnodes,1:numb_cgnodes)
        adjNew<-delEdges_topNew(cinTop,numb_cgnodes, numb_del, pos_del,inNo, list_topCG)
        pTop<-adjNew
      }
      
    }#end if(maxDel>0)
    
    #Now do sign changes, make sure you don't get duplicates
    
    if (sum(pTop)>0)
    {
      edgepTop<-which(pTop!=0)
      pos_changes<-intersect(edgepTop, bothSigns)   
      maxSign<-length(pos_changes)
      if (maxSign>0)
      {
        
        #Try to do sign changes 3 times, stop when you get new topology
        u<-1
        badCond<-TRUE
        while(badCond && (u<=3))
        {
          raUn<-runif(1)
          numb_changes<-as.integer(raUn*maxSign)#this is number of sign changes on inTop
          if (numb_changes > 0)
          { 
            ppTop<-signChanges_top(pTop, numb_cgnodes, numb_changes, pos_changes, inNo, list_topCG)
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

#' Builds the adjacency matrix with probabilities
#' @param data_top: topology of the full network
#' @param gene_list: gene clustering output 
#' @return adj_cgr: adjacency matrix with probabilities
#' @export
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

#' Build the most dense circuit topology (adj. matrix in a vector``)
#' @param adj_matrProb: adjacency matrix with probabilities
#' @param numb_cgnodes: the number of CG nodes
#' @return adj_Matr:  the adj. matrix of the most dense topology as a vector
#' @export
#' @import igraph
initialize_topology<-function(adj_matrProb, numb_cgnodes)
{
#  require(igraph)
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

#' Converting the adjacency matrix of a gene circuit to the topology file (Source/Target/Interaction Type)
#' @param adjMatL: adjacency matrix
#' @param numbG: the number of genes
#' @param colN: a vector of gene names
#' @return data_cgr_top:  a data frame of the circuit topology
#' @export
convAdjTop<-function(adjMatL, numbG, colN)
{
  #colN=vector of gene names=1:numbG by default
  adjMatLin<-unlist(adjMatL)
  nnEntr<-which(adjMatLin!=0)
  numb_conn_cgr<-length(nnEntr)
  
  cgr_top<-NULL
  
  for (i in 1:numb_conn_cgr){
    posN<-nnEntr[i]
    colposin<-posN%%(numbG)
    ###Stopped HERE
    if(colposin==0){
      colpos<-numbG
      ropos<-(posN%/%(numbG))
    }else{
      colpos<-colposin
      ropos<-(posN%/%(numbG))+1
    }
    
    cgr_top<-rbind(cgr_top, c(colN[ropos], colN[colpos], adjMatLin[posN]))  
    
  }
  
  colnames(cgr_top)<-c("Source", "Target", "Type")
  
  data_cgr_top<-as.data.frame(cgr_top)
  
  return(data_cgr_top)
  
}

#' Converting the topology file (Source/Target/Interaction Type) of a gene circuit to the adjacency matrix d
#' @param data_t:  a data frame of the circuit topology
#' @param numbG: the number of genes
#' @param colN: a vector of gene names
#' @return adjMatL: adjacency matrix
#' @export
convTopAdj<-function(data_t, numbG, colN)
{#data_t<-emtdata_top
  #numbG<-22
  #colN<-colnames(data_t)
  adjMatL<-rep(0, numbG*numbG)
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
    adjMatL[(vsou-1)*numbG+vtgt]<-vtype
  }
  
  return(adjMatL)
}

#' Check whether the circuit topology has been sampled or not
#' @param top_list:  a list of sampled topologies
#' @param samp_top: the current topology
#' @param minN: the range of circuit topologies in the list to be searched (minimum index value)
#' @param maxN: the range of circuit topologies in the list to be searched (maximum index value)
#' @return i: the index of the found topology in the list, or 0 if not found
#' @export
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

#' Find an element in a list
#' @param xli:the list
#' @param n: the length of the list
#' @param elem: the element to search
#' @return indli: the index of the found element in the list, or 0 if not found
#' @export
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

#' Circuit sampling
#' @param adj_matrPr: adjacency matrix of the edge probability 
#' @param numb_cgnodes: the nuymber of CG nodes
#' @param inNodes: a list of input nodes
#' @param outNodes: a list of output nodes
#' @param old_top: the last circuit topology
#' @return new_top: the new circuit topology
#' @export
#' @import igraph
inNrepsample_topology<-function(adj_matrPr, numb_cgnodes, inNodes, outNodes, old_top)
{
#  require(igraph)
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

#' Visualize network topology
#' @param tf_links: circuit topology
#' @param height: plot height ("300px")
#' @return network plot
#' @export
#' @import visNetwork
plot_network <- function(tf_links, height = "300px"){
#  require(visNetwork)
  topology=data.frame(as.matrix(tf_links), stringsAsFactors = F)

  node_list <- unique(c(topology[,1], topology[,2]))
  nodes <- data.frame(id = node_list, label = node_list, font.size =30, value=c(rep(1,length(node_list))))
  
  edge_col <- data.frame(c(1,2),c("blue","darkred"))
  colnames(edge_col) <- c("relation", "color")
  arrow_type <- data.frame(c(1,2),c("arrow","circle"))
  colnames(arrow_type) <- c("type", "color")
  edges <- data.frame(from =c(topology[,1]), to = c(topology[,2])
                      , arrows.to.type	=arrow_type$color[c(as.numeric(topology[,3]))]
                      , width = 3
                      , color = edge_col$color[c(as.numeric(topology[,3]))]
  )
  visNetwork(nodes, edges, height = height, width = "100%") %>%
    visEdges(arrows = "to") %>%
    visOptions(manipulation = TRUE) %>%
    visLayout(randomSeed = 123) %>%
    visPhysics(solver = "forceAtlas2Based", stabilization = FALSE)
}
