#' Circuit optimization with Metropolis-Hastings (MH) algorithm
#' @param network_top: topology of the full network
#' @param data: processed gene expression matrix 
#' @param clusterRef: cluster indices of all models
#' @param cenMedRef: cluster centers
#' @param cutOffM: cluster radii
#' @param gene_list: gene clustering output 
#' @param init_top: initial circuit topology
#' @param output: a string of file prefix for saving results ("Results")
#' @param nRepeat: number of repeats of RACIPE simulations for each new circuit topology (5)
#'         A new circuit is simulated by RACIPE nRepeat times for robust score evaluation; 
#'         The scores will then be saved and used in future iterations, when the circuits are sampled again.
#' @param nIter: number of iterations for each simulation (1400)
#' @param modelsCGr: number of RACIPE models to be simulated (10000)
#' @param tempM: temperature for MH (60)
#' @return df: topology of the optimized CG circuit
#' @export
#' @import doParallel
opt_MH <-function(network_top, data, clusterRef, cenMedRef, cutOffM, gene_list, init_top, 
                  output = "Results", nRepeat= 5, nIter = 1400, modelsCGr = 10000, 
                  tempM=60){
#  require(doParallel)
  dataRow = t(data)
  fileAllSamp<-paste0(output, "_tops_allSampled.txt")
  fileAcc<-paste0(output,"_acc.txt")
    
  adjProb<-adjMatCGrProb(network_top, gene_list)
  
  ######################################################################
  numbCGRnodes<-length(gene_list)
  
  bNod<-find_inout(adjProb, numbCGRnodes)
  inNo<-bNod[[1]]
  ouNo<-bNod[[2]]
  ####
  
  list_sampled<-list()
  list_sampled[[1]]<-init_top
  list_accepted<-list_sampled
      
  red_lastSc<-NULL 
  for(s in 1:nRepeat){
      resPar<-inNWsimilarityRefCGr(dataRow, clusterRef, cenMedRef, cutOffM, gene_list, inNo, init_top, modelsCGr) 
      red_lastSc<-rbind(red_lastSc,resPar)
  }
  rMean<-apply(red_lastSc,2,mean)
  scSampV<-rMean[1]
  thrownOut<-rMean[2]
  write.table(t(c(init_top,rMean)), file=fileAcc, append=TRUE, sep="\t", row.names=F, col.names=F)
  write.table(t(c(init_top,rMean)), file=fileAllSamp, append=TRUE, sep="\t", row.names=F, col.names=F)
  
  last_score<-rMean
   
  for(itn in 1:nIter)
  {
    
    nSamp<-length(list_sampled)    
    old_top<-list_accepted[[itn]]
    new_top<-inNrepsample_topology(adjProb, numbCGRnodes, inNo, ouNo, old_top)
    if(sum(new_top)==0)
    {
        new_top<-old_top
    }
    indF<-find_topology(list_sampled, new_top, 1, nSamp)
                
    if(indF==0)
    {
      red_lastSc<-NULL 
      for(s in 1:nRepeat)
      {
       resPar<-inNWsimilarityRefCGr(dataRow, clusterRef, cenMedRef, cutOffM, gene_list, inNo, new_top, modelsCGr) 
       red_lastSc<-rbind(red_lastSc,resPar)
      }
   
      rMean<-apply(red_lastSc,2,mean)
      scSampV<-c(scSampV,rMean[1])
      thrownOut<-c(thrownOut,rMean[2])
      write.table(t(c(new_top,rMean)), file=fileAllSamp, append=TRUE, sep="\t", row.names=F, col.names=F)
      nSamp<-nSamp+1
      list_sampled[[nSamp]]<-new_top
      new_score<-rMean
    }
      
    if(indF!=0)
    {
       new_score<-c(scSampV[indF],thrownOut[indF])
    }    
              
      
    rnUn<-runif(1)
    
    if (rnUn< exp((last_score[1]-new_score[1])/tempM)){
        list_accepted[[itn+1]]<-new_top
        last_score<-new_score  
      }
    else{
        list_accepted[[itn+1]]<-list_accepted[[itn]]
      }   
      
    xout<-c(list_accepted[[itn+1]], last_score)
    write.table(t(xout), file=fileAcc, append=TRUE, sep="\t", row.names=F, col.names=F)
      
  }#end for(itn)
     
allS<-read.table(fileAllSamp,sep="\t")
allS<-as.matrix(allS)
colSc<-numbCGRnodes*numbCGRnodes+1
indMin<-which(allS[,colSc]==min(allS[,colSc]))[1]
bestTop<-allS[indMin,1:(colSc-1)]
df<-convAdjTop(bestTop,numbCGRnodes,1:numbCGRnodes)

return(df)
}

