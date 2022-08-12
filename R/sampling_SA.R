#' Circuit optimization with Simulated Annealing (SA) algorithm
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
#' @param modelsCGr: number of RACIPE models to be simulated (10000)
#' @param maxT: maximum/initial temperature in SA (150)
#' @param threshT: a second temperature in SA, below which SA has a slower temperature decaying rate (40)
#' @param decayRate1: 1st temperature decaying rate (geometrically decaying) (0.8)
#' @param decayRate2: 2nd temperature decaying rate (0.6), until temperature = 1 (current implementation)
#' @param iter_per_temp: number of iterations for each fixed temperature (100)
#' @return df: topology of the optimized CG circuit
#' @export
opt_SA <-function(network_top, data, clusterRef, cenMedRef, cutOffM, gene_list, init_top, 
                 output = "Results", nRepeat= 5, modelsCGr = 10000, 
                 maxT=150, decayRate1=0.8, decayRate2=0.6, threshT=40,iter_per_temp=100){
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
  current_temp<-maxT
  
  while(current_temp>1){
    if (current_temp>threshT){
      decay_rate<-decayRate1
    }else{
      decay_rate<-decayRate2
    }

    numb_iter<-0
    while(numb_iter<iter_per_temp){
    
    nSamp<-length(list_sampled)
    nAcc<-length(list_accepted)    
    old_top<-list_accepted[[nAcc]]
    new_top<-inNrepsample_topology(adjProb, numbCGRnodes, inNo, ouNo, old_top)
    if(sum(new_top)==0){
        new_top<-old_top
    }
    indF<-find_topology(list_sampled, new_top, 1, nSamp)
                
    if(indF==0){
      red_lastSc<-NULL 
      for(s in 1:nRepeat){
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
      
    if(indF!=0){
       new_score<-c(scSampV[indF],thrownOut[indF])
    }    
      
    rnUn<-runif(1)
    
    if (rnUn< exp((last_score[1]-new_score[1])/current_temp)){
        list_accepted[[nAcc+1]]<-new_top
        last_score<-new_score  
      }
    else{
        list_accepted[[nAcc+1]]<-list_accepted[[nAcc]]
      }   
      
    xout<-c(list_accepted[[nAcc+1]], last_score)
    write.table(t(xout), file=fileAcc, append=TRUE, sep="\t", row.names=F, col.names=F)
    
    numb_iter<-numb_iter+1
  }#end while(numb_iter<iter_per_temp)
  
  current_temp<-decay_rate*current_temp

}#end while(current_temp>1)

    
allS<-read.table(fileAllSamp,sep="\t")
allS<-as.matrix(allS)
colSc<-numbCGRnodes*numbCGRnodes+1
indMin<-which(allS[,colSc]==min(allS[,colSc]))[1]
bestTop<-allS[indMin,1:(colSc-1)]
df<-convAdjTop(bestTop,numbCGRnodes,1:numbCGRnodes)

return(df)
}

