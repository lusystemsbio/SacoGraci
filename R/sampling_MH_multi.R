#' Circuit optimization with Metropolis-Hastings (MH) algorithm (multiple threads)
#' @param network_top topology of the full network
#' @param data processed gene expression matrix 
#' @param clusterRef cluster indices of all models
#' @param cenMedRef cluster centers
#' @param cutOffM cluster radii
#' @param gene_list gene clustering output 
#' @param inTopsM a list of all initial circuit topologies
#' @param output a string of file prefix for saving results ("Results")
#' @param nRepeat number of repeats of RACIPE simulations for each new circuit topology (5)
#'         A new circuit is simulated by RACIPE nRepeat times for robust score evaluation; 
#'         The scores will then be saved and used in future iterations, when the circuits are sampled again.
#' @param nIter number of iterations for each simulation (1400)
#' @param modelsCGr number of RACIPE models to be simulated (10000)
#' @param tempM temperature for MH (60)
#' @param numbThr number of requested threads for HPC (40) 
#' @param nSim  number of parallel simulations (20)
#' @return df: topology of the optimized CG circuit
#' @export
#' @import doParallel
#' @import parallel
#' @import foreach
#' @importFrom stats runif
#' @importFrom utils read.table write.table
opt_MH_multi <-function(network_top, data, clusterRef, cenMedRef, cutOffM, gene_list, inTopsM, 
                        output = "Results", nRepeat= 5, nIter = 1400, modelsCGr = 10000, 
                        tempM = 60, numbThr = 40, nSim = 20){
  
  dataRow = t(data)
  fileAllSamp<-paste0(output, "_tops_allSampled.txt")
  fileOutStuff<-paste0(output,"_Stuff.txt")
  
  flnV<-NULL
  sampV<-NULL
  
  for(i in 1:nSim){
    nname<-paste(paste(output,"acc", as.character(i), sep="_"), "txt", sep=".")
    flnV<-c(flnV,nname)
  }
  
  adjProb<-adjMatCGrProb(network_top, gene_list)
  
  ######################################################################
  numbCGRnodes<-length(gene_list)
  
  bNod<-find_inout(adjProb, numbCGRnodes)
  inNo<-bNod[[1]]
  ouNo<-bNod[[2]]
  ####
  redTop<-matrix(inTopsM[1:nSim,],nrow=nSim)
  redIn<-redTop[!duplicated(redTop),]
  if(length(redIn)==(numbCGRnodes*numbCGRnodes)){
    redIn<-matrix(redIn,nrow=1)
  }
  
  uniqueIn<-dim(redIn)[1] #number of unique initial topologies
  
  list_sampled<-list()
  for(i in 1:uniqueIn){
    list_sampled[[i]]<-redIn[i,]
  }
  
  list_accepted<-list()
  for(i in 1:nSim){
    list_accepted[[i]]<-list()
    list_accepted[[i]][[1]]<-inTopsM[i,]
  }
  
  thrownOut<-NULL
  scSampV<-NULL
  
  cl<-makeCluster(numbThr)#takes 52 minutes to run
  registerDoParallel(cl)
  
  #clusterExport(cl, varlist=c('dataRow', 'clusterRef','cenMedRef', 'cutOffM', 'gene_list', 'flnV'))
  
  posM<-NULL 
  for(t in 1:uniqueIn){
    partM<-matrix(rep(redIn[t,],nRepeat),nrow=nRepeat,byrow=T)
    posM<-rbind(posM, partM)
  }
  
  nSimul<-nRepeat*uniqueIn
  nrep<-ceiling(nSimul/numbThr)
  
  red_lastSc<-NULL
  if(nrep>1){
    for(s in 1:(nrep-1)){
      resPar<-foreach(j=1:numbThr, .combine=c, .packages=c('sRACIPE','igraph')) %dopar% 
        inNWsimilarityRefCGr(dataRow, clusterRef, cenMedRef, cutOffM, gene_list, inNo, posM[(s-1)*numbThr+j,], modelsCGr) 
      
      vecRes<-rep(0,numbThr*2)
      for(m in 1:(numbThr*2)){
        vecRes[m]<-resPar[[m]]
      }
      matRes<-matrix(vecRes,nrow=numbThr,ncol=2,byrow=TRUE)
      red_lastSc<-rbind(red_lastSc,matRes)
    }
    if(nSimul>=((nrep-1)*numbThr+1)){
      resPar2<-foreach(j=((nrep-1)*numbThr+1):nSimul, .combine=c, .packages=c('sRACIPE','igraph')) %dopar% 
        inNWsimilarityRefCGr(dataRow, clusterRef, cenMedRef, cutOffM, gene_list, inNo, posM[j,], modelsCGr) 
      
      vecRes<-rep(0,(nSimul-(nrep-1)*numbThr)*2)
      for(m in 1:((nSimul-(nrep-1)*numbThr)*2)){
        vecRes[m]<-resPar2[[m]]
      }
      matRes<-matrix(vecRes,nrow=nSimul-(nrep-1)*numbThr,ncol=2,byrow=TRUE)
      red_lastSc<-rbind(red_lastSc, matRes)
    }
  }else {
    resPar<-foreach(j=1:nSimul, .combine=c, .packages=c('sRACIPE','igraph')) %dopar% 
      inNWsimilarityRefCGr(dataRow, clusterRef, cenMedRef, cutOffM, gene_list, inNo, posM[j,], modelsCGr) 
    
    vecRes<-rep(0,nSimul*2)
    for(m in 1:(nSimul*2)){
      vecRes[m]<-resPar[[m]]
    }
    matRes<-matrix(vecRes,nrow=nSimul,ncol=2,byrow=TRUE)
    red_lastSc<-rbind(red_lastSc,matRes)
  }
  
  if(nRepeat > 1){
    for(t in 1:uniqueIn){
      pMa<-red_lastSc[(nRepeat*t-nRepeat+1):(nRepeat*t),]
      scSampV<-c(scSampV, mean(pMa[,1]))
      thrownOut<-c(thrownOut, mean(pMa[,2]))
    }
  }else{
    for(t in 1:uniqueIn){
      pMa<-red_lastSc[t,]
      scSampV<-c(scSampV,pMa[1])
      thrownOut<-c(thrownOut, pMa[2])
    }
  }
  
  last_scM<-NULL
  
  for(i in 1:nSim){
    indF<-find_topology(list_sampled, inTopsM[i,], 1, uniqueIn)
    last_scM<-rbind(last_scM, c(scSampV[indF],thrownOut[indF]))
    
    xout<-c(inTopsM[i,], last_scM[i,])
    write.table(t(xout), file=flnV[i], append=TRUE, sep="\t", row.names=F, col.names=F)
  }
  
  if(uniqueIn>1){
    xouts<-cbind(rep(0,uniqueIn),cbind(cbind(redIn,scSampV),thrownOut))
  }
  if(uniqueIn==1){
    xouts<-c(0,as.vector(redIn),scSampV,thrownOut)
    xouts<-matrix(xouts,nrow=1)
  }
  write.table(xouts, file=fileAllSamp, append=TRUE, sep="\t", row.names=F, col.names=F)
  
  temp_value<-tempM
  
  for(itn in 1:nIter){
    newTop<-NULL
    indNew<-NULL
    newPosM<-NULL
    nSamp<-length(list_sampled)
    new_scM<-last_scM  
    
    for(i in 1:nSim){
      aiter<-length(list_accepted[[i]])
      old_topol<-list_accepted[[i]][[aiter]]
      new_top<-inNrepsample_topology(adjProb, numbCGRnodes, inNo, ouNo, old_topol)
      indF<-find_topology(list_sampled, new_top, 1, nSamp)
      if(sum(new_top)>0){
        if(indF==0){
          newTop<-rbind(newTop,new_top)
          indNew<-c(indNew,i)
          
        }
        
        newPosM<-rbind(newPosM,new_top)
        
      }#end if (sum(new_top)>0)
      
      if(sum(new_top)==0){
        new_scM[i,]<-last_scM[i,]
        newPosM<-rbind(newPosM, old_topol)
      }
      
    }#end for(i in 1:nSim)
    
    nUpd<-length(indNew)
    
    if(nUpd>0){
      write.table(t(c(itn,indNew)),file=fileOutStuff,append=TRUE, sep="\t", row.names=F, col.names=F)
      
      redM<-newTop[!duplicated(newTop),]
      realSim<-as.integer(length(redM)/(numbCGRnodes*numbCGRnodes))
      #realSim-number of unique new sampled topologies
      if(realSim==1){
        redM<-matrix(redM,nrow=1,byrow=T)
      }
      
      
      posM<-NULL 
      for(t in 1:realSim){
        partM<-matrix(rep(redM[t,],nRepeat),nrow=nRepeat,ncol=dim(redM)[2],byrow=T)
        posM<-rbind(posM, partM)
      }
      
      nSimul<-nRepeat*realSim
      nrep<-ceiling(nSimul/numbThr)
      
      red_newSc<-NULL
      if(nrep>1){
        for(s in 1:(nrep-1)){
          resPar<-foreach(j=1:numbThr, .combine=c, .packages=c('sRACIPE','igraph')) %dopar% 
            inNWsimilarityRefCGr(dataRow, clusterRef, cenMedRef, cutOffM, gene_list, inNo, posM[(s-1)*numbThr+j,], modelsCGr) 
          
          vecRes<-rep(0,numbThr*2)
          for(m in 1:(numbThr*2)){
            vecRes[m]<-resPar[[m]]
          }
          matRes<-matrix(vecRes,nrow=numbThr,ncol=2,byrow=TRUE)
          red_newSc<-rbind(red_newSc,matRes)
        }
        if(nSimul>=((nrep-1)*numbThr+1)){
          resPar2<-foreach(j=((nrep-1)*numbThr+1):nSimul, .combine=c, .packages=c('sRACIPE','igraph')) %dopar% 
            inNWsimilarityRefCGr(dataRow, clusterRef, cenMedRef, cutOffM, gene_list, inNo, posM[j,], modelsCGr) 
          
          vecRes<-rep(0,(nSimul-(nrep-1)*numbThr)*2)
          for(m in 1:((nSimul-(nrep-1)*numbThr)*2)){
            vecRes[m]<-resPar2[[m]]
          }
          matRes<-matrix(vecRes,nrow=nSimul-(nrep-1)*numbThr,ncol=2,byrow=TRUE)
          red_newSc<-rbind(red_newSc, matRes)
        }
      }else{
        resPar<-foreach(j=1:nSimul, .combine=c, .packages=c('sRACIPE','igraph')) %dopar% 
          inNWsimilarityRefCGr(dataRow, clusterRef, cenMedRef, cutOffM, gene_list, inNo, posM[j,], modelsCGr) 
        
        vecRes<-rep(0,nSimul*2)
        for(m in 1:(nSimul*2)){
          vecRes[m]<-resPar[[m]]
        }
        matRes<-matrix(vecRes,nrow=nSimul,ncol=2,byrow=TRUE)
        red_newSc<-rbind(red_newSc,matRes)
      }
      
      if(nRepeat > 1){
        for(t in 1:realSim){
          pMa<-red_newSc[(nRepeat*t-nRepeat+1):(nRepeat*t),]
          scSampV<-c(scSampV, mean(pMa[,1]))
          thrownOut<-c(thrownOut, mean(pMa[,2]))
        }
      }else{
        for(t in 1:realSim){
          pMa<-red_newSc[t,]
          scSampV<-c(scSampV,pMa[1])
          thrownOut<-c(thrownOut, pMa[2])
        }
      }
      
      if(realSim>1){
        xouts<-cbind(rep(itn,realSim),cbind(cbind(redM,scSampV[(nSamp+1):(nSamp+realSim)]),thrownOut[(nSamp+1):(nSamp+realSim)]))
      }
      if(realSim==1){
        xouts<-c(itn,as.vector(redM),scSampV[nSamp+1],thrownOut[nSamp+1])
        xouts<-matrix(xouts,nrow=1)
      }
      
      write.table(xouts, file=fileAllSamp, append=TRUE, sep="\t", row.names=F, col.names=F)
      
      for(t in 1:realSim){
        nSamp<-nSamp+1
        list_sampled[[nSamp]]<-redM[t,]
      }
      
      new_scInM<-NULL
      
      for(i in 1:nUpd){
        indF<-find_topology(list_sampled, newTop[i,], 1, nSamp)
        new_scInM<-rbind(new_scInM, c(scSampV[indF],thrownOut[indF]))
      }
      
      new_scM[indNew,]<-new_scInM  
      
    }#end if(nUpd>0)
    
    nnew<-setdiff(1:nSim,indNew)
    if(length(nnew)>0){
      for(t in nnew){
        indF<-find_topology(list_sampled, newPosM[t,], 1, nSamp)
        new_scM[t,1]<-scSampV[indF]
        new_scM[t,2]<-0
      }
    }
    
    for(i in 1:nSim){  
      rnUn<-runif(1)
      n<-length(list_accepted[[i]])
      if (rnUn< exp((last_scM[i,1]-new_scM[i,1])/temp_value)){
        last_scM[i,]<-new_scM[i,]
        list_accepted[[i]][[n+1]]<-newPosM[i,]  
      }
      else{
        list_accepted[[i]][[n+1]]<-list_accepted[[i]][[n]]
      }   
      
      xout<-c(list_accepted[[i]][[n+1]], last_scM[i,])
      write.table(t(xout), file=flnV[i], append=TRUE, sep="\t", row.names=F, col.names=F)
      
    }
    
  }#end for(itn)
  
  stopCluster(cl)
  
  allS<-read.table(fileAllSamp,sep="\t")
  allS<-as.matrix(allS)
  colSc<-numbCGRnodes*numbCGRnodes+2
  indMin<-which(allS[,colSc]==min(allS[,colSc]))[1]
  bestTop<-allS[indMin,2:(colSc-1)]
  df<-convAdjTop(bestTop,numbCGRnodes,1:numbCGRnodes)
  
  return(df)
}

