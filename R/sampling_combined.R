#### Circuit optimization (combined MH, MA & TE)
#network_top: topology of the full network
# data: processed gene expression matrix
# clusterRef: cluster indices of all models
# cenMedRef: cluster centers
# cutOffM: cluster radii
# gene_list: gene clustering output 
# inTopsM: a list of all initial circuit topologies
# output: a string of file prefix for saving results ("Results")
# nRepeat: number of repeats of RACIPE simulations for each new circuit topology (5)
#          A new circuit is simulated by RACIPE nRepeat times for robust score evaluations; 
#          The scores will then be saved and used in future iterations, when the circuits are sampled again.
# nIter: number of iterations for each simulation 
# modelsCGr: number of RACIPE models to be simulated (10000)
# numbThr: number of HPC threads
# nSim: number of simultaneous simulations 
# (output): (in results files)
opt_combined <- function(network_top, data, clusterRef, cenMedRef, cutOffM, gene_list, inTopsM,
                         output = "Results", nRepeat= 5, nIter = 1400, modelsCGr = 10000, 
                         numbThr = 40, nSim = 20){
  require(doParallel)

  dataRow = t(data)
  fileAllSamp<-paste0(output, "_tops_allSampled.txt")
  fileOutStuff<-paste0(output,"_Stuff.txt")
  fileSampUnsc<-paste0(output,"_tops_SampUnsc.txt")
  
  flnV<-NULL
  sampV<-NULL
  
  for(i in 1:nSim){
    nname<-paste(paste(output,"acc", as.character(i), sep="_"), "txt", sep=".")
    flnV<-c(flnV,nname)
  }
  
  
  adjProb<-adjMatCGrProb(top_ex4, gene_list)
  
  ######################################################################
  numbCGRnodes<-length(gene_list)
  
  bNod<-find_inout(adjProb, numbCGRnodes)
  inNo<-bNod[[1]]
  ouNo<-bNod[[2]]
  ####
  
  redIn<-inTopsM[!duplicated(inTopsM),]
  realIn<-nSim
  
  unscSamp<-cbind(cbind(rep(0,nSim),1:nSim), inTopsM[1:nSim,])
  write.table(unscSamp, file=fileSampUnsc, append=TRUE, sep="\t", row.names=F, col.names=F)
  
  list_sampled<-list()
  for(i in 1:realIn){
    list_sampled[[i]]<-redIn[i,]
  }
  
  list_accepted<-list()
  for(i in 1:nSim){
    list_accepted[[i]]<-list()
    list_accepted[[i]][[1]]<-inTopsM[i,]
  }
  
  thrownOut<-NULL
  scSampV<-NULL
  
  cl<-makeCluster(numbThr)
  registerDoParallel(cl)
  
# clusterExport(cl, varlist=c('data', 'clusterRef', 'cenMedRef', 'cutOffM', 'gene_list', 'flnV'))
  
  posM<-NULL 
  for(t in 1:realIn){
    partM<-matrix(rep(redIn[t,],nRepeat),nrow=nRepeat,ncol=dim(redIn)[2],byrow=T)
    posM<-rbind(posM, partM)
  }
  
  nSimul<-nRepeat*realIn
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
    if(nSimul>((nrep-1)*numbThr+1)){
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
    for(t in 1:realIn){
      pMa<-red_lastSc[(nRepeat*t-nRepeat+1):(nRepeat*t),]
      scSampV<-c(scSampV, mean(pMa[,1]))
      thrownOut<-c(thrownOut, pMa[nRepeat,2])
    }
  }else{
    for(t in 1:realIn){
      pMa<-red_lastSc[t,]
      scSampV<-c(scSampV,pMa[1])
      thrownOut<-c(thrownOut, pMa[2])
    }
  }
  
  last_scM<-NULL
  
  for(i in 1:nSim){
    indF<-find_topology(list_sampled, inTopsM[i,], 1, realIn)
    last_scM<-rbind(last_scM, c(scSampV[indF],thrownOut[indF]))
    
    xout<-c(inTopsM[i,], last_scM[i,])
    write.table(t(xout), file=flnV[i], append=TRUE, sep="\t", row.names=F, col.names=F)
  }
  
  xouts<-cbind(rep(0,realIn),cbind(cbind(redIn,scSampV),thrownOut))
  write.table(xouts, file=fileAllSamp, append=TRUE, sep="\t", row.names=F, col.names=F)
  
  temp_value<-60
  
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
      if(realSim==1){
        redM<-matrix(redM,nrow=1,ncol=numbCGRnodes*numbCGRnodes,byrow=T)
      }
      
      
      unscSamp<-cbind(rep(itn,realSim), redM)
      write.table(unscSamp, file=fileSampUnsc, append=TRUE, sep="\t", row.names=F, col.names=F)
      
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
        if(nSimul>((nrep-1)*numbThr+1)){
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
        resPar<-foreach(j=1:nUpd, .combine=c, .packages=c('sRACIPE','igraph')) %dopar% 
          inNWsimilarityRefCGr(dataRow, clusterRef, cenMedRef, cutOffM, gene_list, inNo, posM[j,], modelsCGr) 
        
        vecRes<-rep(0,nUpd*2)
        for(m in 1:(nUpd*2)){
          vecRes[m]<-resPar[[m]]
        }
        matRes<-matrix(vecRes,nrow=nUpd,ncol=2,byrow=TRUE)
        red_newSc<-rbind(red_newSc,matRes)
      }
      
      if(nRepeat > 1){
        for(t in 1:realIn){
          pMa<-red_newSc[(nRepeat*t-nRepeat+1):(nRepeat*t),]
          scSampV<-c(scSampV, mean(pMa[,1]))
          thrownOut<-c(thrownOut, pMa[nRepeat,2])
        }
      }else{
        for(t in 1:nUpd){
          pMa<-red_newSc[t,]
          scSampV<-c(scSampV,pMa[1])
          thrownOut<-c(thrownOut, pMa[2])
        }
      }
      
      xouts<-cbind(rep(itn,realSim),cbind(cbind(redM,scSampV[(nSamp+1):(nSamp+realSim)]),thrownOut[(nSamp+1):(nSamp+realSim)]))
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
}

