#' Circuit optimization with temperature tempting (TE) algorithm
#' @param network_top: topology of the full network
#' @param data: processed gene expression matrix 
#' @param clusterRef: cluster indices of all models
#' @param cenMedRef: cluster centers
#' @param cutOffM: cluster radii
#' @param gene_list: gene clustering output 
#' @param inTopsM: a list of all initial circuit topologies
#' @param output: a string of file prefix for saving results ("Results")
#' @param nRepeat: number of repeats of RACIPE simulations for each new circuit topology (5)
#'         A new circuit is simulated by RACIPE nRepeat times for robust score evaluation; 
#'         The scores will then be saved and used in future iterations, when the circuits are sampled again.
#' @param modelsCGr: number of RACIPE models to be simulated (10000)
#' @param numbThr: number of requested threads for HPC (40) 
#' @param temp_Array: temperatures for all replicas
#'        Default: a total of 24 replicas with temperatures:
#'        c(1,1.05,1.1,1.15,1.2,1.25,1.3,1.5,2.0,2.5,3.0,3.5,4.0,6,9,11,13,20,28,40,55,70,90,120)
#' @param iter_temp_add: the number of iterations (proposed swaps) during the temperature addition process
#'        Default: three runs of the procedure; (c(50,100,150)) 
#' @param numb_iter_extra: number of extra iterations after the temperature addition process (1100)
#' @param logAlpha: log of the target swap rate (log(0.4))
#' @return df: topology of the optimized CG circuit
#' @export
#' @import doParallel
#' @import parallel
#' @import foreach 
opt_TE<-function(network_top, data, clusterRef, cenMedRef, cutOffM, gene_list, inTopsM,
                 output = "Results", nRepeat = 5, modelsCGr = 10000,
                 numbThr = 2, temp_Array=c(1,1.05,1.1,1.15,1.2,1.25,1.3,1.5,2.0,2.5,3.0,
                                            3.5,4.0,6,9,11,13,20,28,40,55,70,90,120), 
                 iter_temp_add=c(50,100,150), numb_iter_extra=1100, logAlpha=log(0.4)){
  
  dataRow = t(data)
  fileOutTop<-paste0(output, "_tops.txt")
  fileAllSamp<-paste0(output, "_tops_allSampled.txt")
  fileOutErrAcc<-paste0(output,"_scores.txt")
  fileOutStuff<-paste0(output,"_stuff.txt")
  
  adjProb<-adjMatCGrProb(network_top, gene_list)
  
  numbCGRnodes<-length(gene_list)
  
  bNod<-find_inout(adjProb, numbCGRnodes)
  inNo<-bNod[[1]]
  ouNo<-bNod[[2]]
  
#  inTopsM<-gaInitial_gen(network_top, gene_list, 80)
  
  list_sampled<-list()
  scSampV<-NULL
  thrownOut<-NULL
  
  cl<-makeCluster(numbThr)#takes 52 minutes to run
  registerDoParallel(cl)
  
  nChains<-length(temp_Array)
  
  #####CALCULATE SCORES INITIAL TOP'S THEN WRITE THEM OUT
  
  posM<-as.matrix(inTopsM[1:nChains,])
  redIn<-posM[!duplicated(posM),]
  if(length(redIn)==(numbCGRnodes*numbCGRnodes)){
    redIn<-matrix(redIn,nrow=1)
  }
  
  uniqueIn<-dim(redIn)[1] #number of unique initial topologies
  
  for(i in 1:uniqueIn){
    list_sampled[[i]]<-redIn[i,]
  }
  
  calc_posM<-NULL 
  for(t in 1:uniqueIn){
    partM<-matrix(rep(redIn[t,],nRepeat),nrow=nRepeat,byrow=T)
    calc_posM<-rbind(calc_posM, partM)
  }
  
  nSimul<-nRepeat*uniqueIn
  nrep<-ceiling(nSimul/numbThr)
  
  red_lastSc<-NULL
  if(nrep>1){
    for(s in 1:(nrep-1)){
      resPar<-foreach(j=1:numbThr, .combine=c, .packages=c('sRACIPE','igraph')) %dopar% 
        inNWsimilarityRefCGr(dataRow, clusterRef, cenMedRef, cutOffM, gene_list, inNo, calc_posM[(s-1)*numbThr+j,], modelsCGr) 
      
      vecRes<-rep(0,numbThr*2)
      for(m in 1:(numbThr*2)){
        vecRes[m]<-resPar[[m]]
      }
      matRes<-matrix(vecRes,nrow=numbThr,ncol=2,byrow=TRUE)
      red_lastSc<-rbind(red_lastSc,matRes)
    }
    if(nSimul>=((nrep-1)*numbThr+1)){
      resPar2<-foreach(j=((nrep-1)*numbThr+1):nSimul, .combine=c, .packages=c('sRACIPE','igraph')) %dopar% 
        inNWsimilarityRefCGr(dataRow, clusterRef, cenMedRef, cutOffM, gene_list, inNo, calc_posM[j,], modelsCGr) 
      
      vecRes<-rep(0,(nSimul-(nrep-1)*numbThr)*2)
      for(m in 1:((nSimul-(nrep-1)*numbThr)*2)){
        vecRes[m]<-resPar2[[m]]
      }
      matRes<-matrix(vecRes,nrow=nSimul-(nrep-1)*numbThr,ncol=2,byrow=TRUE)
      red_lastSc<-rbind(red_lastSc, matRes)
    }
  }else {
    resPar<-foreach(j=1:nSimul, .combine=c, .packages=c('sRACIPE','igraph')) %dopar% 
      inNWsimilarityRefCGr(dataRow, clusterRef, cenMedRef, cutOffM, gene_list, inNo, calc_posM[j,], modelsCGr) 
    
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
  xout<-NULL
  for(i in 1:nChains){
    indF<-find_topology(list_sampled, inTopsM[i,], 1, uniqueIn)
    last_scM<-rbind(last_scM, c(scSampV[indF],thrownOut[indF]))
    
    xout<-rbind(xout,c(0, inTopsM[i,], last_scM[i,]))
  }
  
  matInd<-matrix(rep(0,2*nChains),nrow=nChains, ncol=2,byrow=T)
  matInd[,2]<-1:nChains
  xout<-cbind(matInd,xout)
  
  write.table(xout, file=fileOutTop, append=TRUE, sep="\t", row.names=F, col.names=F)
  
  
  if(uniqueIn>1){
    xouts<-cbind(rep(0,uniqueIn),cbind(cbind(redIn,scSampV),thrownOut))
  }
  if(uniqueIn==1){
    xouts<-c(0,as.vector(redIn),scSampV,thrownOut)
    xouts<-matrix(xouts,nrow=1)
  }
  write.table(xouts, file=fileAllSamp, append=TRUE, sep="\t", row.names=F, col.names=F)
  
  errAc<-rep(0, 2*nChains+1)
  for(i in 1:nChains){
    errAc[2*i]<-last_scM[i,1]
  }
  write.table(t(errAc), file=fileOutErrAcc, append=TRUE, sep="\t", row.names=F, col.names=F)
  
  
  ######START THE ITERATIONS
  #FOR ADDING TEMPERATURES
  
  for(itAdd in 2:length(iter_temp_add)){
    #itAdd<-1
    
    logSwapA<-rep(0,nChains-1)
    numbSwapsA<-rep(0,nChains+2)
    diffA<-rep(0,nChains)
    errAcc<-rep(0,2*nChains+1)
    
    accRate<-rep(0,nChains)
    nUp<-accRate
    nDown<-nUp
    upDown<-nUp
    permArray<-1:nChains
    
    nUp[1]<-1
    nDown[nChains]<-1
    upDown[1]<-1
    upDown[nChains]<--1
    
    swEven<-seq(from=nChains-1, to=1, by=-1)
    
    numbIt<-iter_temp_add[itAdd]
    
    ###START THE ITERATIONS FOR ADDING TEMPERATURES
    
    for(itn in 1:numbIt){  #itn<-1
      #update the topologies of each chain
      newTop<-NULL
      indNew<-NULL
      nSamp<-length(list_sampled)
      
      new_scM<-last_scM
      newPosM<-posM
      for(i in 1:nChains){
        new_top<-inNrepsample_topology(adjProb, numbCGRnodes, inNo, ouNo, posM[i,])
        indF<-find_topology(list_sampled, new_top, 1, nSamp)
        if(sum(new_top)>0){
          if(indF==0){
            newTop<-rbind(newTop,new_top)
            indNew<-c(indNew,i) 
          }
          
          newPosM[i,]<-new_top
          
        }#end if (sum(new_top)>0)
        
        if(sum(new_top)==0){
          newPosM[i,]<-posM[i,]
          new_scM[i,]<-last_scM[i,]
        }
        
      }#end for(i in 1:nChains)
      
      
      nUpd<-length(indNew)
      
      if(nUpd>0){
        
        write.table(t(c(itn,indNew)),file=fileOutStuff,append=TRUE, sep="\t", row.names=F, col.names=F)
        
        redM<-newTop[!duplicated(newTop),]
        realSim<-as.integer(length(redM)/(numbCGRnodes*numbCGRnodes))
        #realSim-number of unique new sampled topologies
        if(realSim==1){
          redM<-matrix(redM,nrow=1,byrow=T)
        }
        
        calc_posM<-NULL 
        for(t in 1:realSim){
          partM<-matrix(rep(redM[t,],nRepeat),nrow=nRepeat,byrow=T)
          calc_posM<-rbind(calc_posM, partM)
        }
        
        nSimul<-nRepeat*realSim
        nrep<-ceiling(nSimul/numbThr)
        
        red_newSc<-NULL
        if(nrep>1){
          for(s in 1:(nrep-1)){
            resPar<-foreach(j=1:numbThr, .combine=c, .packages=c('sRACIPE','igraph')) %dopar% 
              inNWsimilarityRefCGr(dataRow, clusterRef, cenMedRef, cutOffM, gene_list, inNo, calc_posM[(s-1)*numbThr+j,], modelsCGr) 
            
            vecRes<-rep(0,numbThr*2)
            for(m in 1:(numbThr*2)){
              vecRes[m]<-resPar[[m]]
            }
            matRes<-matrix(vecRes,nrow=numbThr,ncol=2,byrow=TRUE)
            red_newSc<-rbind(red_newSc,matRes)
          }
          if(nSimul>=((nrep-1)*numbThr+1)){
            resPar2<-foreach(j=((nrep-1)*numbThr+1):nSimul, .combine=c, .packages=c('sRACIPE','igraph')) %dopar% 
              inNWsimilarityRefCGr(dataRow, clusterRef, cenMedRef, cutOffM, gene_list, inNo, calc_posM[j,], modelsCGr) 
            
            vecRes<-rep(0,(nSimul-(nrep-1)*numbThr)*2)
            for(m in 1:((nSimul-(nrep-1)*numbThr)*2)){
              vecRes[m]<-resPar2[[m]]
            }
            matRes<-matrix(vecRes,nrow=nSimul-(nrep-1)*numbThr,ncol=2,byrow=TRUE)
            red_newSc<-rbind(red_newSc, matRes)
          }
        }
        if(nrep==1){
          resPar<-foreach(j=1:nSimul, .combine=c, .packages=c('sRACIPE','igraph')) %dopar% 
            inNWsimilarityRefCGr(dataRow, clusterRef, cenMedRef, cutOffM, gene_list, inNo, calc_posM[j,], modelsCGr) 
          
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
        
        if(nUpd>1){
          new_scM[indNew,]<-new_scInM  
        }
        if(nUpd==1){
          new_scM[indNew,]<-as.vector(new_scInM)
        }
        
      }#end if(nUpd>0)
      
      nnew<-setdiff(1:nChains,indNew)
      if(length(nnew)>0){
        for(t in nnew){
          indF<-find_topology(list_sampled, newPosM[t,], 1, nSamp)
          new_scM[t,1]<-scSampV[indF]
          new_scM[t,2]<-thrownOut[indF]
        }
      }
      
      #SEE IF WE ACCEPT NEW SCORES AND THEN DO SWAPS
      
      for(i in 1:nChains){
        if(last_scM[i,1]==new_scM[i,1]){
          accRate[i]<-accRate[i]-1
        }
        rnU<-runif(1)
        if(rnU<exp((last_scM[i,1]-new_scM[i,1])/temp_Array[i])){
          last_scM[i,]<-new_scM[i,]
          posM[i,]<-newPosM[i,]
          accRate[i]<-accRate[i]+1
        }
        
      }#end for(i in 1:nChains)
      
      #WRITE SCORES BEFORE THE SWAPS
      errAcc<-rep(itn,2*nChains+1)
      for(i in 1:nChains){
        errAcc[2*i]<-last_scM[i,1]
        errAcc[2*i+1]<-accRate[i]
      }
      
      write.table(t(errAcc), file=fileOutErrAcc, append=TRUE, sep="\t", row.names=F, col.names=F)
      
      #NOW DO THE SWAPS
      for(t in swEven){
        rnUn<-runif(1)
        swRate<-exp((last_scM[t+1,1]-last_scM[t,1])*(1/temp_Array[t+1]-1/temp_Array[t]))
        if(rnUn<swRate) { #if accept swap
          swErr<-last_scM[t+1,]
          last_scM[t+1,]<-last_scM[t,]
          last_scM[t,]<-swErr
          swPos<-posM[t+1,]
          posM[t+1,]<-posM[t,]
          posM[t,]<-swPos
          
          if((t==1)&&((permArray[2]==1)&&(upDown[1]==-1))){
            numbSwapsA[nChains+2]<-numbSwapsA[nChains+2]+1
          }    
          if(t==1){
            upDown[permArray[2]]<-1
          }
          if((t==nChains-1)&&((permArray[t]==1)&&(upDown[1]==1))){
            numbSwapsA[nChains+1]<-numbSwapsA[nChains+1]+1
          }     
          if(t==nChains-1){
            upDown[permArray[t]]<--1
          }
          if(upDown[permArray[t]]==1){
            nUp[t+1]<-nUp[t+1]+1
          }
          if(upDown[permArray[t+1]]==-1){
            nDown[t]<-nDown[t]+1
          }
          
          swPerm<-permArray[t]
          permArray[t]<-permArray[t+1]
          permArray[t+1]<-swPerm
          
          numbSwapsA[t]<-numbSwapsA[t]+1
        }
        
        if(swRate<1){
          logSwapA[t]<-logSwapA[t]+log(swRate)
        }
      }#end swap chains (1,2),(3,4)...
      
      numbSwapsA[nChains]<-numbSwapsA[nChains]+1 
      
      #### NOW WRITE THE POSITIONS SCORES ARRAYS
      
      matInd<-matrix(rep(itn,2*nChains),nrow=nChains, ncol=2,byrow=T)
      matInd[,2]<-1:nChains
      bi1<-cbind(matInd,posM)
      bi2<-cbind(bi1,last_scM)
      
      write.table(bi2, file=fileOutTop, append=TRUE, sep="\t", row.names=F, col.names=F)
      
      errAcc<-rep(itn,2*nChains+1)
      for(i in 1:nChains){
        errAcc[2*i]<-last_scM[i,1]
        errAcc[2*i+1]<-accRate[i]
      }
      
      write.table(t(errAcc), file=fileOutErrAcc, append=TRUE, sep="\t", row.names=F, col.names=F)
      write.table(t(c(itn,permArray)), file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)
      
      
      if(itn==numbIt){
        write.table(t(numbSwapsA), file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)
        write.table(t(permArray), file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)
        write.table(t(nUp), file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)
        write.table(t(nDown), file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)
        write.table(t(upDown), file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)
      }
      
    }#end for(itn in 1:numbIt)
    
    for (t in 1:(nChains-1)){
      logSwapA[t]<-logSwapA[t]/numbSwapsA[nChains]
      if((nUp[t]+nDown[t])>0){
        diffA[t]<-nUp[t]/(nUp[t]+nDown[t])
      }
    }
    
    write.table("Log Swap rates and diffusion function are", file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)
    write.table(t(logSwapA), file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)
    write.table(t(diffA), file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)
    
    ###ADD THE NEW TEMPERATURES
    rLogs<-rep(0, nChains-1)
    sumR<-0
    newTempA<-NULL
    for(i in 1:(nChains-1)){
      if(numbSwapsA[i]>0){
        rLogs[i]<-log(numbSwapsA[i]/numbSwapsA[nChains])/logAlpha
      }else{
        rLogs[i]<-logSwapA[i]/logAlpha
      }
    }
    
    write.table("RLogs are", file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)
    write.table(t(rLogs), file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)
    
    for(i in 1:(nChains-1)){    
      lrat<-as.integer(sqrt(rLogs[i]))
      sumR<-sumR+lrat
      
      for(m in 1:(lrat+1)){ 
        newT<-temp_Array[i]+(temp_Array[i+1]-temp_Array[i])*(m-1)/(lrat+1)
        newTempA<-c(newTempA, newT)
      }
    }
    
    newTempA<-c(newTempA, temp_Array[nChains])
    temp_Array<-newTempA
    sumR<-as.integer(sumR)
    
    write.table("SumR is", file = fileOutStuff, append = TRUE, sep = "\t", row.names = F, col.names = F)
    write.table(sumR, file = fileOutStuff, append = TRUE, sep = "\t", row.names = F, col.names = F)
    write.table("New temperature array is", file = fileOutStuff, append = TRUE, sep = "\t", row.names = F, col.names = F)
    write.table(t(newTempA), file = fileOutStuff, append = TRUE, sep = "\t", row.names = F, col.names = F)
    
    if (sumR>0){
      if(sumR==1){
        newTops<-matrix(inTopsM[nChains+1,], nrow=1, byrow=T)
      }
      if(sumR>1){
        #sseq<-seq(from=nChains+1, to=nChains+sumR, by=1)
        newTops<-inTopsM[(nChains+1):(nChains+sumR),]
      }
      
      posM<-rbind(posM,newTops)
      
      ##CALCULATE SCORES FOR THE ADDED TOPOLOGIES
      
      inredM<-newTops[!duplicated(newTops),]
      inrealSim<-as.integer(length(inredM)/(numbCGRnodes*numbCGRnodes))
      if(inrealSim==1){
        inredM<-matrix(inredM,nrow=1,byrow=T)
      }
      
      nTops<-NULL
      for(t in 1:inrealSim){
        indF<-find_topology(list_sampled, inredM[t,], 1, nSamp)
        if(indF==0){
          nTops<-rbind(nTops,inredM[t,])           
        }
      }
      
      realSim<-as.integer(length(nTops)/(numbCGRnodes*numbCGRnodes))
      if(realSim==1){
        nTops<-matrix(nTops,nrow=1,byrow=T)
      }
      
      if(realSim>0){
        calc_posM<-NULL 
        for(t in 1:realSim){
          partM<-matrix(rep(nTops[t,],nRepeat),nrow=nRepeat,byrow=T)
          calc_posM<-rbind(calc_posM, partM)
        }
        
        nSimul<-nRepeat*realSim
        nrep<-ceiling(nSimul/numbThr)
        
        red_newSc<-NULL
        if(nrep>1){
          for(s in 1:(nrep-1)){
            resPar<-foreach(j=1:numbThr, .combine=c, .packages=c('sRACIPE','igraph')) %dopar% 
              inNWsimilarityRefCGr(dataRow, clusterRef, cenMedRef, cutOffM, gene_list, inNo, calc_posM[(s-1)*numbThr+j,], modelsCGr) 
            
            vecRes<-rep(0,numbThr*2)
            for(m in 1:(numbThr*2)){
              vecRes[m]<-resPar[[m]]
            }
            matRes<-matrix(vecRes,nrow=numbThr,ncol=2,byrow=TRUE)
            red_newSc<-rbind(red_newSc,matRes)
          }
          if(nSimul>=((nrep-1)*numbThr+1)){
            resPar2<-foreach(j=((nrep-1)*numbThr+1):nSimul, .combine=c, .packages=c('sRACIPE','igraph')) %dopar% 
              inNWsimilarityRefCGr(dataRow, clusterRef, cenMedRef, cutOffM, gene_list, inNo, calc_posM[j,], modelsCGr) 
            
            vecRes<-rep(0,(nSimul-(nrep-1)*numbThr)*2)
            for(m in 1:((nSimul-(nrep-1)*numbThr)*2)){
              vecRes[m]<-resPar2[[m]]
            }
            matRes<-matrix(vecRes,nrow=nSimul-(nrep-1)*numbThr,ncol=2,byrow=TRUE)
            red_newSc<-rbind(red_newSc, matRes)
          }
        }else{
          resPar<-foreach(j=1:nSimul, .combine=c, .packages=c('sRACIPE','igraph')) %dopar% 
            inNWsimilarityRefCGr(dataRow, clusterRef, cenMedRef, cutOffM, gene_list, inNo, calc_posM[j,], modelsCGr) 
          
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
          xouts<-cbind(rep(itn,realSim),cbind(cbind(nTops,scSampV[(nSamp+1):(nSamp+realSim)]),thrownOut[(nSamp+1):(nSamp+realSim)]))
        }
        if(realSim==1){
          xouts<-c(itn,as.vector(nTops),scSampV[nSamp+1],thrownOut[nSamp+1])
          xouts<-matrix(xouts,nrow=1)
        }
        
        write.table(xouts, file=fileAllSamp, append=TRUE, sep="\t", row.names=F, col.names=F)
        
        for(t in 1:realSim){
          nSamp<-nSamp+1
          list_sampled[[nSamp]]<-nTops[t,]
        }
        
      }#end if(realSim>0)
      
      last_scM<-NULL
      ntop<-dim(posM)[1]
      for(t in 1:ntop)
      {
        indF<-find_topology(list_sampled, posM[t,], 1, nSamp)
        last_scM<-rbind(last_scM, c(scSampV[indF],thrownOut[indF]))
      }
      
      nChains<-nChains+sumR
      
    }#end(if sumR>0)
    
    write.table("New number of temperatures are", file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)
    write.table(t(nChains), file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)
    
    
  }#end for (itAdd in 1:length(iter_temp_add))
  
  #############################################
  #######DO THE LAST RUN
  
  for(itn in 1:numb_iter_extra){
    #update the topologies of each chain
    newTop<-NULL
    indNew<-NULL
    nSamp<-length(list_sampled)
    
    new_scM<-last_scM
    newPosM<-posM
    for(i in 1:nChains){
      new_top<-inNrepsample_topology(adjProb, numbCGRnodes, inNo, ouNo, posM[i,])
      indF<-find_topology(list_sampled, new_top, 1, nSamp)
      if(sum(new_top)>0){
        if(indF==0){
          newTop<-rbind(newTop,new_top)
          indNew<-c(indNew,i)   
        }
        newPosM[i,]<-new_top
        
      }#end if (sum(new_top)>0)
      
      if(sum(new_top)==0){
        newPosM[i,]<-posM[i,]
        new_scM[i,]<-last_scM[i,]
      }
      
    }#end for(i in 1:nChains)
    
    nUpd<-length(indNew)
    
    if(nUpd>0){
      
      write.table(t(c(itn,indNew)),file=fileOutStuff,append=TRUE, sep="\t", row.names=F, col.names=F)
      
      redM<-newTop[!duplicated(newTop),]
      realSim<-as.integer(length(redM)/(numbCGRnodes*numbCGRnodes))
      #realSim-number of unique new sampled topologies
      if(realSim==1){
        redM<-matrix(redM,nrow=1,byrow=T)
      }
      
      calc_posM<-NULL 
      for(t in 1:realSim){
        partM<-matrix(rep(redM[t,],nRepeat),nrow=nRepeat,byrow=T)
        calc_posM<-rbind(calc_posM, partM)
      }
      
      nSimul<-nRepeat*realSim
      nrep<-ceiling(nSimul/numbThr)
      
      red_newSc<-NULL
      if(nrep>1){
        for(s in 1:(nrep-1)){
          resPar<-foreach(j=1:numbThr, .combine=c, .packages=c('sRACIPE','igraph')) %dopar% 
            inNWsimilarityRefCGr(dataRow, clusterRef, cenMedRef, cutOffM, gene_list, inNo, calc_posM[(s-1)*numbThr+j,], modelsCGr) 
          
          vecRes<-rep(0,numbThr*2)
          for(m in 1:(numbThr*2)){
            vecRes[m]<-resPar[[m]]
          }
          matRes<-matrix(vecRes,nrow=numbThr,ncol=2,byrow=TRUE)
          red_newSc<-rbind(red_newSc,matRes)
        }
        if(nSimul>=((nrep-1)*numbThr+1)){
          resPar2<-foreach(j=((nrep-1)*numbThr+1):nSimul, .combine=c, .packages=c('sRACIPE','igraph')) %dopar% 
            inNWsimilarityRefCGr(dataRow, clusterRef, cenMedRef, cutOffM, gene_list, inNo, calc_posM[j,], modelsCGr) 
          
          vecRes<-rep(0,(nSimul-(nrep-1)*numbThr)*2)
          for(m in 1:((nSimul-(nrep-1)*numbThr)*2)){
            vecRes[m]<-resPar2[[m]]
          }
          matRes<-matrix(vecRes,nrow=nSimul-(nrep-1)*numbThr,ncol=2,byrow=TRUE)
          red_newSc<-rbind(red_newSc, matRes)
        }
      }else{
        resPar<-foreach(j=1:nSimul, .combine=c, .packages=c('sRACIPE','igraph')) %dopar% 
          inNWsimilarityRefCGr(dataRow, clusterRef, cenMedRef, cutOffM, gene_list, inNo, calc_posM[j,], modelsCGr) 
        
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
      
      if(nUpd>1){
        new_scM[indNew,]<-new_scInM  
      }
      if(nUpd==1)
      {
        new_scM[indNew,]<-as.vector(new_scInM)
      }
      
    }#end if(nUpd>0)
    
    nnew<-setdiff(1:nChains,indNew)
    if(length(nnew)>0){
      for(t in nnew){
        indF<-find_topology(list_sampled, newPosM[t,], 1, nSamp)
        new_scM[t,1]<-scSampV[indF]
        new_scM[t,2]<-thrownOut[indF]
      }
    }
    
    #SEE IF WE ACCEPT NEW SCORES AND THEN DO SWAPS
    
    for( i in 1:nChains){
      #rnUn< min(1, exp((last_sc[1]-new_sc[1])/current_temp)
      if(last_scM[i,1]==new_scM[i,1]){
        accRate[i]<-accRate[i]-1
      }
      rnU<-runif(1)
      if(rnU<exp((last_scM[i,1]-new_scM[i,1])/temp_Array[i])){
        last_scM[i,]<-new_scM[i,]
        posM[i,]<-newPosM[i,]
        accRate[i]<-accRate[i]+1
      }
      
    } 
    
    #WRITE SCORES BEFORE THE SWAPS
    errAcc<-rep(itn,2*nChains+1)
    for(i in 1:nChains){
      errAcc[2*i]<-last_scM[i,1]
      errAcc[2*i+1]<-accRate[i]
    }
    
    write.table(t(errAcc), file=fileOutErrAcc, append=TRUE, sep="\t", row.names=F, col.names=F)
    
    #NOW DO THE SWAPS
    for(t in swEven){
      rnUn<-runif(1)
      swRate<-exp((last_scM[t+1,1]-last_scM[t,1])*(1/temp_Array[t+1]-1/temp_Array[t]))
      if(rnUn<swRate) {    #if accept swap
        swErr<-last_scM[t+1,]
        last_scM[t+1,]<-last_scM[t,]
        last_scM[t,]<-swErr
        swPos<-posM[t+1,]
        posM[t+1,]<-posM[t,]
        posM[t,]<-swPos
        
        
        if((t==1)&&((permArray[2]==1)&&(upDown[1]==-1))){
          numbSwapsA[nChains+2]<-numbSwapsA[nChains+2]+1
        }    
        if(t==1){
          upDown[permArray[2]]<-1
        }
        if((t==nChains-1)&&((permArray[t]==1)&&(upDown[1]==1))){
          numbSwapsA[nChains+1]<-numbSwapsA[nChains+1]+1
        }     
        if(t==nChains-1){
          upDown[permArray[t]]<--1
        }
        if(upDown[permArray[t]]==1){
          nUp[t+1]<-nUp[t+1]+1
        }
        if(upDown[permArray[t+1]]==-1){
          nDown[t]<-nDown[t]+1
        }
        
        swPerm<-permArray[t]
        permArray[t]<-permArray[t+1]
        permArray[t+1]<-swPerm
        
        numbSwapsA[t]<-numbSwapsA[t]+1
      }
      
      if(swRate<1){
        logSwapA[t]<-logSwapA[t]+log(swRate)
      }
    }#end swap chains (1,2),(3,4)...
    
    numbSwapsA[nChains]<-numbSwapsA[nChains]+1 
    
    #### NOW WRITE THE POSITIONS SCORES ARRAYS
    
    matInd<-matrix(rep(itn,2*nChains),nrow=nChains, ncol=2,byrow=T)
    matInd[,2]<-1:nChains
    bi1<-cbind(matInd,posM)
    bi2<-cbind(bi1,last_scM)
    
    write.table(bi2, file=fileOutTop, append=TRUE, sep="\t", row.names=F, col.names=F)
    
    errAcc<-rep(itn,2*nChains+1)
    for(i in 1:nChains){
      errAcc[2*i]<-last_scM[i,1]
      errAcc[2*i+1]<-accRate[i]
    }
    
    write.table(t(errAcc), file=fileOutErrAcc, append=TRUE, sep="\t", row.names=F, col.names=F)
    write.table(t(c(itn,permArray)), file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)
    
    if(itn==numb_iter_extra){
      write.table(t(numbSwapsA), file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)
      write.table(t(permArray), file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)
      write.table(t(nUp), file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)
      write.table(t(nDown), file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)
      write.table(t(upDown), file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)
    }
    
  }#end for(itn in 1:numb_iter_extra)
  
  
  for (t in 1:(nChains-1)){
    logSwapA[t]<-logSwapA[t]/numbSwapsA[nChains]
    if((nUp[t]+nDown[t])>0){
      diffA[t]<-nUp[t]/(nUp[t]+nDown[t])
    }
  }
  
  write.table("Log Swap rates and diffusion function are", file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)
  write.table(t(logSwapA), file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)
  write.table(t(diffA), file=fileOutStuff, append=TRUE, sep="\t", row.names=F, col.names=F)
  
  stopCluster(cl)
  
  allS<-read.table(fileAllSamp,sep="\t")
  allS<-as.matrix(allS)
  colSc<-numbCGRnodes*numbCGRnodes+2
  indMin<-which(allS[,colSc]==min(allS[,colSc]))[1]
  bestTop<-allS[indMin,2:(colSc-1)]
  df<-convAdjTop(bestTop,numbCGRnodes,1:numbCGRnodes)
  
  return(df)
  
}
