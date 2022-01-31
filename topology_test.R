# In this script we perform the binary topology tests in subsets of defined networks
# In the RMarkdown file ("binary.Rmd" in the folder "analysis/results") statistical analysis is performed for the subset.
# The "aggregation_results.R" script is called within the "RMarkdown" file to aggregate results from subsets that were separately tested here
# Packages
library(bipartite)
library(doSNOW)
library(vegan)
library(stringr)
setwd("~/compound topology") # Define here the home folder
# Functions
source("analysis/functions/PosteriorProb.R")
source("analysis/functions/RestNullModel.R")
## This code is prepared for parallel computation during null model analysis
N_CORES=35 ## Number of parallel cores
N=3 ### Null matrices per core

# File IDs = list of network ids to be analyzed
# Names of file IDS = "IDS_CODE.txt"
# We did not run all networks together, but in three subsets, indicated in the files IDS_set1.txt, IDS_set2.txt and IDS_set3.txt (see one of these for the format)
## Code for the subset
CODE="set3"
##
X= read.table(paste("analysis/IDS_",CODE,".txt",sep=""), header = T, sep="\t")
ids=X$ID
L=length(ids)
#loading networks
NETS=list()
for (i in 1:L){
  if(str_detect(X$File[i],pattern = ".txt")){ #txt files
    NETS[[ids[i]]]<-read.table(paste("Networks/",X$File[i],sep = ""))}
  if(str_detect(X$File[i],pattern = ".csv")){ #csv files
    NETS[[ids[i]]]<-read.table(paste("Networks/",X$File[i],sep = ""), sep = ",")}
  #empty rows and columns
  Rs=rowSums(NETS[[i]])
  Cs=colSums(NETS[[i]])
  if(any(Rs==0)){print(paste(ids[[i]], "had", sum(Rs==0),"empty rows"))}
  if(any(Cs==0)){print(paste(ids[[i]], "had", sum(Cs==0),"empty columns"))}
  NETS[[i]]=NETS[[i]][Rs>0,Cs>0]
  #binary matrix
  if(any(NETS[[ids[i]]]!=0&NETS[[ids[i]]]!=1)){
    print(paste(ids[i], "was not binary"))
    NETS[[ids[i]]][NETS[[ids[i]]]>0]=1
  }
}
rm(Cs,Rs)
### size and connectance
X$rows=NA
X$cols=NA
X$connectance=NA
for (i in 1:L){
  X$rows[i]=nrow(NETS[[ids[i]]])
  X$cols[i]=ncol(NETS[[ids[i]]])
  X$connectance[i]=sum(NETS[[ids[i]]])/(X$rows[i]*X$cols[i])
}

### nestedness analysis ####
NODFcomplete=list()
NODFsummary=data.frame(NODF=rep(NA,L),equi_null_mean=rep(NA,L),equi_null_sd=rep(NA,L),equi_z=rep(NA,L),prop_null_mean=rep(NA,L),prop_null_sd=rep(NA,L),prop_z=rep(NA,L))
for (i in 1:L){
  print(Sys.time())
  print(paste(ids[i],"Network",i,"in",L,": nestedness analysis start"))
  # observed nestedness
  NODFcomplete[[i]]=list(ID=ids[i])
  NODF=nestednodf(as.matrix(NETS[[ids[i]]]))$statistic[[3]]
  NODFcomplete[[i]]$NODF=NODFsummary$NODF[i]=NODF
  # null models: equiprobable
  print("equi_null start")
  Mprob_equi=matrix(1/(nrow(NETS[[ids[i]]])*ncol(NETS[[ids[i]]])),nrow=nrow(NETS[[ids[i]]]),ncol=ncol(NETS[[ids[i]]]))
  cl <- makeCluster(N_CORES, outfile="")
  registerDoSNOW(cl)
  nodf.equi_nulls=foreach(j=1:N_CORES,.packages =c("vegan"))%dopar%{
    equi_nulls=RestNullModel(as.matrix(NETS[[ids[i]]]),Pij.Prob = Mprob_equi,Numbernulls = N)
    N_nulls1=list()
    for(nn in 1:N){
      N_nulls1[[nn]]=nestednodf(equi_nulls[[nn]]$NullMatrix)
      if(j==nn){print(paste(i,"    ",ids[i], "   nest_equi_null:",nn," in N: ",N,"   ",Sys.time(),sep=""))}
    }
    N_nulls=list(N_nulls1,equi_nulls)
    return(N_nulls)
  }
  NODFcomplete[[i]]$equi_nullnodf= nodfequi_nulls = as.numeric(sapply(nodf.equi_nulls, function (y) c(sapply(y[[1]],function (x) x$statistic[[3]]))))
  equi_nulls=sapply (nodf.equi_nulls, function (y) c(y[[2]])) 
  NODFcomplete[[i]]$equi_nulls=equi_nulls
  NODFcomplete[[i]]$equi_null_mean = NODFsummary$equi_null_mean[i] = mean(nodfequi_nulls)
  NODFcomplete[[i]]$equi_null_sd = NODFsummary$equi_null_sd[i] = sd(nodfequi_nulls)
  NODFcomplete[[i]]$equi_zscore = NODFsummary$equi_z[i] = ((NODF-NODFsummary$equi_null_mean[i])/NODFsummary$equi_null_sd[i])
  rm(equi_nulls,nodf.equi_nulls)
  # null models: proportional
  print("prop_null start")
  Mprob_prop=PosteriorProb(as.matrix(NETS[[ids[i]]]), R.partitions = 1:nrow(NETS[[ids[i]]]), C.partitions = 1:ncol(NETS[[ids[i]]]), Prior.Pij = "degreeprob", Conditional.level ="matrix")
  nodf.prop_nulls=foreach(j=1:N_CORES,.packages =c("vegan"))%dopar%{
    prop_nulls=RestNullModel(as.matrix(NETS[[ids[i]]]),Pij.Prob = Mprob_prop,Numbernulls = N)
    N_nulls1=list()
    for(nn in 1:N){
      N_nulls1[[nn]]=nestednodf(prop_nulls[[nn]]$NullMatrix)
      if(j==nn){print(paste(i,"    ",ids[i], "   nest_prop_null:",nn," in N: ",N,"   ",Sys.time(),sep=""))}
    }
    N_nulls=list(N_nulls1,prop_nulls)
    return(N_nulls)
  }
  stopCluster(cl)
  NODFcomplete[[i]]$prop_nullnodf= nodfprop_nulls = as.numeric(sapply(nodf.prop_nulls, function (y) c(sapply(y[[1]],function (x) x$statistic[[3]]))))
  prop_nulls=sapply (nodf.prop_nulls, function (y) c(y[[2]])) 
  NODFcomplete[[i]]$prop_nulls=prop_nulls
  NODFcomplete[[i]]$prop_null_mean = NODFsummary$prop_null_mean[i] = mean(nodfprop_nulls)
  NODFcomplete[[i]]$prop_null_sd = NODFsummary$prop_null_sd[i] = sd(nodfprop_nulls)
  NODFcomplete[[i]]$prop_zscore = NODFsummary$prop_z[i] = ((NODF-NODFsummary$prop_null_mean[i])/NODFsummary$prop_null_sd[i])
  rm(prop_nulls,nodf.prop_nulls)
}
NODFsummary$equi_ci_up= sapply(NODFcomplete, function(x) quantile(x$equi_nullnodf,0.975))
NODFsummary$equi_ci_down= sapply(NODFcomplete, function(x) quantile(x$equi_nullnodf,0.025))
NODFsummary$prop_ci_up= sapply(NODFcomplete, function(x) quantile(x$prop_nullnodf,0.975))
NODFsummary$prop_ci_down= sapply(NODFcomplete, function(x) quantile(x$prop_nullnodf,0.025))
names(NODFcomplete)=rownames(NODFsummary)=ids
save(NODFcomplete,NODFsummary, file=paste("analysis/files/nestedness_",CODE,".RData", sep=""))

### modularity analysis ####
MODcomplete=list()
MODsummary=data.frame(modularity=rep(NA,L),nmod=rep(NA,L),equi_null_mean=rep(NA,L),equi_null_sd=rep(NA,L),equi_z=rep(NA,L),equi_null_ci_up=rep(NA,L),equi_null_ci_down=rep(NA,L),prop_null_mean=rep(NA,L),prop_null_sd=rep(NA,L),prop_z=rep(NA,L),prop_null_ci_up=rep(NA,L),prop_null_ci_down=rep(NA,L))
# recovering null models from nestedness analysis - same models
load(paste("analysis/files/nestedness_",CODE,".RData", sep=""))
# equiprobable
equi_nulls_L=list() #equinulls for all networks
for (T in 1:L){
  equi_nulls_L[[T]]=NODFcomplete[[T]]$equi_nulls
}
# proportional
prop_nulls_L=list() #propnulls for all networks
for (T in 1:L){
  prop_nulls_L[[T]]=NODFcomplete[[T]]$prop_nulls
}
rm(NODFcomplete,NODFsummary,T)
# loop for modularity analysis
for (i in 1:L){ 
  print(Sys.time())
  print(paste(ids[i], "Network", i, "in", L))
  # observed modularity
  MODcomplete[[i]]=list()
  if((X$rows[i]+X$cols[i])>=600){MOD=LPA_wb_plus(as.matrix(NETS[[ids[i]]]))}
  if((X$rows[i]+X$cols[i])<600){MOD=DIRT_LPA_wb_plus(as.matrix(NETS[[ids[i]]]))}
  MODcomplete[[i]]$modularity=MODsummary$modularity[i]=MOD$modularity
  # observed number of modules
  row_modules=MOD$Row_labels
  col_modules=MOD$Col_labels
  nmod=length(unique(row_modules))
  MODcomplete[[i]]$nmod=MODsummary$nmod[i]=nmod
  # observed partitions
  MODcomplete[[i]]$row_modules=row_modules
  MODcomplete[[i]]$col_modules=col_modules
  # null models: equiprobable
  print(paste("Modularity: equi_null start",Sys.time()))
  equi_nulls=equi_nulls_L[[i]]
  cl <- makeCluster(N_CORES, outfile="")
  registerDoSNOW(cl)
  mod.equi_nulls=foreach(j=1:N_CORES,.packages =c("bipartite"))%dopar%{
    equi_nulls_N=equi_nulls[((j-1)*N+1):(j*N)] # N equinulls to run in each core
    N_nulls=list()
    for(nn in 1:N){
      if((X$rows[i]+X$cols[i])>=600){
        N_nulls[[nn]]=LPA_wb_plus(equi_nulls_N[[nn]]$NullMatrix)}
      if((X$rows[i]+X$cols[i])<600){
        N_nulls[[nn]]=DIRT_LPA_wb_plus(equi_nulls_N[[nn]]$NullMatrix)}
      if((j%%N)==(nn%%N)){print(paste(i,"    ",ids[i], "   mod_equi_null:",nn," in N: ",N,"   ",Sys.time(),sep=""))}
    }
    N_nulls
  }  
  MODcomplete[[i]]$equi_nullmod= like.equi_nulls = as.numeric(sapply(mod.equi_nulls, function (y) c(sapply(y,function (x) x$modularity))))
  MODcomplete[[i]]$equi_null_mean = MODsummary$equi_null_mean[i] = mean(like.equi_nulls)
  MODsummary$equi_null_ci_up[i]=quantile(like.equi_nulls,0.975)
  MODsummary$equi_null_ci_down[i]=quantile(like.equi_nulls,0.025)
  MODcomplete[[i]]$equi_null_sd = MODsummary$equi_null_sd[i] = sd(like.equi_nulls)
  MODcomplete[[i]]$equi_zscore = MODsummary$equi_z[i] = ((MOD$modularity-MODsummary$equi_null_mean[i])/MODsummary$equi_null_sd[i])
  rm(equi_nulls)
  # null models: proportional
  print(paste("Modularity: prop_null start",Sys.time()))
  prop_nulls=prop_nulls_L[[i]]
  mod.prop_nulls=foreach(j=1:N_CORES,.packages =c("bipartite"))%dopar%{
    prop_nulls_N=prop_nulls[((j-1)*N+1):(j*N)] # N propnulls to run in each core
    N_nulls=list()
    for(nn in 1:N){
      if((X$rows[i]+X$cols[i])>=600){
        N_nulls[[nn]]=LPA_wb_plus(prop_nulls_N[[nn]]$NullMatrix)}
      if((X$rows[i]+X$cols[i])<600){
        N_nulls[[nn]]=DIRT_LPA_wb_plus(prop_nulls_N[[nn]]$NullMatrix)}
      if((j%%N)==(nn%%N)){print(paste(i,"    ",ids[i], "   mod_prop_null:",nn," in N: ",N,"   ",Sys.time(),sep=""))}
    }
    N_nulls
  }
  stopCluster(cl)
  MODcomplete[[i]]$prop_nullmod= like.prop_nulls = as.numeric(sapply(mod.prop_nulls, function (y) c(sapply(y,function (x) x$modularity))))
  MODcomplete[[i]]$prop_null_mean = MODsummary$prop_null_mean[i] = mean(like.prop_nulls)
  MODsummary$prop_null_ci_up[i]=quantile(like.prop_nulls,0.975)
  MODsummary$prop_null_ci_down[i]=quantile(like.prop_nulls,0.025)
  MODcomplete[[i]]$prop_null_sd = MODsummary$prop_null_sd[i] = sd(like.prop_nulls)
  MODcomplete[[i]]$prop_zscore = MODsummary$prop_z[i] = ((MOD$modularity-MODsummary$prop_null_mean[i])/MODsummary$prop_null_sd[i])
  rm(prop_nulls)
}
names(MODcomplete)=rownames(MODsummary)=ids
save(MODcomplete,MODsummary, file=paste("analysis/files/modularity_",CODE,".RData", sep=""))
rm(MOD,mod.equi_nulls,mod.prop_nulls,equi_nulls_L,prop_nulls_L,i,nmod,like.equi_nulls,like.prop_nulls,row_modules,col_modules)
#### modular networks ####
ISMOD=(MODsummary$modularity>MODsummary$prop_null_ci_up)|(MODsummary$modularity>MODsummary$equi_null_ci_up)
### Nestedness SM and DM ####
COMPOUNDcomplete=list()
COMPOUNDsummary=data.frame(NODF_SM=rep(NA,L),equi_null_SM_mean=rep(NA,L),equi_null_SM_sd=rep(NA,L),equi_SM_z=rep(NA,L),equi_SM_ci_up=rep(NA,L),equi_SM_ci_down=rep(NA,L),prop_null_SM_mean=rep(NA,L),prop_null_SM_sd=rep(NA,L),prop_SM_z=rep(NA,L),prop_SM_ci_up=rep(NA,L),prop_SM_ci_down=rep(NA,L),NODF_DM=rep(NA,L),equi_null_DM_mean=rep(NA,L),equi_null_DM_sd=rep(NA,L),equi_DM_z=rep(NA,L),equi_DM_ci_up=rep(NA,L),equi_DM_ci_down=rep(NA,L),prop_null_DM_mean=rep(NA,L),prop_null_DM_sd=rep(NA,L),prop_DM_z=rep(NA,L),prop_DM_ci_up=rep(NA,L),prop_DM_ci_down=rep(NA,L))
for (i in 1:L){ 
  print(Sys.time())
  print(paste(ids[i],"Network",sum(ISMOD[1:i]), " in ", sum(ISMOD), ": nestedness SM and DM analysis start"))
  COMPOUNDcomplete[[i]]=list(ID=ids[i])
  # observed nestedness
  NEST_SM_DM= nest.smdm(NETS[[ids[i]]],constraints = c(MODcomplete[[i]]$row_modules,MODcomplete[[i]]$col_modules),weighted = F,decreasing = "fill", sort = T)
  COMPOUNDcomplete[[i]]$NODF_SM=COMPOUNDsummary$NODF_SM[i]=NEST_SM_DM$NODF_SM_matrix
  COMPOUNDcomplete[[i]]$NODF_DM=COMPOUNDsummary$NODF_DM[i]=NEST_SM_DM$NODF_DM_matrix
  # null models (Only for modular networks)
  if(!ISMOD[i]){print(paste(ids[i]," non modular network"))}
  if(ISMOD[i]){
    # null models: equiprobable
    print("equi_null start")
    Mprob_equi=PosteriorProb(as.matrix(NETS[[ids[i]]]), R.partitions = MODcomplete[[i]]$row_modules, C.partitions = MODcomplete[[i]]$col_modules,Prior.Pij = "equiprobable", Conditional.level ="areas")
    cl <- makeCluster(N_CORES, outfile="")
    registerDoSNOW(cl)
    nodf_SM_DM.equi_nulls=foreach(j=1:N_CORES,.packages =c("vegan","bipartite"))%dopar%{
      equi_nulls=RestNullModel(as.matrix(NETS[[ids[i]]]),Pij.Prob = Mprob_equi,Numbernulls = N,byarea = T,R.partitions = MODcomplete[[i]]$row_modules, C.partitions = MODcomplete[[i]]$col_modules)
      N_nulls1=list()
      for(nn in 1:N){
        N_nulls1[[nn]]=nest.smdm(equi_nulls[[nn]]$NullMatrix,constraints = c(MODcomplete[[i]]$row_modules,MODcomplete[[i]]$col_modules),weighted = F,decreasing = "fill", sort = T)
        if((j%%N)==nn){print(paste(i,"    ",ids[i], "   nest_SM_DM_equi_null:",nn," in N: ",N,"   ",Sys.time(),sep=""))}
      }
      N_nulls=list(N_nulls1,equi_nulls)
      return(N_nulls)
    }
    equi_nulls=sapply (nodf_SM_DM.equi_nulls, function (y) c(y[[2]])) 
    COMPOUNDcomplete[[i]]$equi_nulls=equi_nulls
    #SM
    COMPOUNDcomplete[[i]]$equi_nullnodf_SM= nodf_SM_equi_nulls = as.numeric(sapply(nodf_SM_DM.equi_nulls, function (y) c(sapply(y[[1]],function (x) x$NODF_SM_matrix))))
    COMPOUNDcomplete[[i]]$equi_null_SM_mean= COMPOUNDsummary$equi_null_SM_mean[i] = mean(nodf_SM_equi_nulls)
    COMPOUNDcomplete[[i]]$equi_null_SM_sd = COMPOUNDsummary$equi_null_SM_sd[i] = sd(nodf_SM_equi_nulls)
    COMPOUNDcomplete[[i]]$equi_SM_zscore = COMPOUNDsummary$equi_SM_z[i] = ((COMPOUNDsummary$NODF_SM[i]-COMPOUNDsummary$equi_null_SM_mean[i])/COMPOUNDsummary$equi_null_SM_sd[i])
    COMPOUNDsummary$equi_SM_ci_up[i]=quantile(nodf_SM_equi_nulls,0.975, na.rm=T)
    COMPOUNDsummary$equi_SM_ci_down[i]=quantile(nodf_SM_equi_nulls,0.025, na.rm=T)
    #DM
    COMPOUNDcomplete[[i]]$equi_nullnodf_DM= nodf_DM_equi_nulls = as.numeric(sapply(nodf_SM_DM.equi_nulls, function (y) c(sapply(y[[1]],function (x) x$NODF_DM_matrix))))
    COMPOUNDcomplete[[i]]$equi_null_DM_mean= COMPOUNDsummary$equi_null_DM_mean[i] = mean(nodf_DM_equi_nulls)
    COMPOUNDcomplete[[i]]$equi_null_DM_sd = COMPOUNDsummary$equi_null_DM_sd[i] = sd(nodf_DM_equi_nulls)
    COMPOUNDcomplete[[i]]$equi_DM_zscore = COMPOUNDsummary$equi_DM_z[i] = ((COMPOUNDsummary$NODF_DM[i]-COMPOUNDsummary$equi_null_DM_mean[i])/COMPOUNDsummary$equi_null_DM_sd[i])
    COMPOUNDsummary$equi_DM_ci_up[i]=quantile(nodf_DM_equi_nulls,0.975, na.rm=T)
    COMPOUNDsummary$equi_DM_ci_down[i]=quantile(nodf_DM_equi_nulls,0.025, na.rm=T)
    rm(equi_nulls,nodf_SM_DM.equi_nulls,nodf_DM_equi_nulls,nodf_SM_equi_nulls, Mprob_equi)
    
    # null models: proportional
    print("prop_null start")
    Mprob_prop=PosteriorProb(as.matrix(NETS[[ids[i]]]), R.partitions = MODcomplete[[i]]$row_modules, C.partitions = MODcomplete[[i]]$col_modules,Prior.Pij = "degreeprob", Conditional.level ="areas")
    nodf_SM_DM.prop_nulls=foreach(j=1:N_CORES,.packages =c("vegan","bipartite"))%dopar%{
      prop_nulls=RestNullModel(as.matrix(NETS[[ids[i]]]),Pij.Prob = Mprob_prop,Numbernulls = N,byarea = T,R.partitions = MODcomplete[[i]]$row_modules, C.partitions = MODcomplete[[i]]$col_modules)
      N_nulls1=list()
      for(nn in 1:N){
        N_nulls1[[nn]]=nest.smdm(prop_nulls[[nn]]$NullMatrix,constraints = c(MODcomplete[[i]]$row_modules,MODcomplete[[i]]$col_modules),weighted = F,decreasing = "fill", sort = T)
        if((j%%N)==nn){print(paste(i,"    ",ids[i], "   nest_SM_DM_prop_null:",nn," in N: ",N,"   ",Sys.time(),sep=""))}
      }
      N_nulls=list(N_nulls1,prop_nulls)
      return(N_nulls)
    }
    stopCluster(cl)
    rm(cl)
    prop_nulls=sapply (nodf_SM_DM.prop_nulls, function (y) c(y[[2]])) 
    COMPOUNDcomplete[[i]]$prop_nulls=prop_nulls
    #SM
    COMPOUNDcomplete[[i]]$prop_nullnodf_SM= nodf_SM_prop_nulls = as.numeric(sapply(nodf_SM_DM.prop_nulls, function (y) c(sapply(y[[1]],function (x) x$NODF_SM_matrix))))
    COMPOUNDcomplete[[i]]$prop_null_SM_mean= COMPOUNDsummary$prop_null_SM_mean[i] = mean(nodf_SM_prop_nulls)
    COMPOUNDcomplete[[i]]$prop_null_SM_sd = COMPOUNDsummary$prop_null_SM_sd[i] = sd(nodf_SM_prop_nulls)
    COMPOUNDcomplete[[i]]$prop_SM_zscore = COMPOUNDsummary$prop_SM_z[i] = ((COMPOUNDsummary$NODF_SM[i]-COMPOUNDsummary$prop_null_SM_mean[i])/COMPOUNDsummary$prop_null_SM_sd[i])
    COMPOUNDsummary$prop_SM_ci_up[i]=quantile(nodf_SM_prop_nulls,0.975, na.rm=T)
    COMPOUNDsummary$prop_SM_ci_down[i]=quantile(nodf_SM_prop_nulls,0.025, na.rm=T)
    #DM
    COMPOUNDcomplete[[i]]$prop_nullnodf_DM= nodf_DM_prop_nulls = as.numeric(sapply(nodf_SM_DM.prop_nulls, function (y) c(sapply(y[[1]],function (x) x$NODF_DM_matrix))))
    COMPOUNDcomplete[[i]]$prop_null_DM_mean= COMPOUNDsummary$prop_null_DM_mean[i] = mean(nodf_DM_prop_nulls)
    COMPOUNDcomplete[[i]]$prop_null_DM_sd = COMPOUNDsummary$prop_null_DM_sd[i] = sd(nodf_DM_prop_nulls)
    COMPOUNDcomplete[[i]]$prop_DM_zscore = COMPOUNDsummary$prop_DM_z[i] = ((COMPOUNDsummary$NODF_DM[i]-COMPOUNDsummary$prop_null_DM_mean[i])/COMPOUNDsummary$prop_null_DM_sd[i])
    COMPOUNDsummary$prop_DM_ci_up[i]=quantile(nodf_DM_prop_nulls,0.975,na.rm=T)
    COMPOUNDsummary$prop_DM_ci_down[i]=quantile(nodf_DM_prop_nulls,0.025,na.rm=T)
    rm(prop_nulls,nodf_SM_DM.prop_nulls,nodf_DM_prop_nulls,nodf_SM_prop_nulls, Mprob_prop)
  }
  rm(NEST_SM_DM)
}
names(COMPOUNDcomplete)=rownames(COMPOUNDsummary)=ids
save(COMPOUNDcomplete,COMPOUNDsummary, file=paste("analysis/files/compound_",CODE,".RData", sep=""))
## This script creates 3 files with results of tests in the folder "analysis/files/"
## 1 modularity_CODE.RData
## 2 Nestedness_CODE.RData 
## 3 Compound tests_CODE.RData