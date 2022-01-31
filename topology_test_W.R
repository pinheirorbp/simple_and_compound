# In this script we perform the weighted topology tests in subsets of defined networks
# In the RMarkdown file ("weighted.Rmd" in the folder "analysis/results") statistical analysis is performed for the subset.
# The "aggregation_results_W.R" script is called within the "RMarkdown" file to aggregate results from subsets that were separately tested here
### Preparing files for analysis ####
# Packages
library(bipartite)
library(doSNOW)
library(vegan)
library(stringr)
setwd("~/compound topology")# Define here the home folder
# Functions
source("analysis/functions/PosteriorProb.R")
source("analysis/functions/RestNullModel.R")
## This code is prepared for parallel computation during null model analysis
N_CORES=15 ## Number of parallel cores
N=7 ### Null matrices per core
# File IDs = list of network ids to be analyzed
# Names of file IDS = "IDS_CODE.txt"
# For the weighted networks we included all networks in the list of IDS_set1W.txt (see for format)
## Code for the subset
CODE="set1W"
# File IDs = list of network ids to be analyzed
X= read.table(paste("analysis/IDS_",CODE,".txt",sep=""), header = T, sep="\t")
ids=X$ID
L=length(ids)
#loading binary networks
BIN_NETS=list()
for (i in 1:L){
  if(str_detect(X$Binary_file[i],pattern = ".txt")){ #txt files
    BIN_NETS[[ids[i]]]<-read.table(paste("Networks/",X$Binary_file[i],sep = ""))}
  if(str_detect(X$Binary_file[i],pattern = ".csv")){ #csv files
    BIN_NETS[[ids[i]]]<-read.table(paste("Networks/",X$Binary_file[i],sep = ""), sep = ",")}
  #empty rows and columns
  Rs=rowSums(BIN_NETS[[i]])
  Cs=colSums(BIN_NETS[[i]])
  if(any(Rs==0)){print(paste("binary",ids[[i]], "had", sum(Rs==0),"empty rows"))}
  if(any(Cs==0)){print(paste("binary",ids[[i]], "had", sum(Cs==0),"empty columns"))}
  BIN_NETS[[i]]=BIN_NETS[[i]][Rs>0,Cs>0]
  #binary matrix
  if(any(BIN_NETS[[ids[i]]]!=0&BIN_NETS[[ids[i]]]!=1)){
    print(paste("binary version of",ids[i], "was not really binary"))
    BIN_NETS[[ids[i]]][BIN_NETS[[ids[i]]]>0]=1
  }
}
#loading weighted networks
NETS=list()
for (i in 1:L){
  if(str_detect(X$Weighted_file[i],pattern = ".txt")){ #txt files
    NETS[[ids[i]]]<-read.table(paste("Networks/",X$Weighted_file[i],sep = ""))}
  if(str_detect(X$Weighted_file[i],pattern = ".csv")){ #csv files
    NETS[[ids[i]]]<-read.table(paste("Networks/",X$Weighted_file[i],sep = ""), sep = ",")}
  #empty rows and columns
  Rs=rowSums(NETS[[i]])
  Cs=colSums(NETS[[i]])
  if(any(Rs==0)){print(paste("weighted",ids[[i]], "had", sum(Rs==0),"empty rows"))}
  if(any(Cs==0)){print(paste("weighted",ids[[i]], "had", sum(Cs==0),"empty columns"))}
  NETS[[i]]=NETS[[i]][Rs>0,Cs>0]
  #comparison with binary matrix
  bin2_NET= NETS[[i]]
  bin2_NET[bin2_NET>0]=1
  if((dim(bin2_NET)[1]>dim(bin2_NET)[2])+(dim(BIN_NETS[[i]])[1]>dim(BIN_NETS[[i]])[2])==1){#transpose networks (rows/columns)
    bin2_NET=t(bin2_NET)
    NETS[[i]]=t(NETS[[i]])
  } 
  if(nrow(BIN_NETS[[i]])!=nrow(NETS[[i]])){ # difference in rows
    print(paste(ids[[i]], "  Rows in binary matrix:",nrow(BIN_NETS[[i]])," Rows in weighted matrix: ",nrow(NETS[[i]])))}
  if(ncol(BIN_NETS[[i]])!=ncol(NETS[[i]])){ # difference in cols
    print(paste(ids[[i]], "  Cols in binary matrix:",ncol(BIN_NETS[[i]])," Cols in weighted matrix: ",ncol(NETS[[i]])))}
  conn_BIN_NET=round(sum(BIN_NETS[[i]])/(nrow(BIN_NETS[[i]])*ncol(BIN_NETS[[i]])),digits = 2)
  conn_bin2_NET=round(sum(bin2_NET)/(nrow(bin2_NET)*ncol(bin2_NET)),digits = 2)
  if(conn_BIN_NET!=conn_bin2_NET){ # difference in connectance 
    print(paste(ids[[i]],"  Connectance in binary matrix:",conn_BIN_NET,"Connectance in weighted matrix:",conn_bin2_NET))}
}
## Non integer networks: e.g., prevalence or relative frequencies
net_non_integer=logical()
for (i in 1:L){ 
  net_non_integer[i]=any(NETS[[i]]%%1!=0)
  if(net_non_integer[i]){
    print(paste(ids[[i]], "has non integer values"))
    print(paste(ids[[i]], "rounded with 1 decimal place and multiplied by 10"))
    print("Values < 0.1 converted to 0.1 before this procedure")
    NETS[[i]][NETS[[i]]<0.1&NETS[[i]]>0]=0.1
    NETS[[i]]=round(NETS[[i]],digits = 1)*10
  }
}
rm(Cs,Rs, bin2_NET, conn_BIN_NET,conn_bin2_NET,i)
### size and connectance
X$rows=NA
X$cols=NA
X$connectance=NA
for (i in 1:L){
  X$rows[i]=nrow(NETS[[ids[i]]])
  X$cols[i]=ncol(NETS[[ids[i]]])
  X$connectance[i]=sum(NETS[[ids[i]]]>0)/(X$rows[i]*X$cols[i])
}
### Weighted Nestedness Analysis ####
WNODAcomplete=list()
WNODAsummary=data.frame(WNODA=rep(NA,L),equi_null_mean=rep(NA,L),equi_null_sd=rep(NA,L),equi_z=rep(NA,L),prop_null_mean=rep(NA,L),prop_null_sd=rep(NA,L),prop_z=rep(NA,L))
for (i in c(11,13,15,16)){#c(1:10,12,14,17:L) #c(11,13,15,16)
  print(Sys.time())
  print(paste(ids[i],"Network",i,"in",L,":weighted nestedness analysis start"))
  # observed nestedness
  WNODAcomplete[[i]]=list(ID=ids[i])
  WNODA=nest.smdm(NETS[[ids[i]]],weighted = T,decreasing = "abund")$WNODAmatrix
  WNODAcomplete[[i]]$WNODA=WNODAsummary$WNODA[i]=WNODA
  # null models: equiprobable
  print("equi_null start")
  NR=nrow(NETS[[ids[i]]])
  NC=ncol(NETS[[ids[i]]])
  TOTAL=sum(NETS[[ids[i]]])
  Mprob_equi=matrix(1/(NR*NC),nrow=NR,ncol=NC)
  cl <- makeCluster(N_CORES, outfile="")
  registerDoSNOW(cl)
  wnoda.equi_nulls=foreach(j=1:N_CORES,.packages =c("vegan","bipartite"))%dopar%{
    equi_nulls=list()
    for(nn in 1:N){
      MAT=matrix(0,NR,NC)
      INT=sample(1:(NR*NC),TOTAL,replace=T, Mprob_equi)
      INT_TAB=table(INT)
      MAT[as.numeric(names(INT_TAB))]=INT_TAB
      equi_nulls[[nn]]=MAT
    }
    N_nulls1=list()
    for(nn in 1:N){
      N_nulls1[[nn]]=nest.smdm(equi_nulls[[nn]],weighted = T,decreasing = "abund")$WNODAmatrix
      if((j%%N)==nn){print(paste(i,"    ",ids[i], "   nest_equi_null:",nn," in N: ",N,"   ",Sys.time(),sep=""))}
    }
    N_nulls=list(N_nulls1,equi_nulls)
    return(N_nulls)
  }
  WNODAcomplete[[i]]$equi_nullwnoda= wnodaequi_nulls = as.numeric(sapply(wnoda.equi_nulls, function (y) c(sapply(y[[1]],function (x) x))))
  equi_nulls=sapply (wnoda.equi_nulls, function (y) c(y[[2]])) 
  WNODAcomplete[[i]]$equi_nulls=equi_nulls
  WNODAcomplete[[i]]$equi_null_mean = WNODAsummary$equi_null_mean[i] = mean(wnodaequi_nulls)
  WNODAcomplete[[i]]$equi_null_sd = WNODAsummary$equi_null_sd[i] = sd(wnodaequi_nulls)
  WNODAcomplete[[i]]$equi_zscore = WNODAsummary$equi_z[i] = ((WNODA-WNODAsummary$equi_null_mean[i])/WNODAsummary$equi_null_sd[i])
  rm(equi_nulls,wnoda.equi_nulls)
  # null models: proportional
  print("prop_null start")
  Mprob_prop=PosteriorProb(as.matrix(NETS[[ids[i]]]), R.partitions = 1:NR, C.partitions = 1:NC, Prior.Pij = "degreeprob", Conditional.level ="matrix")
  wnoda.prop_nulls=foreach(j=1:N_CORES,.packages =c("vegan"))%dopar%{
    prop_nulls=list()
    for(nn in 1:N){
      MAT=matrix(0,NR,NC)
      INT=sample(1:(NR*NC),TOTAL,replace=T, Mprob_prop)
      INT_TAB=table(INT)
      MAT[as.numeric(names(INT_TAB))]=INT_TAB
      prop_nulls[[nn]]=MAT
    }
    N_nulls1=list()
    for(nn in 1:N){
      N_nulls1[[nn]]=nest.smdm(prop_nulls[[nn]],weighted = T,decreasing = "abund")$WNODAmatrix
      if((j%%N)==nn){print(paste(i,"    ",ids[i], "   prop_equi_null:",nn," in N: ",N,"   ",Sys.time(),sep=""))}
    }
    N_nulls=list(N_nulls1,prop_nulls)
    return(N_nulls)
  }
  stopCluster(cl)
  WNODAcomplete[[i]]$prop_nullwnoda= wnodaprop_nulls = as.numeric(sapply(wnoda.prop_nulls, function (y) c(sapply(y[[1]],function (x) x))))
  prop_nulls=sapply (wnoda.prop_nulls, function (y) c(y[[2]])) 
  WNODAcomplete[[i]]$prop_nulls=prop_nulls
  WNODAcomplete[[i]]$prop_null_mean = WNODAsummary$prop_null_mean[i] = mean(wnodaprop_nulls)
  WNODAcomplete[[i]]$prop_null_sd = WNODAsummary$prop_null_sd[i] = sd(wnodaprop_nulls)
  WNODAcomplete[[i]]$prop_zscore = WNODAsummary$prop_z[i] = ((WNODA-WNODAsummary$prop_null_mean[i])/WNODAsummary$prop_null_sd[i])
  rm(prop_nulls,wnoda.prop_nulls)
}
WNODAsummary$equi_ci_up= sapply(WNODAcomplete, function(x) quantile(x$equi_nullwnoda,0.975))
WNODAsummary$equi_ci_down= sapply(WNODAcomplete, function(x) quantile(x$equi_nullwnoda,0.025))
WNODAsummary$prop_ci_up= sapply(WNODAcomplete, function(x) quantile(x$prop_nullwnoda,0.975))
WNODAsummary$prop_ci_down= sapply(WNODAcomplete, function(x) quantile(x$prop_nullwnoda,0.025))
WNODAsummary$prop_ci_up2= sapply(WNODAcomplete, function(x) quantile(x$prop_nullwnoda,1))
WNODAsummary$prop_ci_down2= sapply(WNODAcomplete, function(x) quantile(x$prop_nullwnoda,0.05))
names(WNODAcomplete)=rownames(WNODAsummary)=ids
save(WNODAcomplete,WNODAsummary, file=paste("analysis/files/nestedness_",CODE,".RData", sep=""))
### modularity analysis ####
WMODcomplete=list()
WMODsummary=data.frame(modularity=rep(NA,L),nmod=rep(NA,L),equi_null_mean=rep(NA,L),equi_null_sd=rep(NA,L),equi_z=rep(NA,L),equi_null_ci_up=rep(NA,L),equi_null_ci_down=rep(NA,L),prop_null_mean=rep(NA,L),prop_null_sd=rep(NA,L),prop_z=rep(NA,L),prop_null_ci_up=rep(NA,L),prop_null_ci_down=rep(NA,L))
# recovering null models from nestedness analysis - same models
load(paste("analysis/files/nestedness_",CODE,".RData", sep=""))
# equiprobable
equi_nulls_L=list() #equinulls for all networks
for (T in 1:L){
  equi_nulls_L[[T]]=WNODAcomplete[[T]]$equi_nulls
}
# proportional
prop_nulls_L=list() #propnulls for all networks
for (T in 1:L){
  prop_nulls_L[[T]]=WNODAcomplete[[T]]$prop_nulls
}
rm(WNODAcomplete,WNODAsummary,T)
# loop for modularity analysis
for (i in c(11,15)){ #c(1:10,12,14,17:L) #c(13,16) #c(11,15)(LPA)
  print(Sys.time())
  print(paste(ids[i], "Network", i, "in", L))
  # observed modularity
  WMODcomplete[[i]]=list()
  #WMOD=LPA_wb_plus(as.matrix(NETS[[ids[i]]])) # LPA
  WMOD=DIRT_LPA_wb_plus(as.matrix(NETS[[ids[i]]])) # DIRT LPA
  WMODcomplete[[i]]$modularity=WMODsummary$modularity[i]=WMOD$modularity
  # observed number of modules
  row_modules=WMOD$Row_labels
  col_modules=WMOD$Col_labels
  nmod=length(unique(row_modules))
  WMODcomplete[[i]]$nmod=WMODsummary$nmod[i]=nmod
  # observed partitions
  WMODcomplete[[i]]$row_modules=row_modules
  WMODcomplete[[i]]$col_modules=col_modules
  # null models: equiprobable
  print(paste("Modularity: equi_null start",Sys.time()))
  equi_nulls=equi_nulls_L[[i]]
  cl <- makeCluster(N_CORES, outfile="")
  registerDoSNOW(cl)
  mod.equi_nulls=foreach(j=1:N_CORES,.packages =c("bipartite"))%dopar%{
    equi_nulls_N=equi_nulls[((j-1)*N+1):(j*N)] # N equinulls to run in each core
    N_nulls=list()
    for(nn in 1:N){
      #N_nulls[[nn]]=LPA_wb_plus(equi_nulls_N[[nn]]) # LPA
      N_nulls[[nn]]=DIRT_LPA_wb_plus(equi_nulls_N[[nn]]) # DIRT LPA
      if((j%%N)==nn){print(paste(i,"    ",ids[i], "   mod_equi_null:",nn," in N: ",N,"   ",Sys.time(),sep=""))}
    }
    N_nulls
  }  
  WMODcomplete[[i]]$equi_nullmod= like.equi_nulls = as.numeric(sapply(mod.equi_nulls, function (y) c(sapply(y,function (x) x$modularity))))
  WMODcomplete[[i]]$equi_null_mean = WMODsummary$equi_null_mean[i] = mean(like.equi_nulls)
  WMODsummary$equi_null_ci_up[i]=quantile(like.equi_nulls,0.975)
  WMODsummary$equi_null_ci_down[i]=quantile(like.equi_nulls,0.025)
  WMODcomplete[[i]]$equi_null_sd = WMODsummary$equi_null_sd[i] = sd(like.equi_nulls)
  WMODcomplete[[i]]$equi_zscore = WMODsummary$equi_z[i] = ((WMOD$modularity-WMODsummary$equi_null_mean[i])/WMODsummary$equi_null_sd[i])
  rm(equi_nulls)
  # null models: proportional
  print(paste("Modularity: prop_null start",Sys.time()))
  prop_nulls=prop_nulls_L[[i]]
  mod.prop_nulls=foreach(j=1:N_CORES,.packages =c("bipartite"))%dopar%{
    prop_nulls_N=prop_nulls[((j-1)*N+1):(j*N)] # N propnulls to run in each core
    N_nulls=list()
    for(nn in 1:N){
      #N_nulls[[nn]]=LPA_wb_plus(prop_nulls_N[[nn]]) # LPA
      N_nulls[[nn]]=DIRT_LPA_wb_plus(prop_nulls_N[[nn]]) # DIRT LPA
      if((j%%N)==nn){print(paste(i,"    ",ids[i], "   mod_prop_null:",nn," in N: ",N,"   ",Sys.time(),sep=""))}
    }
    N_nulls
  }
  stopCluster(cl)
  WMODcomplete[[i]]$prop_nullmod= like.prop_nulls = as.numeric(sapply(mod.prop_nulls, function (y) c(sapply(y,function (x) x$modularity))))
  WMODcomplete[[i]]$prop_null_mean = WMODsummary$prop_null_mean[i] = mean(like.prop_nulls)
  WMODsummary$prop_null_ci_up[i]=quantile(like.prop_nulls,0.975)
  WMODsummary$prop_null_ci_down[i]=quantile(like.prop_nulls,0.025)
  WMODcomplete[[i]]$prop_null_sd = WMODsummary$prop_null_sd[i] = sd(like.prop_nulls)
  WMODcomplete[[i]]$prop_zscore = WMODsummary$prop_z[i] = ((WMOD$modularity-WMODsummary$prop_null_mean[i])/WMODsummary$prop_null_sd[i])
  rm(prop_nulls)
}
names(WMODcomplete)=rownames(WMODsummary)=ids
save(WMODcomplete,WMODsummary, file=paste("analysis/files/modularity_",CODE,".RData", sep=""))
rm(WMOD,mod.equi_nulls,mod.prop_nulls,equi_nulls_L,prop_nulls_L,i,nmod,like.equi_nulls,like.prop_nulls,row_modules,col_modules)
#### modular networks ####
ISWMOD=(WMODsummary$modularity>WMODsummary$prop_null_ci_up)|(WMODsummary$modularity>WMODsummary$equi_null_ci_up)
WMODTESTED=!is.na(ISWMOD)
ISWMOD[is.na(ISWMOD)]=FALSE
#### nestedness SM DM ####
WCOMPOUNDcomplete=list()
WCOMPOUNDsummary=data.frame(WNODA_SM=rep(NA,L),equi_null_SM_mean=rep(NA,L),equi_null_SM_sd=rep(NA,L),equi_SM_z=rep(NA,L),equi_SM_ci_up=rep(NA,L),equi_SM_ci_down=rep(NA,L),prop_null_SM_mean=rep(NA,L),prop_null_SM_sd=rep(NA,L),prop_SM_z=rep(NA,L),prop_SM_ci_up=rep(NA,L),prop_SM_ci_down=rep(NA,L),prop_SM_ci_up2=rep(NA,L),prop_SM_ci_down2=rep(NA,L),WNODA_DM=rep(NA,L),equi_null_DM_mean=rep(NA,L),equi_null_DM_sd=rep(NA,L),equi_DM_z=rep(NA,L),equi_DM_ci_up=rep(NA,L),equi_DM_ci_down=rep(NA,L),prop_null_DM_mean=rep(NA,L),prop_null_DM_sd=rep(NA,L),prop_DM_z=rep(NA,L),prop_DM_ci_up=rep(NA,L),prop_DM_ci_down=rep(NA,L),prop_DM_ci_up2=rep(NA,L),prop_DM_ci_down2=rep(NA,L))
for (i in c(1:10,12,14,17:L)){#c(1:10,12,14,17:L) 
  if(!WMODTESTED[i]){WCOMPOUNDcomplete[[i]]=NULL}
  if(WMODTESTED[i]){
    print(Sys.time())
    print(paste(ids[i],"Network",sum(ISWMOD[1:i]), " in ", sum(ISWMOD), ": nestedness SM and DM analysis start"))
    WCOMPOUNDcomplete[[i]]=list(ID=ids[i])
    # observed nestedness
    row_modules=WMODcomplete[[i]]$row_modules
    col_modules=WMODcomplete[[i]]$col_modules
    NEST_SM_DM= nest.smdm(NETS[[ids[i]]],constraints = c(row_modules,col_modules),weighted = T,decreasing = "abund", sort = T)
    WCOMPOUNDcomplete[[i]]$WNODA_SM=WCOMPOUNDsummary$WNODA_SM[i]=NEST_SM_DM$WNODA_SM_matrix
    WCOMPOUNDcomplete[[i]]$WNODA_DM=WCOMPOUNDsummary$WNODA_DM[i]=NEST_SM_DM$WNODA_DM_matrix
  }
  # null models (Only for modular networks)
  if(!ISWMOD[i]){print(paste(ids[i]," non modular or not tested yet"))}
  if(ISWMOD[i]){
    # null models: equiprobable
    print("equi_null start")
    Mprob_equi=PosteriorProb(as.matrix(NETS[[ids[i]]]), R.partitions = row_modules, C.partitions = col_modules,Prior.Pij = "equiprobable", Conditional.level ="areas")
    NR=length(row_modules)
    NC=length(col_modules)
    MAT=matrix(0,NR,NC)
    M_areas=matrix(paste(rep(col_modules,each=NR),rep(row_modules,times=NC),sep="X"),nrow = NR,ncol = NC) #areas: modules and adjacencies between modules
    Names_areas=unique(as.character(M_areas))
    cl <- makeCluster(N_CORES, outfile="")
    registerDoSNOW(cl)
    wnoda_SM_DM.equi_nulls=foreach(j=1:N_CORES,.packages =c("vegan","bipartite"))%dopar%{
      equi_nulls=list()
      for(nn in 1:N){ # null matrices
        for(Nam1 in Names_areas){ # conserving frequencies of interactions in each area
          TF=M_areas==Nam1
          TOT_area=sum(NETS[[ids[i]]][TF])
          if(TOT_area>0){
            INT=sample(1:sum(TF),TOT_area,replace=T, Mprob_equi[TF])
            INT_TAB=table(INT)
            MAT[TF][as.numeric(names(INT_TAB))]=INT_TAB
          }
        }
        equi_nulls[[nn]]=MAT
      }
      N_nulls1=list()
      for(nn in 1:N){
        R_TF=rowSums(equi_nulls[[nn]])>0 # removing empty rows...
        C_TF=colSums(equi_nulls[[nn]])>0 # ... and columns
        N_nulls1[[nn]]=nest.smdm(equi_nulls[[nn]][R_TF,C_TF],constraints = c(WMODcomplete[[i]]$row_modules[R_TF],WMODcomplete[[i]]$col_modules[C_TF]),weighted = T,decreasing = "abund", sort = T)
        if((j%%N)==nn){print(paste(i,"    ",ids[i], "   nest_SM_DM_equi_null:",nn," in N: ",N,"   ",Sys.time(),sep=""))}
      }
      N_nulls=list(N_nulls1,equi_nulls)
      return(N_nulls)
    }
    equi_nulls=sapply (wnoda_SM_DM.equi_nulls, function (y) c(y[[2]])) 
    WCOMPOUNDcomplete[[i]]$equi_nulls=equi_nulls
    #SM
    WCOMPOUNDcomplete[[i]]$equi_nullwnoda_SM= wnoda_SM_equi_nulls = as.numeric(sapply(wnoda_SM_DM.equi_nulls, function (y) c(sapply(y[[1]],function (x) x$WNODA_SM_matrix))))
    WCOMPOUNDcomplete[[i]]$equi_null_SM_mean= WCOMPOUNDsummary$equi_null_SM_mean[i] = mean(wnoda_SM_equi_nulls)
    WCOMPOUNDcomplete[[i]]$equi_null_SM_sd = WCOMPOUNDsummary$equi_null_SM_sd[i] = sd(wnoda_SM_equi_nulls)
    WCOMPOUNDcomplete[[i]]$equi_SM_zscore = WCOMPOUNDsummary$equi_SM_z[i] = ((WCOMPOUNDsummary$WNODA_SM[i]-WCOMPOUNDsummary$equi_null_SM_mean[i])/WCOMPOUNDsummary$equi_null_SM_sd[i])
    WCOMPOUNDsummary$equi_SM_ci_up[i]=quantile(wnoda_SM_equi_nulls,0.975, na.rm=T)
    WCOMPOUNDsummary$equi_SM_ci_down[i]=quantile(wnoda_SM_equi_nulls,0.025, na.rm=T)
    #DM
    WCOMPOUNDcomplete[[i]]$equi_nullwnoda_DM= wnoda_DM_equi_nulls = as.numeric(sapply(wnoda_SM_DM.equi_nulls, function (y) c(sapply(y[[1]],function (x) x$WNODA_DM_matrix))))
    WCOMPOUNDcomplete[[i]]$equi_null_DM_mean= WCOMPOUNDsummary$equi_null_DM_mean[i] = mean(wnoda_DM_equi_nulls)
    WCOMPOUNDcomplete[[i]]$equi_null_DM_sd = WCOMPOUNDsummary$equi_null_DM_sd[i] = sd(wnoda_DM_equi_nulls)
    WCOMPOUNDcomplete[[i]]$equi_DM_zscore = WCOMPOUNDsummary$equi_DM_z[i] = ((WCOMPOUNDsummary$WNODA_DM[i]-WCOMPOUNDsummary$equi_null_DM_mean[i])/WCOMPOUNDsummary$equi_null_DM_sd[i])
    WCOMPOUNDsummary$equi_DM_ci_up[i]=quantile(wnoda_DM_equi_nulls,0.975, na.rm=T)
    WCOMPOUNDsummary$equi_DM_ci_down[i]=quantile(wnoda_DM_equi_nulls,0.025, na.rm=T)
    rm(equi_nulls,wnoda_SM_DM.equi_nulls,wnoda_DM_equi_nulls,wnoda_SM_equi_nulls, Mprob_equi)
    
    # null models: proportional
    print("prop_null start")
    Mprob_prop=PosteriorProb(as.matrix(NETS[[ids[i]]]), R.partitions =row_modules, C.partitions = col_modules,Prior.Pij = "degreeprob", Conditional.level ="areas")
    wnoda_SM_DM.prop_nulls=foreach(j=1:N_CORES,.packages =c("vegan","bipartite"))%dopar%{
      prop_nulls=list()
      for(nn in 1:N){ # null matrices
        for(Nam1 in Names_areas){ # conserving frequencies of interactions in each area
          TF=M_areas==Nam1
          TOT_area=sum(NETS[[ids[i]]][TF])
          if(TOT_area>0){
            INT=sample(1:sum(TF),TOT_area,replace=T, Mprob_prop[TF])
            INT_TAB=table(INT)
            MAT[TF][as.numeric(names(INT_TAB))]=INT_TAB
          }
        }
        prop_nulls[[nn]]=MAT
      }
      N_nulls1=list()
      for(nn in 1:N){
        R_TF=rowSums(prop_nulls[[nn]])>0 # removing empty rows...
        C_TF=colSums(prop_nulls[[nn]])>0 # ... and columns
        N_nulls1[[nn]]=nest.smdm(prop_nulls[[nn]][R_TF,C_TF],constraints = c(WMODcomplete[[i]]$row_modules[R_TF],WMODcomplete[[i]]$col_modules[C_TF]),weighted = T,decreasing = "abund", sort = T)
        if((j%%N)==nn){print(paste(i,"    ",ids[i], "   nest_SM_DM_prop_null:",nn," in N: ",N,"   ",Sys.time(),sep=""))}
      }
      N_nulls=list(N_nulls1,prop_nulls)
      return(N_nulls)
    }
    stopCluster(cl)
    rm(cl)
    prop_nulls=sapply (wnoda_SM_DM.prop_nulls, function (y) c(y[[2]])) 
    WCOMPOUNDcomplete[[i]]$prop_nulls=prop_nulls
    #SM
    WCOMPOUNDcomplete[[i]]$prop_nullwnoda_SM= wnoda_SM_prop_nulls = as.numeric(sapply(wnoda_SM_DM.prop_nulls, function (y) c(sapply(y[[1]],function (x) x$WNODA_SM_matrix))))
    WCOMPOUNDcomplete[[i]]$prop_null_SM_mean= WCOMPOUNDsummary$prop_null_SM_mean[i] = mean(wnoda_SM_prop_nulls)
    WCOMPOUNDcomplete[[i]]$prop_null_SM_sd = WCOMPOUNDsummary$prop_null_SM_sd[i] = sd(wnoda_SM_prop_nulls)
    WCOMPOUNDcomplete[[i]]$prop_SM_zscore = WCOMPOUNDsummary$prop_SM_z[i] = ((WCOMPOUNDsummary$WNODA_SM[i]-WCOMPOUNDsummary$prop_null_SM_mean[i])/WCOMPOUNDsummary$prop_null_SM_sd[i])
    WCOMPOUNDsummary$prop_SM_ci_up[i]=quantile(wnoda_SM_prop_nulls,0.975, na.rm=T)
    WCOMPOUNDsummary$prop_SM_ci_down[i]=quantile(wnoda_SM_prop_nulls,0.025, na.rm=T)
    WCOMPOUNDsummary$prop_SM_ci_up2[i]=quantile(wnoda_SM_prop_nulls,1, na.rm=T)
    WCOMPOUNDsummary$prop_SM_ci_down2[i]=quantile(wnoda_SM_prop_nulls,0.05, na.rm=T)
    #DM
    WCOMPOUNDcomplete[[i]]$prop_nullwnoda_DM= wnoda_DM_prop_nulls = as.numeric(sapply(wnoda_SM_DM.prop_nulls, function (y) c(sapply(y[[1]],function (x) x$WNODA_DM_matrix))))
    WCOMPOUNDcomplete[[i]]$prop_null_DM_mean= WCOMPOUNDsummary$prop_null_DM_mean[i] = mean(wnoda_DM_prop_nulls)
    WCOMPOUNDcomplete[[i]]$prop_null_DM_sd = WCOMPOUNDsummary$prop_null_DM_sd[i] = sd(wnoda_DM_prop_nulls)
    WCOMPOUNDcomplete[[i]]$prop_DM_zscore = WCOMPOUNDsummary$prop_DM_z[i] = ((WCOMPOUNDsummary$WNODA_DM[i]-WCOMPOUNDsummary$prop_null_DM_mean[i])/WCOMPOUNDsummary$prop_null_DM_sd[i])
    WCOMPOUNDsummary$prop_DM_ci_up[i]=quantile(wnoda_DM_prop_nulls,0.975,na.rm=T)
    WCOMPOUNDsummary$prop_DM_ci_down[i]=quantile(wnoda_DM_prop_nulls,0.025,na.rm=T)
    WCOMPOUNDsummary$prop_DM_ci_up2[i]=quantile(wnoda_DM_prop_nulls,1,na.rm=T)
    WCOMPOUNDsummary$prop_DM_ci_down2[i]=quantile(wnoda_DM_prop_nulls,0.05,na.rm=T)
    rm(prop_nulls,wnoda_SM_DM.prop_nulls,wnoda_DM_prop_nulls,wnoda_SM_prop_nulls, Mprob_prop)
  }
  rm(NEST_SM_DM)
}
names(WCOMPOUNDcomplete)=rownames(WCOMPOUNDsummary)=ids
save(WCOMPOUNDcomplete,WCOMPOUNDsummary, file=paste("analysis/files/compound_",CODE,".RData", sep=""))
## This script creates 3 files with results of tests in the folder "analysis/files/"
## 1 modularity_CODE.RData
## 2 Nestedness_CODE.RData 
## 3 Compound tests_CODE.RData