### Necessary loadings, packages and inputs ####
library(FSA)
library(stringr)
library(bipartite)
options(digits=3)
##
CODES=c("set1","set2","set3")
NAME_OUTPUT=paste(CODES,collapse = "_")
COLORS1=c(1,"#1E46D9","#F80631","#00D400","#FFA100","#7A05C5")
COLORS2=c("#AF2E45","#ffd143","#353082","#4EA02A")
TABLE_RESULTS=read.table(paste("results/",NAME_OUTPUT,".txt",sep=""), sep="\t", header=T, row.names = 1)
X=read.table(paste("results/","NETDATA_",NAME_OUTPUT,".txt",sep=""), sep="\t", header=T, row.names = 1)
rownames(TABLE_RESULTS)=X$ID
load(paste("results/MODULES_",NAME_OUTPUT,".RData",sep=""))
## loading networks
ids=X$ID
L=length(ids)
NETS=list()
for (i in 1:L){
  if(str_detect(X$File[i],pattern = ".txt")){ #txt files
    NETS[[ids[i]]]<-read.table(paste("networks/",X$File[i],sep = ""))}
  if(str_detect(X$File[i],pattern = ".csv")){ #csv files
    NETS[[ids[i]]]<-read.table(paste("networks/",X$File[i],sep = ""), sep = ",")}
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
##
ROB_RESULTS=TABLE_RESULTS[,c(1,4:8)]
#### Removal of rows####

# All networks by decreasing degrees
# ROB_RESULTS$rob_low_deg
ROB_RESULTS$rob_low_deg=NA
for (Ri in 1:nrow(ROB_RESULTS)){
  net=NETS[[Ri]]
  OR_DEG=order(rowSums(net), decreasing = T)
  SE=second.extinct(net, participant = "lower", method="external", ext.row = OR_DEG)
  ROB_RESULTS$rob_low_deg[Ri]=robustness(SE)
}
# All networks: random sequences (1000 simulations)
# RANDOMROB_RESULTS_L= robustness for all simulations
# ROB_RESULTS$rob_low_ran= means
RANDOMROB_RESULTS_L=list()
ROB_RESULTS$rob_low_ran=NA
for ( Ri in 1:nrow(ROB_RESULTS)){
    print(names(NETS)[Ri])
    net=NETS[[Ri]]
    rr_results=numeric()
    for (IT in 1:1000){
      if(IT%%100==0){print(IT)}
      OR_ran=sample(1:nrow(net))
      SE=second.extinct(net, participant = "lower", method="external", ext.row = OR_ran)
      rr_results[IT]=robustness(SE)
    }
    RANDOMROB_RESULTS_L[[Ri]]=rr_results
    ROB_RESULTS$rob_low_ran[Ri]=mean(rr_results)
  } 

# Only compound networks (1000 simulations)
# Random order of modules
# Randomly removes species from the first module until the module is fully excluded, then moves to the next and so on
# COMPOUND1ROB_RESULTS_L= robustness for all simulations
# ROB_RESULTS$rob_low_mod= means
COMPOUND1ROB_RESULTS_L=list()
ROB_RESULTS$rob_low_mod=NA
for ( Ri in 1:nrow(ROB_RESULTS)){
  if(!is.na(ROB_RESULTS$topology[Ri])){
    if(ROB_RESULTS$topology[Ri]=="compound"){
      print(names(NETS)[Ri])
      net=NETS[[Ri]]
      rr_results=numeric()
      mod1=MODULES[[Ri]]$rows
      for (IT in 1:1000){
        if(IT%%100==0){print(IT)}
        OR_ran=numeric()
        n2=1:nrow(net)
        n3=rep(TRUE,nrow(net))
        mod2=sample(unique(mod1))
        MO=1
        for (n1 in n2){
          n3=(mod1==mod2[MO]&!is.element(n2,OR_ran))
          n4=ifelse(sum(n3)==1,n2[n3],sample(n2[n3],1))
          OR_ran[n1]=n4
          if(sum(n3)==1){MO=MO+1}
        }
        SE=second.extinct(net, participant = "lower", method="external", ext.row = OR_ran)
        rr_results[IT]=robustness(SE)
      }
      COMPOUND1ROB_RESULTS_L[[Ri]]=rr_results
      ROB_RESULTS$rob_low_mod[Ri]=mean(rr_results)
    }
  }
}
# Only compound networks (1000 simulations)
# Avoid excluding modules
# In each round, randonly removes a species from the largest module(s)
# COMPOUND2ROB_RESULTS_L= robustness for all simulations
# ROB_RESULTS$rob_low_mod2= means
COMPOUND2ROB_RESULTS_L=list()
ROB_RESULTS$rob_low_mod2=NA
for ( Ri in 1:nrow(ROB_RESULTS)){
  if(!is.na(ROB_RESULTS$topology[Ri])){
    if(ROB_RESULTS$topology[Ri]=="compound"){
      print(names(NETS)[Ri])
      net=NETS[[Ri]]
      rr_results=numeric()
      mod1=MODULES[[Ri]]$rows
      for (IT in 1:1000){
        if(IT%%100==0){print(IT)}
        OR_ran=numeric()
        n2=1:nrow(net)
        n3=rep(TRUE,nrow(net))
        for (n1 in n2){
          n3=(!is.element(n2,OR_ran))
          MAX1=max(table(mod1[n3]))
          mod2=names(table(mod1[n3]))[table(mod1[n3])==MAX1]
          n3=(is.element(mod1,mod2)&!is.element(n2,OR_ran))
          n4=ifelse(sum(n3)==1,n2[n3],sample(n2[n3],1))
          OR_ran[n1]=n4
        }
        SE=second.extinct(net, participant = "lower", method="external", ext.row = OR_ran)
        rr_results[IT]=robustness(SE)
      }
      COMPOUND2ROB_RESULTS_L[[Ri]]=rr_results
      ROB_RESULTS$rob_low_mod2[Ri]=mean(rr_results)
    }
  }
}
save(ROB_RESULTS,RANDOMROB_RESULTS_L,COMPOUND1ROB_RESULTS_L,COMPOUND2ROB_RESULTS_L, file=paste("results/robustness_",NAME_OUTPUT,".RData",sep=""))

#### Removal of columns ####

# All networks by decreasing degrees
# ROB_RESULTS$rob_high_deg
ROB_RESULTS$rob_high_deg=NA
for (Ri in 1:nrow(ROB_RESULTS)){
  net=NETS[[Ri]]
  OR_DEG=order(colSums(net), decreasing = T)
  SE=second.extinct(net, participant = "higher", method="external", ext.col = OR_DEG)
  ROB_RESULTS$rob_high_deg[Ri]=robustness(SE)
}
# All networks: random sequences (1000 simulations)
# RANDOMROB_RESULTS_H= robustness for all simulations
# ROB_RESULTS$rob_high_ran= means
RANDOMROB_RESULTS_H=list()
ROB_RESULTS$rob_high_ran=NA
for ( Ri in 1:nrow(ROB_RESULTS)){
    print(names(NETS)[Ri])
    net=NETS[[Ri]]
    rr_results=numeric()
    for (IT in 1:1000){
      if(IT%%100==0){print(IT)}
      OR_ran=sample(1:ncol(net))
      SE=second.extinct(net, participant = "higher", method="external", ext.col = OR_ran)
      rr_results[IT]=robustness(SE)
    }
    RANDOMROB_RESULTS_H[[Ri]]=rr_results
    ROB_RESULTS$rob_high_ran[Ri]=mean(rr_results)
} 
# Only compound networks (1000 simulations)
# Random order of modules
# Randomly removes species from the first module until the module is fully excluded, then moves to the next and so on
# COMPOUND1ROB_RESULTS_H= robustness for all simulations
# ROB_RESULTS$rob_high_mod= means
COMPOUND1ROB_RESULTS_H=list()
ROB_RESULTS$rob_high_mod=NA
for ( Ri in 1:nrow(ROB_RESULTS)){
  if(!is.na(ROB_RESULTS$topology[Ri])){
    if(ROB_RESULTS$topology[Ri]=="compound"){
      print(names(NETS)[Ri])
      net=NETS[[Ri]]
      rr_results=numeric()
      mod1=MODULES[[Ri]]$cols
      for (IT in 1:1000){
        if(IT%%100==0){print(IT)}
        OR_ran=numeric()
        n2=1:ncol(net)
        n3=rep(TRUE,ncol(net))
        mod2=sample(unique(mod1))
        MO=1
        for (n1 in n2){
          n3=(mod1==mod2[MO]&!is.element(n2,OR_ran))
          n4=ifelse(sum(n3)==1,n2[n3],sample(n2[n3],1))
          OR_ran[n1]=n4
          if(sum(n3)==1){MO=MO+1}
        }
        SE=second.extinct(net, participant = "higher", method="external", ext.col = OR_ran)
        rr_results[IT]=robustness(SE)
      }
      COMPOUND1ROB_RESULTS_H[[Ri]]=rr_results
      ROB_RESULTS$rob_high_mod[Ri]=mean(rr_results)
    }
  }
}
save(ROB_RESULTS,RANDOMROB_RESULTS_L,COMPOUND1ROB_RESULTS_L,COMPOUND2ROB_RESULTS_L,RANDOMROB_RESULTS_H,COMPOUND1ROB_RESULTS_H, file=paste("results/robustness_",NAME_OUTPUT,".RData",sep=""))
# Only compound networks (1000 simulations)
# Avoid excluding modules
# In each round, randonly removes a species from the largest module(s)
# COMPOUND2ROB_RESULTS_H= robustness for all simulations
# ROB_RESULTS$rob_high_mod2= means
COMPOUND2ROB_RESULTS_H=list()
ROB_RESULTS$rob_high_mod2=NA
for ( Ri in 1:nrow(ROB_RESULTS)){
  if(!is.na(ROB_RESULTS$topology[Ri])){
    if(ROB_RESULTS$topology[Ri]=="compound"){
      print(names(NETS)[Ri])
      net=NETS[[Ri]]
      rr_results=numeric()
      mod1=MODULES[[Ri]]$cols
      for (IT in 1:1000){
        if(IT%%100==0){print(IT)}
        OR_ran=numeric()
        n2=1:ncol(net)
        n3=rep(TRUE,ncol(net))
        for (n1 in n2){
          n3=(!is.element(n2,OR_ran))
          MAX1=max(table(mod1[n3]))
          mod2=names(table(mod1[n3]))[table(mod1[n3])==MAX1]
          n3=(is.element(mod1,mod2)&!is.element(n2,OR_ran))
          n4=ifelse(sum(n3)==1,n2[n3],sample(n2[n3],1))
          OR_ran[n1]=n4
        }
        SE=second.extinct(net, participant = "higher", method="external", ext.col = OR_ran)
        rr_results[IT]=robustness(SE)
      }
      COMPOUND2ROB_RESULTS_H[[Ri]]=rr_results
      ROB_RESULTS$rob_high_mod2[Ri]=mean(rr_results)
    }
  }
}
save(ROB_RESULTS,RANDOMROB_RESULTS_L,COMPOUND1ROB_RESULTS_L,COMPOUND2ROB_RESULTS_L,RANDOMROB_RESULTS_H,COMPOUND1ROB_RESULTS_H,COMPOUND2ROB_RESULTS_H, file=paste("results/robustness_",NAME_OUTPUT,".RData",sep=""))