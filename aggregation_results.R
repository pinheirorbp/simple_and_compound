# This script is called by the RMarkdown file ("binary.Rmd" in the folder "results"). 
#However, they might be used separated. In this case, define the codes for the subsets for aggregation in an object CODES
# CODES= c("set1","set2","set3")
### packages ####
library(bipartite)
library(stringr)
### Loading results ####
L= length(CODES)
load(paste("files/modularity_",CODES[1],".RData",sep=""))
load(paste("files/nestedness_",CODES[1],".RData",sep=""))
load(paste("files/compound_",CODES[1],".RData",sep=""))
rm(COMPOUNDcomplete)
### significance and topology ####
nodf_equi_sig=NODFsummary$NODF>NODFsummary$equi_ci_up
mod_equi_sig=MODsummary$modularity>MODsummary$equi_null_ci_up
nodf_equi_SM_sig=COMPOUNDsummary$NODF_SM>COMPOUNDsummary$equi_SM_ci_up
nodf_equi_DM_sig=COMPOUNDsummary$NODF_DM>COMPOUNDsummary$equi_DM_ci_up
nodf_prop_sig=NODFsummary$NODF>NODFsummary$prop_ci_up
mod_prop_sig=MODsummary$modularity>MODsummary$prop_null_ci_up
nodf_prop_SM_sig=COMPOUNDsummary$NODF_SM>COMPOUNDsummary$prop_SM_ci_up
nodf_prop_DM_sig=COMPOUNDsummary$NODF_DM>COMPOUNDsummary$prop_DM_ci_up
# topology: equiprobable
topology_equi=rep(NA,nrow(MODsummary))
topology_equi[!mod_equi_sig&nodf_equi_sig]="nested"
topology_equi[mod_equi_sig&nodf_equi_SM_sig]="compound"
topology_equi[mod_equi_sig&!nodf_equi_SM_sig]="pure modular"
topology_equi[!mod_equi_sig&!nodf_equi_sig]="unstructured"
# topology: proportional
topology_prop=rep(NA,nrow(MODsummary))
topology_prop[!mod_prop_sig&nodf_prop_sig]="nested"
topology_prop[mod_prop_sig&nodf_prop_SM_sig]="compound"
topology_prop[mod_prop_sig&!nodf_prop_SM_sig]="pure modular"
topology_prop[!mod_prop_sig&!nodf_prop_sig]="unstructured"
# topology: proportional for modularity, equiprobable for nestedness
topology=rep(NA,nrow(MODsummary))
topology[!mod_prop_sig&nodf_equi_sig]="nested"
topology[mod_prop_sig&nodf_equi_SM_sig]="compound"
topology[mod_prop_sig&!nodf_equi_SM_sig]="pure modular"
topology[!mod_prop_sig&!nodf_equi_sig]="unstructured"
### loading networks ####
X= read.table(paste("IDS_",CODES[1],".txt",sep=""), header = T, sep="\t")
ids=X$ID
L1=length(ids)
NETS=list()
for (i in 1:L1){
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
rm(Cs,Rs)
### dimensions and connectance ####
X$rows=NA
X$cols=NA
X$connectance=NA
X$nint=NA
for (i in 1:L1){
  X$rows[i]=nrow(NETS[[ids[i]]])
  X$cols[i]=ncol(NETS[[ids[i]]])
  X$nint[i]=sum(NETS[[ids[i]]])
  X$connectance[i]=X$nint[i]/(X$rows[i]*X$cols[i])
}
linkage_density=(X$nint/(X$rows+X$cols))
min_dim=X$rows>=MIN&X$cols>=MIN
### distortions in null models ####
mean_equi_nullconnec=numeric()
mean_equi_nullconnec_P=numeric()
sd_equi_nullconnec=numeric()
for (i in 1:L1){
  if(!is.null(NODFcomplete[[i]])){
    connec=numeric()
    connec2=numeric()
    for (j in 1:105){
      connec[j]=sum(NODFcomplete[[i]]$equi_nulls[[j]]$NullMatrix)/(nrow(NODFcomplete[[i]]$equi_nulls[[j]]$NullMatrix)*ncol(NODFcomplete[[i]]$equi_nulls[[j]]$NullMatrix))
      connec2[j]=connec[j]/X$connectance[i]
    }
    mean_equi_nullconnec[i]=mean(connec, na.rm=T)
    sd_equi_nullconnec[i]=sd(connec, na.rm=T)
    mean_equi_nullconnec_P[i]=mean(connec2, na.rm=T)
  }
}
mean_prop_nullconnec=numeric()
mean_prop_nullconnec_P=numeric()
sd_prop_nullconnec=numeric()
for (i in 1:L1){
  if(!is.null(NODFcomplete[[i]])){
    connec=numeric()
    connec2=numeric()
    for (j in 1:105){
      connec[j]=sum(NODFcomplete[[i]]$prop_nulls[[j]]$NullMatrix)/(nrow(NODFcomplete[[i]]$prop_nulls[[j]]$NullMatrix)*ncol(NODFcomplete[[i]]$prop_nulls[[j]]$NullMatrix))
      connec2[j]=connec[j]/X$connectance[i]
    }
    mean_prop_nullconnec[i]=mean(connec, na.rm=T)
    sd_prop_nullconnec[i]=sd(connec, na.rm=T)
    mean_prop_nullconnec_P[i]=mean(connec2, na.rm=T)
  }
}
### plot nested ####
for (i in 1:L1){
  par(mar=c(4.2,2,2,2))
  png(filename = paste("plots/nest_matrix/binary/",ids[i],".png",sep = ""),width = 14, height = nrow(NETS[[i]])/ncol(NETS[[i]])*10+6.2, units = "cm", res=300)
  plotmatrix(sortmatrix(as.matrix(NETS[[i]]),topology = "nested"), binary=T,xlab = names(NETS)[i])
  if(!min_dim[i]){
    mtext(paste("[Less than",MIN,"species in each dimension]"), side=1, line=1, col="red", cex=1)}
  if(min_dim[i]&linkage_density[i]<1){
    mtext("[Linkage density <1]", side=1, line=1, col="red", cex=1)
  }
  dev.off()
  dev.off()
}
### plot compound topology ####
ISMOD=MODsummary$modularity>MODsummary$prop_null_ci_up
ISMOD[is.na(ISMOD)]=FALSE
for (i in 1:L1){
  if(ISMOD[i]&linkage_density[i]>=1&min_dim[i]){
    par(mar=c(4.2,2,2,2))
    png(filename = paste("plots/compound_matrix/binary/",ids[i],".png",sep = ""),width = 14, height = nrow(NETS[[i]])/ncol(NETS[[i]])*10+6.2, units = "cm", res=300)
    plotmatrix(sortmatrix(as.matrix(NETS[[i]]),topology = "compound",row_partitions = MODcomplete[[i]]$row_modules, col_partitions = MODcomplete[[i]]$col_modules,mod_similarity = T), binary=T,xlab = names(NETS)[i],border = T)
    dev.off()
    dev.off()
  }
}
#####
NETDATA=X
TABLE_RESULTS=cbind(topology, topology_equi, topology_prop,rows=X$rows,cols=X$cols,nint=X$nint,connectance=X$connectance,interaction_type=X$Interaction, NODFsummary, MODsummary, COMPOUNDsummary, nodf_equi_sig, mod_equi_sig, nodf_equi_SM_sig, nodf_equi_DM_sig, nodf_prop_sig, mod_prop_sig, nodf_prop_SM_sig, nodf_prop_DM_sig,mean_equi_nullconnec,mean_equi_nullconnec_P,sd_equi_nullconnec, mean_prop_nullconnec, mean_prop_nullconnec_P,sd_prop_nullconnec)
rownames(TABLE_RESULTS)=NETDATA$ID
rm(MODcomplete,MODsummary,NODFcomplete,NODFsummary,COMPOUNDsummary)

if(L>1){
  for(ii in 2:L){
    load(paste("files/modularity_",CODES[ii],".RData",sep=""))
    load(paste("files/nestedness_",CODES[ii],".RData",sep=""))
    load(paste("files/compound_",CODES[ii],".RData",sep=""))
    rm(COMPOUNDcomplete)
    ### significance and topology ####
    nodf_equi_sig=NODFsummary$NODF>NODFsummary$equi_ci_up
    mod_equi_sig=MODsummary$modularity>MODsummary$equi_null_ci_up
    nodf_equi_SM_sig=COMPOUNDsummary$NODF_SM>COMPOUNDsummary$equi_SM_ci_up
    nodf_equi_DM_sig=COMPOUNDsummary$NODF_DM>COMPOUNDsummary$equi_DM_ci_up
    nodf_prop_sig=NODFsummary$NODF>NODFsummary$prop_ci_up
    mod_prop_sig=MODsummary$modularity>MODsummary$prop_null_ci_up
    nodf_prop_SM_sig=COMPOUNDsummary$NODF_SM>COMPOUNDsummary$prop_SM_ci_up
    nodf_prop_DM_sig=COMPOUNDsummary$NODF_DM>COMPOUNDsummary$prop_DM_ci_up
    # topology: equiprobable
    topology_equi=rep(NA,nrow(MODsummary))
    topology_equi[!mod_equi_sig&nodf_equi_sig]="nested"
    topology_equi[mod_equi_sig&nodf_equi_SM_sig]="compound"
    topology_equi[mod_equi_sig&!nodf_equi_SM_sig]="pure modular"
    topology_equi[!mod_equi_sig&!nodf_equi_sig]="other"
    # topology: proportional
    topology_prop=rep(NA,nrow(MODsummary))
    topology_prop[!mod_prop_sig&nodf_prop_sig]="nested"
    topology_prop[mod_prop_sig&nodf_prop_SM_sig]="compound"
    topology_prop[mod_prop_sig&!nodf_prop_SM_sig]="pure modular"
    topology_prop[!mod_prop_sig&!nodf_prop_sig]="other"
    # topology: proportional for modularity, equiprobable for nestedness
    topology=rep(NA,nrow(MODsummary))
    topology[!mod_prop_sig&nodf_equi_sig]="nested"
    topology[mod_prop_sig&nodf_equi_SM_sig]="compound"
    topology[mod_prop_sig&!nodf_equi_SM_sig]="pure modular"
    topology[!mod_prop_sig&!nodf_equi_sig]="other"
    ### loading networks ####
    X= read.table(paste("IDS_",CODES[ii],".txt",sep=""), header = T, sep="\t")
    ids=X$ID
    L1=length(ids)
    NETS=list()
    for (i in 1:L1){
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
    rm(Cs,Rs)
    X$rows=NA
    X$cols=NA
    X$connectance=NA
    X$nint
    for (i in 1:L1){
      X$rows[i]=nrow(NETS[[ids[i]]])
      X$cols[i]=ncol(NETS[[ids[i]]])
      X$nint[i]=sum(NETS[[ids[i]]])
      X$connectance[i]=X$nint[i]/(X$rows[i]*X$cols[i])
    }
    linkage_density=X$nint/(X$rows+X$cols)
    min_dim=X$rows>=MIN&X$cols>=MIN
    ### distortions in null models ####
    mean_equi_nullconnec=numeric()
    mean_equi_nullconnec_P=numeric()
    sd_equi_nullconnec=numeric()
    for (i in 1:L1){
      if(!is.null(NODFcomplete[[i]])){
        connec=numeric()
        connec2=numeric()
        for (j in 1:105){
          connec[j]=sum(NODFcomplete[[i]]$equi_nulls[[j]]$NullMatrix)/(nrow(NODFcomplete[[i]]$equi_nulls[[j]]$NullMatrix)*ncol(NODFcomplete[[i]]$equi_nulls[[j]]$NullMatrix))
          connec2[j]=connec[j]/X$connectance[i]
        }
        mean_equi_nullconnec[i]=mean(connec, na.rm=T)
        sd_equi_nullconnec[i]=sd(connec, na.rm=T)
        mean_equi_nullconnec_P[i]=mean(connec2, na.rm=T)
      }
    }
    mean_prop_nullconnec=numeric()
    mean_prop_nullconnec_P=numeric()
    sd_prop_nullconnec=numeric()
    for (i in 1:L1){
      if(!is.null(NODFcomplete[[i]])){
        connec=numeric()
        connec2=numeric()
        for (j in 1:105){
          connec[j]=sum(NODFcomplete[[i]]$prop_nulls[[j]]$NullMatrix)/(nrow(NODFcomplete[[i]]$prop_nulls[[j]]$NullMatrix)*ncol(NODFcomplete[[i]]$prop_nulls[[j]]$NullMatrix))
          connec2[j]=connec[j]/X$connectance[i]
        }
        mean_prop_nullconnec[i]=mean(connec, na.rm=T)
        sd_prop_nullconnec[i]=sd(connec, na.rm=T)
        mean_prop_nullconnec_P[i]=mean(connec2, na.rm=T)
      }
    }
    ### plot nested ####
    for (i in 1:L1){
      par(mar=c(4.2,2,2,2))
      png(filename = paste("plots/nest_matrix/binary/",ids[i],".png",sep = ""),width = 14, height = nrow(NETS[[i]])/ncol(NETS[[i]])*10+6.2, units = "cm", res=300)
      plotmatrix(sortmatrix(as.matrix(NETS[[i]]),topology = "nested"), binary=T,xlab = names(NETS)[i])
      if(!min_dim[i]){
        mtext(paste("Less than",MIN,"species in each dimension"), side=1, line=1, col="red", cex=.5)}
      if(min_dim[i]&linkage_density[i]<1){
        mtext("Linkage density <1", side=1, line=1, col="red", cex=.5)
      }
      dev.off()
      dev.off()
    }
    ### plot compound topology ####
    ISMOD=(MODsummary$modularity>MODsummary$prop_null_ci_up)
    ISMOD[is.na(ISMOD)]=FALSE
    for (i in 1:L1){
      if(ISMOD[i]&linkage_density[i]>=1&min_dim[i]){
        par(mar=c(4.2,2,2,2))
        png(filename = paste("plots/compound_matrix/binary/",ids[i],".png",sep = ""),width = 14, height = nrow(NETS[[i]])/ncol(NETS[[i]])*10+6.2, units = "cm", res=300)
        plotmatrix(sortmatrix(as.matrix(NETS[[i]]),topology = "compound",row_partitions = MODcomplete[[i]]$row_modules, col_partitions = MODcomplete[[i]]$col_modules,mod_similarity = T), binary=T,xlab = names(NETS)[i],border = T)
        dev.off()
        dev.off()
      }
    }
    # binding tables
    TAB1=cbind(topology, topology_equi, topology_prop,rows=X$rows,cols=X$cols,nint=X$nint,connectance=X$connectance,interaction_type=X$Interaction, NODFsummary, MODsummary, COMPOUNDsummary, nodf_equi_sig, mod_equi_sig, nodf_equi_SM_sig, nodf_equi_DM_sig, nodf_prop_sig, mod_prop_sig, nodf_prop_SM_sig, nodf_prop_DM_sig,mean_equi_nullconnec,mean_equi_nullconnec_P,sd_equi_nullconnec, mean_prop_nullconnec, mean_prop_nullconnec_P,sd_prop_nullconnec)
    rm(MODcomplete,MODsummary,NODFcomplete,NODFsummary,COMPOUNDsummary)
    TABLE_RESULTS=rbind(TABLE_RESULTS,TAB1)
    NETDATA=rbind(NETDATA,X)
    rownames(TABLE_RESULTS)=NETDATA$ID
  }}
write.table(TABLE_RESULTS,file =paste("results/",NAME_OUTPUT,".txt",sep=""),sep = "\t",row.names = T)
write.table(NETDATA,file =paste("results/","NETDATA_",NAME_OUTPUT,".txt",sep=""),sep = "\t",row.names = T)
rm(list = as.character(ls())[!is.element(as.character(ls()),c("CODES","NAME_OUTPUT","MIN"))])