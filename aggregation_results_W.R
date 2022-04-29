# This script is called within the RMarkdown file ("weighted.Rmd" in the folder "results"). 
#However, they might be used separated. In this case, define the codes for the subsets for aggregation in an object CODES
# CODES= c("set1W")
### packages ####
library(bipartite)
library(stringr)
### Loading results ####
L= length(CODES)
load(paste("files/modularity_",CODES[1],".RData",sep=""))
load(paste("files/nestedness_",CODES[1],".RData",sep=""))
load(paste("files/compound_",CODES[1],".RData",sep=""))
### significance and topology ####
wnoda_equi_sig=WNODAsummary$WNODA>WNODAsummary$equi_ci_up
wmod_equi_sig=WMODsummary$modularity>WMODsummary$equi_null_ci_up
wnoda_equi_SM_sig=WCOMPOUNDsummary$WNODA_SM>WCOMPOUNDsummary$equi_SM_ci_up
wnoda_equi_DM_sig=WCOMPOUNDsummary$WNODA_DM>WCOMPOUNDsummary$equi_DM_ci_up
wnoda_prop_equal=WNODAsummary$WNODA>WNODAsummary$prop_ci_down2
wmod_prop_sig=WMODsummary$modularity>WMODsummary$prop_null_ci_up
wnoda_prop_SM_equal=WCOMPOUNDsummary$WNODA_SM>WCOMPOUNDsummary$prop_SM_ci_down2
wnoda_prop_DM_equal=WCOMPOUNDsummary$WNODA_DM>WCOMPOUNDsummary$prop_DM_ci_down2
# topology
topology=rep(NA,nrow(WMODsummary))
topology[!wmod_prop_sig&wnoda_equi_sig]="nested"
topology[wmod_prop_sig&wnoda_equi_SM_sig]="compound"
topology[wmod_prop_sig&!wnoda_equi_SM_sig]="pure modular"
topology[!wmod_prop_sig&!wnoda_equi_sig]="unstructured"
### loading networks ####
X= read.table(paste("IDS_",CODES[1],".txt",sep=""), header = T, sep="\t")
ids=X$ID
L1=length(ids)
# binary version
BIN_NETS=list()
for (i in 1:L1){
  if(str_detect(X$Binary_file[i],pattern = ".txt")){ #txt files
    BIN_NETS[[ids[i]]]<-read.table(paste("networks/",X$Binary_file[i],sep = ""))}
  if(str_detect(X$Binary_file[i],pattern = ".csv")){ #csv files
    BIN_NETS[[ids[i]]]<-read.table(paste("networks/",X$Binary_file[i],sep = ""), sep = ",")}
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
# weighted network
NETS=list()
X$WB_inconsistencies=""
for (i in 1:L1){
  if(str_detect(X$Weighted_file[i],pattern = ".txt")){ #txt files
    NETS[[ids[i]]]<-read.table(paste("networks/",X$Weighted_file[i],sep = ""))}
  if(str_detect(X$Weighted_file[i],pattern = ".csv")){ #csv files
    NETS[[ids[i]]]<-read.table(paste("networks/",X$Weighted_file[i],sep = ""), sep = ",")}
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
    print(paste(ids[[i]], "  Rows in binary matrix:",nrow(BIN_NETS[[i]])," Rows in weighted matrix: ",nrow(NETS[[i]])))
    X$WB_inconsistencies[i]=paste(X$WB_inconsistencies[i],"R")}
  if(ncol(BIN_NETS[[i]])!=ncol(NETS[[i]])){ # difference in cols
    print(paste(ids[[i]], "  Cols in binary matrix:",ncol(BIN_NETS[[i]])," Cols in weighted matrix: ",ncol(NETS[[i]])))
    X$WB_inconsistencies[i]=paste(X$WB_inconsistencies[i],"C")}
  conn_BIN_NET=round(sum(BIN_NETS[[i]])/(nrow(BIN_NETS[[i]])*ncol(BIN_NETS[[i]])),digits = 2)
  conn_bin2_NET=round(sum(bin2_NET)/(nrow(bin2_NET)*ncol(bin2_NET)),digits = 2)
  if(conn_BIN_NET!=conn_bin2_NET){ # difference in connectance 
    print(paste(ids[[i]],"  Connectance in binary matrix:",conn_BIN_NET,"Connectance in weighted matrix:",conn_bin2_NET))
    X$WB_inconsistencies[i]=paste(X$WB_inconsistencies[i],"connec")}
}
### dimensions and connectance ####
X$rows=NA
X$cols=NA
X$connectance=NA
X$nint=NA
X$Wint=NA
X$net_non_integer=NA
for (i in 1:L1){
  X$rows[i]=nrow(NETS[[ids[i]]])
  X$cols[i]=ncol(NETS[[ids[i]]])
  X$nint[i]=sum(NETS[[ids[i]]]>0)
  X$connectance[i]=X$nint[i]/(X$rows[i]*X$cols[i])
  X$Wint[i]=sum(NETS[[ids[i]]])
  X$net_non_integer[i]=any(NETS[[ids[i]]]%%1!=0)
}
min_dim=X$rows>=MIN&X$cols>=MIN
linkage_density=(X$nint/(X$rows+X$cols))
### plot nested ####
for (i in 1:L1){
  par(mar=c(4.2,2,2,2))
  png(filename = paste("plots/nest_matrix/weighted/",ids[i],".png",sep = ""),width = 14, height = nrow(NETS[[i]])/ncol(NETS[[i]])*10+6.2, units = "cm", res=300)
  plotmatrix(sortmatrix(as.matrix(NETS[[i]]),topology = "nested",sort_by = "weights"), binary=F,xlab = names(NETS)[i])
  if(!min_dim[i]){
    mtext(paste("[Less than",MIN,"species in each dimension]"), side=1, line=1, col="red", cex=1)}
  if(min_dim[i]&linkage_density[i]<1){
    mtext("[Linkage density <1]", side=1, line=1, col="red", cex=1)
  }
  dev.off()
  dev.off()
}
### plot compound topology ####
ISWMOD=WMODsummary$modularity>WMODsummary$prop_null_ci_up
for (i in 1:L1){
  if(ISWMOD[i]){
    par(mar=c(4.2,2,2,2))
    png(filename = paste("plots/compound_matrix/weighted/",ids[i],".png",sep = ""),width = 14, height = nrow(NETS[[i]])/ncol(NETS[[i]])*10+6.2, units = "cm", res=300)
    plotmatrix(sortmatrix(as.matrix(NETS[[i]]),topology = "compound",row_partitions = WMODcomplete[[i]]$row_modules, col_partitions = WMODcomplete[[i]]$col_modules,mod_similarity = T,sort_by = "weights"), binary=F,xlab = names(NETS)[i],border = T)
    dev.off()
    dev.off()
  }
}
NETDATA=X
TABLE_RESULTS=cbind(topology, rows=X$rows,cols=X$cols,nint=X$nint,Wint=X$Wint,WB_inconsistencies=X$WB_inconsistencies, net_non_interger=X$net_non_integer,connectance=X$connectance,interaction_type=X$Interaction, WNODAsummary, WMODsummary, WCOMPOUNDsummary, wnoda_equi_sig, wmod_equi_sig, wnoda_equi_SM_sig, wnoda_equi_DM_sig, wnoda_prop_equal, wmod_prop_sig, wnoda_prop_SM_equal,wnoda_prop_DM_equal)
rownames(TABLE_RESULTS)=NETDATA$ID
rm(WMODcomplete,WMODsummary,WNODAcomplete,WNODAsummary,WCOMPOUNDsummary,WCOMPOUNDcomplete)

if(L>1){
  for(ii in 2:L){
    load(paste("files/modularity_",CODES[ii],".RData",sep=""))
    load(paste("files/nestedness_",CODES[ii],".RData",sep=""))
    load(paste("files/compound_",CODES[ii],".RData",sep=""))
    ### significance and topology ####
    wnoda_equi_sig=WNODAsummary$WNODA>WNODAsummary$equi_ci_up
    wmod_equi_sig=WMODsummary$modularity>WMODsummary$equi_null_ci_up
    wnoda_equi_SM_sig=WCOMPOUNDsummary$WNODA_SM>WCOMPOUNDsummary$equi_SM_ci_up
    wnoda_equi_DM_sig=WCOMPOUNDsummary$WNODA_DM>WCOMPOUNDsummary$equi_DM_ci_up
    wnoda_prop_equal=WNODAsummary$WNODA>WNODAsummary$prop_ci_down2
    wmod_prop_sig=WMODsummary$modularity>WMODsummary$prop_null_ci_up
    wnoda_prop_SM_equal=WCOMPOUNDsummary$WNODA_SM>WCOMPOUNDsummary$prop_SM_ci_down2
    wnoda_prop_DM_equal=WCOMPOUNDsummary$WNODA_DM>WCOMPOUNDsummary$prop_DM_ci_down2
    # topology
    topology=rep(NA,nrow(WMODsummary))
    topology[!wmod_prop_sig&wnoda_equi_sig&wnoda_prop_equal]="perfectly nested"
    topology[!wmod_prop_sig&wnoda_equi_sig&!wnoda_prop_equal]="significant nestedness"
    topology[wmod_prop_sig&wnoda_equi_SM_sig&wnoda_prop_SM_equal]="perfectly compound"
    topology[wmod_prop_sig&wnoda_equi_SM_sig&!wnoda_prop_SM_equal]="significant compound"
    topology[wmod_prop_sig&!wnoda_equi_SM_sig]="pure modular"
    topology[!wmod_prop_sig&!wnoda_equi_sig]="other"
    ### loading networks ####
    X= read.table(paste("IDS_",CODES[1],".txt",sep=""), header = T, sep="\t")
    ids=X$ID
    L1=length(ids)
    # binary version
    BIN_NETS=list()
    for (i in 1:L1){
      if(str_detect(X$Binary_file[i],pattern = ".txt")){ #txt files
        BIN_NETS[[ids[i]]]<-read.table(paste("networks/",X$Binary_file[i],sep = ""))}
      if(str_detect(X$Binary_file[i],pattern = ".csv")){ #csv files
        BIN_NETS[[ids[i]]]<-read.table(paste("networks/",X$Binary_file[i],sep = ""), sep = ",")}
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
    # weighted network
    NETS=list()
    X$WB_inconsistencies=""
    for (i in 1:L1){
      if(str_detect(X$Weighted_file[i],pattern = ".txt")){ #txt files
        NETS[[ids[i]]]<-read.table(paste("networks/",X$Weighted_file[i],sep = ""))}
      if(str_detect(X$Weighted_file[i],pattern = ".csv")){ #csv files
        NETS[[ids[i]]]<-read.table(paste("networks/",X$Weighted_file[i],sep = ""), sep = ",")}
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
        print(paste(ids[[i]], "  Rows in binary matrix:",nrow(BIN_NETS[[i]])," Rows in weighted matrix: ",nrow(NETS[[i]])))
        X$WB_inconsistencies[i]=paste(X$WB_inconsistencies[i],"R")}
      if(ncol(BIN_NETS[[i]])!=ncol(NETS[[i]])){ # difference in cols
        print(paste(ids[[i]], "  Cols in binary matrix:",ncol(BIN_NETS[[i]])," Cols in weighted matrix: ",ncol(NETS[[i]])))
        X$WB_inconsistencies[i]=paste(X$WB_inconsistencies[i],"C")}
      conn_BIN_NET=round(sum(BIN_NETS[[i]])/(nrow(BIN_NETS[[i]])*ncol(BIN_NETS[[i]])),digits = 2)
      conn_bin2_NET=round(sum(bin2_NET)/(nrow(bin2_NET)*ncol(bin2_NET)),digits = 2)
      if(conn_BIN_NET!=conn_bin2_NET){ # difference in connectance 
        print(paste(ids[[i]],"  Connectance in binary matrix:",conn_BIN_NET,"Connectance in weighted matrix:",conn_bin2_NET))
        X$WB_inconsistencies[i]=paste(X$WB_inconsistencies[i],"connec")}
    }
    ### dimensions and connectance ####
    X$rows=NA
    X$cols=NA
    X$connectance=NA
    X$nint=NA
    X$Wint=NA
    X$net_non_integer=NA
    for (i in 1:L1){
      X$rows[i]=nrow(NETS[[ids[i]]])
      X$cols[i]=ncol(NETS[[ids[i]]])
      X$nint[i]=sum(NETS[[ids[i]]]>0)
      X$connectance[i]=X$nint[i]/(X$rows[i]*X$cols[i])
      X$Wint[i]=sum(NETS[[ids[i]]])
      X$net_non_integer[i]=any(NETS[[ids[i]]]%%1!=0)
    }
    min_dim=X$rows>=MIN&X$cols>=MIN
    linkage_density=(X$nint/(X$rows+X$cols))
    ### plot nested ####
    for (i in 1:L1){
      par(mar=c(4.2,2,2,2))
      png(filename = paste("plots/nest_matrix/weighted/",ids[i],".png",sep = ""),width = 14, height = nrow(NETS[[i]])/ncol(NETS[[i]])*10+6.2, units = "cm", res=300)
      plotmatrix(sortmatrix(as.matrix(NETS[[i]]),topology = "nested",sort_by = "weights"), binary=F,xlab = names(NETS)[i])
      if(!min_dim[i]){
        mtext(paste("[Less than",MIN,"species in each dimension]"), side=1, line=1, col="red", cex=1)}
      if(min_dim[i]&linkage_density[i]<1){
        mtext("[Linkage density <1]", side=1, line=1, col="red", cex=1)
      }
      dev.off()
      dev.off()
    }
    ### plot compound topology ####
    ISWMOD=WMODsummary$modularity>WMODsummary$prop_null_ci_up
    for (i in 1:L1){
      if(ISWMOD[i]){
        par(mar=c(4.2,2,2,2))
        png(filename = paste("plots/compound_matrix/weighted/",ids[i],".png",sep = ""),width = 14, height = nrow(NETS[[i]])/ncol(NETS[[i]])*10+6.2, units = "cm", res=300)
        plotmatrix(sortmatrix(as.matrix(NETS[[i]]),topology = "compound",row_partitions = WMODcomplete[[i]]$row_modules, col_partitions = WMODcomplete[[i]]$col_modules,mod_similarity = T,sort_by = "weights"), binary=F,xlab = names(NETS)[i],border = T)
        dev.off()
        dev.off()
      }
    }
    # binding tables
    TAB1=cbind(topology, rows=X$rows,cols=X$cols,nint=X$nint,Wint=X$Wint,WB_inconsistencies=X$WB_inconsistencies, net_non_interger=X$net_non_integer,connectance=X$connectance,interaction_type=X$Interaction, WNODAsummary, WMODsummary, WCOMPOUNDsummary, wnoda_equi_sig, wmod_equi_sig, wnoda_equi_SM_sig, wnoda_equi_DM_sig, wnoda_prop_equal, wmod_prop_sig, wnoda_prop_SM_equal,wnoda_prop_DM_equal)
    rm(WMODcomplete,WMODsummary,WNODAcomplete,WNODAsummary,WCOMPOUNDsummary,WCOMPOUNDcomplete)
    TABLE_RESULTS=rbind(TABLE_RESULTS,TAB1)
    NETDATA=rbind(NETDATA,X)
    rownames(TABLE_RESULTS)=NETDATA$ID
  }}
write.table(TABLE_RESULTS,file =paste("results/",NAME_OUTPUT,".txt",sep=""),sep = "\t",row.names = T)
write.table(NETDATA,file =paste("results/","NETDATA_",NAME_OUTPUT,".txt",sep=""),sep = "\t",row.names = T)
rm(list = as.character(ls())[!is.element(as.character(ls()),c("CODES","NAME_OUTPUT","MIN"))])