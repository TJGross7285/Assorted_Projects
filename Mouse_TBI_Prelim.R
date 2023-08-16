library(dplyr)

metab_pos<-read.csv("mapstone_a884_metabolomics_POS.csv",header=FALSE)
metab_pos_abunds<-as.data.frame(t(metab_pos[-c(1),-c(1:2)]))
colnames(metab_pos_abunds)<-paste(seq(1,dim(metab_pos_abunds)[2]),"POS_metab",sep=".")
metab_pos_IDs<-as.data.frame(t(metab_pos[1,-c(1:2)]))
colnames(metab_pos_IDs)<-"RodentID"

metab_neg<-read.csv("mapstone_a884_metabolomics_NEG.csv",header=FALSE)
metab_neg_abunds<-as.data.frame(t(metab_neg[-c(1),-c(1:2)]))
colnames(metab_neg_abunds)<-paste(seq(1,dim(metab_neg_abunds)[2]),"NEG_metab",sep=".")
metab_neg_IDs<-as.data.frame(t(metab_neg[1,-c(1:2)]))
colnames(metab_neg_IDs)<-"RodentID"

lipid_pos<-read.csv("mapstone_a884_lipidomics_POS.csv",header=FALSE)
lipid_pos_abunds<-as.data.frame(t(lipid_pos[-c(1),-c(1:2)]))
colnames(lipid_pos_abunds)<-paste(seq(1,dim(lipid_pos_abunds)[2]),"POS_lipid",sep=".")
lipid_pos_IDs<-as.data.frame(t(lipid_pos[1,-c(1:2)]))
colnames(lipid_pos_IDs)<-"RodentID"

lipid_neg<-read.csv("mapstone_a884_lipidomics_NEG.csv",header=FALSE)
lipid_neg_abunds<-as.data.frame(t(lipid_neg[-c(1),-c(1:2)]))
colnames(lipid_neg_abunds)<-paste(seq(1,dim(lipid_neg_abunds)[2]),"NEG_lipid",sep=".")
lipid_neg_IDs<-as.data.frame(t(lipid_neg[1,-c(1:2)]))
colnames(lipid_neg_IDs)<-"RodentID"

all.equal(metab_pos_IDs,metab_neg_IDs,lipid_pos_IDs,lipid_neg_IDs)

all<-cbind(metab_pos_IDs,metab_pos_abunds,metab_neg_abunds,lipid_pos_abunds,lipid_neg_abunds)
index<-colnames(all) %in% c("1359.NEG_lipid","1210.NEG_lipid","1236.NEG_metab","1088.NEG_lipid","1829.NEG_lipid","758.POS_lipid")


metab_pos_feature<-metab_pos[-c(1),1:2]
metab_neg_feature<-metab_neg[-c(1),1:2]
lipid_pos_feature<-lipid_pos[-c(1),1:2]
lipid_neg_feature<-lipid_neg[-c(1),1:2]
key<-cbind(as.data.frame(colnames(all[,-c(1)])),rbind(metab_pos_feature,metab_neg_feature,lipid_pos_feature,lipid_neg_feature))
colnames(key)<-c("FeatureID","M/Z","RT")
vec<-as.factor(c("C","C","C","C","A","B","D","C","D","B","D","A","B","A","D","B","A","D","D","B","B","A","A","C","B","A","D","C"))

subset<-cbind(metab_pos_IDs,all[,index==TRUE])%>%filter(metab_pos_IDs!="QC")
subset<-cbind(vec,subset)%>%write.csv(file="Mouse_TBI_USE.csv")