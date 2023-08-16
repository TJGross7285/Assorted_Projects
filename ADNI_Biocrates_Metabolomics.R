####################################(4/25/2020) Biocrates ADNI
library(dplyr)
library(impute)

###consolidate/write ADNI pheno data to csv *****DO NOT OVERWRITE***** 
pheno<-read.csv("DXSUM_PDXCONV_ADNIALL.csv",check.names=FALSE)[,c(1:6,10:11)]%>%filter(VISCODE=="bl"|VISCODE2=="bl")%>%filter(Phase!="ADNI3")
write.csv(pheno,file="fix_pheno.csv")

######################
######################
######################


###Read in ADNI Biocrates Data
Z<-read.csv("ADMCDUKEP180UPLCADNI2GO.csv",na.strings="",check.names=FALSE) 
Y<-read.csv("ADMCDUKEP180FIAADNI2GO.csv",na.strings="",check.names=FALSE)
all.equal(Z$RID,Y$RID)

###Filter abundances by LOD  
total<-cbind(Z[-c(67)],Y[,-c(1:25,167)])
index<-apply(apply(total,2,is.na),2,sum)/dim(total)[1]
use<-total[,index<.333]

###Join clinical data and abundances 
pheno<-read.csv("fix_pheno.csv",check.names=FALSE)
total_join<-dplyr::inner_join(pheno,use,by="RID")%>%filter(DXMASTER!=4 & DXMASTER!=7)
write.csv(total_join,file="FOR_ANALYSIS_5_13.csv")

pheno<-as.numeric(total_join$DXMASTER)
levels(pheno)<-list("NC"=c("1"),"MCI"=c("2"),"AD"=c("3"))
edata<-t(total_join[,-c(1:33)])
edata<-t(impute.knn(edata,k = 10,rng.seed=362436069)$data)%>%as.data.frame()

######Cut Features:`C16.1.OH`
all_metabs_cross<-dplyr::transmute(edata,
	`lysoPC.a.C18.1`,
	`lysoPC.a.C18.2`,
	`PC.aa.C34.4`,
	`PC.aa.C36.6`,
	`PC.aa.C38.3`,
	`PC.aa.C40.5`,
	`PC.aa.C40.6`,
	`PC.aa.C34.4/lysoPC.a.C18.1`=`PC.aa.C34.4`/`lysoPC.a.C18.1`,
	`PC.aa.C34.4/lysoPC.a.C18.2`=`PC.aa.C34.4`/`lysoPC.a.C18.2`,
	`PC.aa.C36.5/lysoPC.a.C18.2`=`PC.aa.C36.5`/`lysoPC.a.C18.2`, 
	`PC.aa.C36.6/lysoPC.a.C18.1`=`PC.aa.C36.6`/`lysoPC.a.C18.1`,
	`PC.aa.C36.6/lysoPC.a.C18.2`=`PC.aa.C36.6`/`lysoPC.a.C18.2`,
	`C0`,
	`C2`,
	`C10`,
	`C12`,
	`C14.2`,
	`C14.1`,
	`C18.1`,
	`C16.1`,
	`C18`,
	`PC.ae.C36.2`,
	`PC.ae.C40.3`, 
	`PC.ae.C42.4`,
	`PC.ae.C44.4`,
	`SM..OH..C14.1`,
	`SM.C16.0`,
	`SM.C20.2`,
	`Val`,
	`PC.aa.C34.2`,
	`lysoPC.a.C18.0`,
	`Gly`,
	`lysoPC.a.C16.1`,
	`PC.ae.C42.3`,
	`C14.2.OH`,
	`PC.ae.C32.2`,
	`PC.ae.C34.0`,
	`PC.ae.C36.0`,
	`Phe`,
	`PC.ae.C34.1`,
	`PC.aa.C38.4`,
	`PC.aa.C42.2`,
	`PC.ae.C40.5`,
	`PC.ae.C38.5`,
	`Thr`)







#### Set up outcome vector and Caret control objects for classification
model<-cbind(pheno,edata)%>%filter(pheno!=2)
model_pheno<-as.factor(model$pheno)
levels(model_pheno)<-list("NC"=c("1"),"AD"=c("3"))
cctrl1<- trainControl(method = "repeatedcv", repeats = 10, number = 10, classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE)

set.seed(122)
test<-train(model[,-c(1)],model_pheno, method="rf", trControl = cctrl1, metric = "ROC")



A<-bcorsis(edata,as.numeric(pheno),method="interaction")
F<-as.numeric(A$ix)
pdf(file="ADNI.pdf")
for(i in 1:length(F)){
 	vioplot(edata[,F[i]]~pheno,main=colnames(edata[i]))
}
dev.off()