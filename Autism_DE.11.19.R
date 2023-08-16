########################################### autism 



pheno_sheet<-read.csv("Autism_Pheno_7_27.csv",check.names=FALSE)%>%as.data.frame()
colnames(pheno_sheet)[1]<-"SubjectID"
abunds_sheet<-read.csv("Autism_Metab_7_27.csv",check.names=FALSE,na.strings=c("0",".","NA","na"))%>%filter(method!="GCMS")%>%select(-c(1:5))%>%log2()%>%t()%>%as.data.frame()
abunds_sheet<-cbind(rownames(abunds_sheet), abunds_sheet)
colnames(abunds_sheet)[1]<-"SubjectID"
join<-dplyr::inner_join(pheno_sheet, abunds_sheet, by="SubjectID")
colnames(join)[4:6]<-c("SpecDx","BroadDx","Sex")

abunds_sheet<-read.csv("Autism_Metab_7_27.csv",check.names=FALSE,na.strings=c("0",".","NA","na"))%>%filter(method!="GCMS")%>%select(c(1:5))
lookup<-cbind(colnames(join)[10:dim(join)[2]],abunds_sheet)
colnames(lookup)[1]<-"FeatureID"

####Pull concentration data, log2 transform, and impute NA/Inf values 
abunds<-join[,-c(1:9)]
pheno<-join[,c(1:9)]
abunds[is.na(abunds)==TRUE]<-NA 

library(impute)
T_abunds<-t(abunds)
edata<-impute.knn(T_abunds,k = 10, rowmax = 1.0, colmax = 1.0, maxp = 1500, rng.seed=362436069)$data
rownames(edata)<-colnames(abunds)

####Calculate SVs 
library(sva)
mod<-model.matrix(~as.factor(SpecDx)+as.factor(Batch)+as.numeric(Age)+as.numeric(IQ)+as.factor(Sex),data=pheno)
mod0<-model.matrix(~as.factor(Batch)+as.numeric(Age)+as.numeric(IQ)+as.factor(Sex),data=pheno)
n.sv_BE<-num.sv(edata,mod,seed=122)
svobj_BE<-sva(edata,mod,mod0,n.sv=n.sv_BE)
n.sv_LEEK<-num.sv(edata,mod,seed=122,method="leek")
svobj_LEEK<-sva(edata,mod,mod0,n.sv=n.sv_LEEK)

design_table<-cbind(as.data.frame(pheno),as.data.frame(svobj_BE$sv))
colnames(design_table)[10:25]<-c("SV1","SV2","SV3","SV4","SV5","SV6","SV7","SV8","SV9","SV10","SV11","SV12","SV13","SV14","SV15","SV16")

####Set up accessory objects for DE incorporating pre-post timepoints
design1<-model.matrix(~0+SpecDx+SV1+SV2+SV3+SV4+SV5+SV6+SV7+SV8+SV9+SV10+SV11+SV12+SV13+SV14+SV15+SV16,data=design_table)
colnames(design1)[1:3]<-c("HFA",
						  "LFA",
						  "TYP")

####Fit linear model for pairwise contrasts 
arrayw<-arrayWeights(edata, design=design1)
fit1<-lmFit(edata,design1,weights=arrayw)
cm1 <- makeContrasts(`HFA-LFA` = HFA-LFA,levels=design1)
cm2<- makeContrasts(`HFA-TYP` = HFA-TYP, levels=design1)
cm3<- makeContrasts(`LFA-TYP` = LFA-TYP, levels=design1)
	
fit1_F1 <- contrasts.fit(fit1, cm1)
fit1_F1 <- eBayes(fit1_F1,trend=TRUE)

fit2_F2 <- contrasts.fit(fit1, cm2)
fit2_F2 <- eBayes(fit2_F2,trend=TRUE)

fit3_F3 <- contrasts.fit(fit1, cm3)
fit3_F3 <- eBayes(fit3_F3,trend=TRUE)


####Output final DE table and write to file
T.1<-topTableF(fit1_F1,adjust="BH",number=100000)
T.1<-cbind(rownames(T.1),T.1)
colnames(T.1)[1]<-"FeatureID"

T.2<-topTableF(fit2_F2,adjust="BH",number=100000)
T.2<-cbind(rownames(T.2),T.2)
colnames(T.2)[1]<-"FeatureID"

T.3<-topTableF(fit3_F3,adjust="BH",number=100000)
T.3<-cbind(rownames(T.3),T.3)
colnames(T.3)[1]<-"FeatureID"


final_HFA_LFA<-dplyr::inner_join(lookup,T.1,by="FeatureID")
final_HFA_TYP<-dplyr::inner_join(lookup,T.2,by="FeatureID")
final_LFA_TYP<-dplyr::inner_join(lookup,T.3,by="FeatureID")



write.table(final_HFA_LFA%>%filter(ESImode=="ESIpos")%>%select("mz","RT","P.Value","HFA.LFA"),
			sep="\t",file="ASD_POS_HFA_LFA.txt",row.names=FALSE)
write.table(final_HFA_LFA%>%filter(ESImode=="ESIneg")%>%select("mz","RT","P.Value","HFA.LFA"),
			sep="\t",file="ASD_NEG_HFA_LFA.txt",row.names=FALSE)

write.table(final_HFA_TYP%>%filter(ESImode=="ESIpos")%>%select("mz","RT","P.Value","HFA.TYP"),
			sep="\t",file="ASD_POS_HFA_TYP.txt",row.names=FALSE)
write.table(final_HFA_TYP%>%filter(ESImode=="ESIneg")%>%select("mz","RT","P.Value","HFA.TYP"),
			sep="\t",file="ASD_NEG_HFA_TYP.txt",row.names=FALSE)

write.table(final_LFA_TYP%>%filter(ESImode=="ESIpos")%>%select("mz","RT","P.Value","LFA.TYP"),
			sep="\t",file="ASD_POS_LFA_TYP.txt",row.names=FALSE)
write.table(final_LFA_TYP%>%filter(ESImode=="ESIneg")%>%select("mz","RT","P.Value","LFA.TYP"),
			sep="\t",file="ASD_NEG_LFA_TYP.txt",row.names=FALSE)


mummichog -f ASD_POS_HFA_LFA.txt -o ASD_POS_HFA_LFA -m positive 
mummichog -f ASD_NEG_HFA_LFA.txt -o ASD_NEG_HFA_LFA -m negative  

mummichog -f ASD_POS_HFA_TYP.txt -o ASD_POS_HFA_TYP -m positive  
mummichog -f ASD_NEG_HFA_TYP.txt -o ASD_NEG_HFA_TYP -m negative 

mummichog -f ASD_POS_LFA_TYP.txt -o ASD_POS_LFA_TYP -m positive 
mummichog -f ASD_NEG_LFA_TYP.txt -o ASD_NEG_LFA_TYP -m negative 









model formula 
