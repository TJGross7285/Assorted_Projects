{\rtf1\ansi\ansicpg1252\cocoartf2511
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;\red0\green0\blue0;}
{\*\expandedcolortbl;;\csgenericrgb\c0\c0\c0;}
\margl1440\margr1440\vieww11840\viewh8360\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 ############################ (4/11) Evaluate Simon's 24-Metabolite Findings \
#################Generate 24-metabolite data set for RAS\
setwd("~/Desktop/Simon_Findings")\
\
####Load in entire RAS data set metabolite space\
load("diff.mets.Rochester.RData")\
\
####Generate Discovery data set containing only metabolites in the 24-Metabolite panel\
discovery_transposed<-as.data.frame(t(mat.discovery))\
simon_ras_discovery_24<-dplyr::select(discovery_transposed,`PC ae C40:6`,`PC aa C40:1`,`PC aa C38:6`,`PC aa C38:0`,`PC aa C38:0`,`PC aa C36:6`,`lysoPC a C18:2`, `C3`,`PC ae C36:4`,`C10:2`,`C9`,`PC ae C42:1`,`PC aa C38:3`,`C5`,`ADMA`,`Asn`,`PC aa C34:4`,`C18:1-OH`, `PC ae C34:0`,`C5-OH (C3-DC-M)`,`PC aa C40:5`,`PC aa C32:0`,`C16:2`,`C12:1`,`C10:1`)\
matx1<-t(simon_ras_discovery_24)\
\
####Generate Validation data set containing only metabolites in the 24-Metabolite panel\
validation_transposed<-as.data.frame(t(mat.validation))\
simon_ras_validation_24<-dplyr::select(validation_transposed,`PC ae C40:6`,`PC aa C40:1`,`PC aa C38:6`,`PC aa C38:0`,`PC aa C38:0`,`PC aa C36:6`,`lysoPC a C18:2`, `C3`,`PC ae C36:4`,`C10:2`,`C9`,`PC ae C42:1`,`PC aa C38:3`,`C5`,`ADMA`,`Asn`,`PC aa C34:4`,`C18:1-OH`, `PC ae C34:0`,`C5-OH (C3-DC-M)`,`PC aa C40:5`,`PC aa C32:0`,`C16:2`,`C12:1`,`C10:1`)\
matx2<-t(simon_ras_validation_24)\
\
####Rename Disease State factor vector to be consistent with "metab10.RData" object \
gid1<-gid.discovery\
gid2<-gid.validation\
\
####Save RData object\
save(matx1,matx2,gid1,gid2, file="metab24New_4_11.RData")\
\
############################\
############################*****IN NEW R SESSION*****\
\
###################Evaluate Simon's Code for the 24-Metabolite Model\
\
setwd("~/Desktop/Simon_Findings")\
library(hlr)\
library(pROC)\
\
load("metab24New_4_11.RData")\
\
\pard\tx528\tx1056\tx1584\tx2112\tx2640\tx3168\tx3696\tx4224\tx4752\tx5280\tx5808\tx6337\tx6865\tx7393\tx7921\tx8449\tx8977\tx9505\tx10033\tx10561\tx11089\tx11617\tx12145\tx12674\tx13202\tx13730\tx14258\tx14786\tx15314\tx15842\tx16370\tx16898\tx17426\tx17954\tx18483\tx19011\tx19539\tx20067\tx20595\tx21123\tx21651\tx22179\tx22707\tx23235\tx23763\tx24291\tx24820\tx25348\tx25876\tx26404\tx26932\tx27460\tx27988\tx28516\tx29044\tx29572\tx30100\tx30628\tx31157\tx31685\tx32213\tx32741\tx33269\tx33797\li80\fi-80\pardirnatural\partightenfactor0
\cf2 #gid1 is the group id for the patients in the discovery set, 1: Control, 3: Converter-pre\
idx <- c((1:length(gid1))[gid1==1],(1:length(gid1))[gid1==3])\
\
#metabolite level matrix for the patients in the discovery set, row names of the matrix are the metabolite names\
matx1 <- matx1[,idx]\
\
dptmp <- gid1[idx]\
dptmp1 <- ifelse(dptmp==1,0,1)\
\
#gid2 is the group id for the patients in the validation set, 1: Control, 3: Converter-pre, row names of the matrix are the metabolite names\
idx <- c((1:length(gid2))[gid2==1],(1:length(gid2))[gid2==3])\
\
#metabolite level matrix for the patients in the validation set\
matx2 <- matx2[,idx]\
\
dptmp <- gid2[idx]\
dptmp2 <- ifelse(dptmp==1,0,1)\
\
#####*****WEMEL APPEARS TO BE USING RANDOM NUMBER GENERATION; SEED INSERTED BY T.G. (OTHERWISE SIMON'S CODE)*****\
set.seed(122)\
out <- WEMEL(t(matx1),t(matx1),dptmp1,delta=0.001)\
phat1 <- coef(out[[2]])[1]+t(matx1)%*%matrix(coef(out[[2]])[-1],nrow(matx1),1)\
phat_new = exp(phat1)/(1+exp(phat1))\
dy <- factor(dptmp1)\
\
pdf("roc24Confirmation.pdf",paper="letter",width=8,height=8,bg="transparent")\
\
\
rocobj.nc.conv.d<- plot.roc(dy, phat_new,ylab="True positive rate", xlab="False positive rate", main="ROC curve for the 24 metabolites for the discovery set\\n NC vs CONVERTERpre", percent=TRUE,  ci=TRUE,print.auc=TRUE,legacy.axes=TRUE)\
\
ciobj <- ci.se(rocobj.nc.conv.d, specificities=seq(0, 100, 5))\
plot(ciobj, type="shape", col="#1c61b6AA")\
\
phat1 <- coef(out[[2]])[1]+t(matx2)%*%matrix(coef(out[[2]])[-1],nrow(matx1),1)\
phat_new = exp(phat1)/(1+exp(phat1))\
dy <- factor(dptmp2)\
\
rocobj.nc.conv.v<- plot.roc(dy, phat_new,ylab="True positive rate", xlab="False positive rate", main="ROC curve for the 24 metabolites for the validation set\\n NC vs CONVERTERpre", percent=0.16,  ci=TRUE,print.auc=TRUE,legacy.axes=TRUE)\
\
ciobj <- ci.se(rocobj.nc.conv.v, specificities=seq(0, 100, 5))\
plot(ciobj, type="shape", col="#1c61b6AA")\
\
dev.off()\
}