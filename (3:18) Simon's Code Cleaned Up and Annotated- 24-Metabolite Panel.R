{\rtf1\ansi\ansicpg1252\cocoartf2511
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww10800\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 ######################(3/18) Simon's Code Cleaned Up and Annotated: 24-Metabolite Panel\
library(pROC)\
library(hlr)\
\
####Load 24-Metabolite RAS Data\
load("metab24.RData")\
\
####Extract only Discovery data set samples from total 24-Metabolite data set\
### gid1 is the group id for the patients in the discovery set, 1: Control; 3: Converter\
idx_Discovery <- c((1:length(gid1))[gid1==1],(1:length(gid1))[gid1==3])\
\
####Generate metabolite matrix for Discovery portion of total data set\
###metabolite level matrix for the patients in the discovery set, row names of the matrix are the metabolite names\
matx1_discovery <- matx1[,idx_Discovery]\
\
####Generate disease state factor for samples in Discovery portion of total data set \
dptmp_Discovery<- gid1[idx_Discovery]\
dptmp1 <- ifelse(dptmp_Discovery==1,0,1)\
\
####Evaluate WEMEL model in HLR on Discovery data w/ delta constant fixed at 0.001\
24_WEMEL <- WEMEL(t(matx1_discovery), t(matx1_discovery), dptmp1, delta=0.001)\
\
####Derive logits and probabilities for samples in Discovery sample\
phat1_discovery <- coef(24_WEMEL[[2]])[1]+t(matx1_discovery)%*%matrix(coef(24_WEMEL[[2]])[-1],nrow(matx1_discovery),1)\
phat_new_discovery<- as.vector(exp(phat1_discovery)/(1+exp(phat1_discovery)))\
\
####Set up factor vector for Discovery data \
dy_discovery<-as.factor(factor(dptmp1))\
levels(dy_discovery)<-c("Control", "Case")\
\
####Evaluate WEMEL ROC statistics and plot ROC curve for Discovery data\
discovery_log_roc<-roc(dy_discovery, phat_new_discovery, ci=TRUE, auc=TRUE)\
rocobj_discovery<- plot.roc(dy_discovery, phat_new_discovery, ylab="True Positive Rate", xlab="False Positive Rate", main="24-Metabolite Panel: Discovery Sample (HLR)\\n Control vs. Converter", percent=TRUE,  ci=TRUE, print.auc=TRUE, legacy.axes=TRUE)\
ciobj_discovery <- ci.se(rocobj_discovery, specificities=seq(0, 100, .5))\
plot(ciobj_discovery, type="shape", col="#1c61b6AA")\
\
####Set up Caret object\
library(caret)\
cctrl1 <- trainControl(method = "repeatedcv", repeats=100, number = 10, classProbs = TRUE, summaryFunction = twoClassSummary)\
\
####Evaluate 24-metabolite panel in Discovery sample using svmRadial\
set.seed(122)\
svm_discovery_metabolites<-as.data.frame(t(matx1_discovery))\
svm.24metab<-train(svm_discovery_metabolites, dy_discovery, method="svmRadial", family="binomial", trControl = cctrl1, metric = "ROC")\
predict.svm.discovery.24<-predict(svm.24metab, t(matx1_discovery), type="prob")\
\
\
####Evaluate rbfSVM ROC statistics and plot ROC curve for Discovery data\
roc.svm.discovery.24<-roc(dy_discovery, predict.svm.discovery.24$Case, ci=TRUE, auc=TRUE)\
rocobj_discovery_svm<- plot.roc(dy_discovery, predict.svm.discovery.24$Case, ylab="True Positive Rate", xlab="False Positive Rate", main="24-Metabolite Panel (rbfSVM): Discovery Sample\\n Control vs. Converter", percent=TRUE,  ci=TRUE, print.auc=TRUE, legacy.axes=TRUE)\
ciobj_discovery_svm<- ci.se(rocobj_discovery_svm, specificities=seq(0, 100, .5))\
plot(ciobj_discovery_svm, type="shape", col="#1c61b6AA")\
\
####Extract only Validation data set samples from total data set\
###gid2 is the group id for the patients in the validation set, 1: Control, 3: Converter-pre, row names of the matrix are the metabolite names\
idx_Validation <- c((1:length(gid2))[gid2==1],(1:length(gid2))[gid2==3])\
\
####Generate metabolite matrix for Validation portion of total data set\
###metabolite level matrix for the patients in the validation set\
matx2 <- matx2[,idx_Validation]\
\
####Generate disease state factor for samples in Validation portion of total data set \
dptmp_Validation<- gid2[idx_Validation]\
dptmp2 <- ifelse(dptmp_Validation==1,0,1)\
\
####Evaluate fit of WEMEL model in held-out validation sample\
phat1_validation<- coef(24_WEMEL[[2]])[1]+t(matx2)%*%matrix(coef(24_WEMEL[[2]])[-1],nrow(matx1),1)\
phat_new_validation <- as.vector(exp(phat1_validation)/(1+exp(phat1_validation)))\
\
####Set up factor vector for Validation data \
dy_validation <- as.factor(factor(dptmp2))\
levels(dy_validation)<-c("Control","Case")\
\
####Evaluate fit of rbfSVM model in held-out Validation sample \
svm_validation_metabolites<-as.data.frame(t(matx2))\
predict.svm.validation.24<-predict(svm.24metab, svm_validation_metabolites, type="prob")\
roc.svm.validation.24<-roc(dy_validation, predict.svm.validation.24$Case, ci=TRUE, auc=TRUE)\
\
####Evaluate WEMEL ROC statistics and plot ROC curve for Validation data\
validation_roc_WEMEL<-roc(dy_validation, phat_new_validation, ci=TRUE, auc=TRUE)\
rocobj_validation_WEMEL<- plot.roc(dy_validation, phat_new_validation, ylab="True Positive Rate", xlab="False Positive Rate", main="24-Metabolite Panel (HLR): Validation Sample\\n Control vs. Converter", percent=TRUE,  ci=TRUE, print.auc=TRUE, legacy.axes=TRUE)\
ciobj_validation_WEMEL<- ci.se(rocobj_validation_WEMEL, specificities=seq(0, 100, .5))\
plot(ciobj_validation_WEMEL, type="shape", col="#1c61b6AA")\
\
####Evaluate rbfSVM ROC statistics and plot ROC curve for Validation data\
roc.svm.validation.24<-roc(dy_validation, predict.svm.validation.24$Case, ci=TRUE, auc=TRUE)\
rocobj_validation_svm<- plot.roc(dy_validation, predict.svm.validation.24$Case, ylab="True Positive Rate", xlab="False Positive Rate", main="24-Metabolite Panel (rbfSVM): Validation Sample\\n Control vs. Converter", percent=TRUE,  ci=TRUE, print.auc=TRUE, legacy.axes=TRUE)\
ciobj_validation_svm<- ci.se(rocobj_validation_svm, specificities=seq(0, 100, .5))\
plot(ciobj_validation_svm, type="shape", col="#1c61b6AA")\
\
####Save 24-Metabolite WEMEL model and as an RData object\
save(24_WEMEL, svm.24metab, file="Simons24MetabModel.RData")\
}