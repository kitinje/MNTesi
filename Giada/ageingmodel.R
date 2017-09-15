# AGING PROJECT

setwd("~/Desktop/MRC/AGE/")
dir()

load("preselection.set")
load("preTRAIN.danes_bmiq.bval")
load("hm450")

rownames(hm450) = hm450$probe_id

t2  = which(hm450$channel=='Both') 
t1g  = which(hm450$channel=='Grn')
t1r = which(hm450$channel=='Red')

preTRAIN = preTRAIN.danes_bmiq.bval[rownames(hm450),rownames(preselection.set)]

NEURO.qn = preTRAIN[,preselection.set$BATCH=="GSE41826_NEURONS"]
boxplot(NEURO.qn)

t1g.goldstd <- preprocessCore::normalize.quantiles.determine.target(x=NEURO.qn[t1g,])
t1r.goldstd <- preprocessCore::normalize.quantiles.determine.target(x=NEURO.qn[t1r,])
t2.goldstd <- preprocessCore::normalize.quantiles.determine.target(x=NEURO.qn[t2,])

levels(factor(preselection.set$BATCH))

NEURO.qn = preTRAIN[,preselection.set$BATCH=="GSE41826_NEURONS"]
NEURO.qn[t1g,] = preprocessCore::normalize.quantiles.use.target(NEURO.qn[t1g,], target=t1g.goldstd)
NEURO.qn[t1r,] = preprocessCore::normalize.quantiles.use.target(NEURO.qn[t1r,], target=t1r.goldstd)
NEURO.qn[t2,] = preprocessCore::normalize.quantiles.use.target(NEURO.qn[t2,], target=t2.goldstd)

BONE.qn = preTRAIN[,preselection.set$BATCH=="BONE"]
BONE.qn[t1g,] = preprocessCore::normalize.quantiles.use.target(BONE.qn[t1g,], target=t1g.goldstd)
BONE.qn[t1r,] = preprocessCore::normalize.quantiles.use.target(BONE.qn[t1r,], target=t1r.goldstd)
BONE.qn[t2,] = preprocessCore::normalize.quantiles.use.target(BONE.qn[t2,], target=t2.goldstd)

EMTAB1866.qn = preTRAIN[,preselection.set$BATCH=="EMTAB1866"]
EMTAB1866.qn[t1g,] = preprocessCore::normalize.quantiles.use.target(EMTAB1866.qn[t1g,], target=t1g.goldstd)
EMTAB1866.qn[t1r,] = preprocessCore::normalize.quantiles.use.target(EMTAB1866.qn[t1r,], target=t1r.goldstd)
EMTAB1866.qn[t2,] = preprocessCore::normalize.quantiles.use.target(EMTAB1866.qn[t2,], target=t2.goldstd)

GSE41826_GLIA.qn = preTRAIN[,preselection.set$BATCH=="GSE41826_GLIA"]
GSE41826_GLIA.qn[t1g,] = preprocessCore::normalize.quantiles.use.target(GSE41826_GLIA.qn[t1g,], target=t1g.goldstd)
GSE41826_GLIA.qn[t1r,] = preprocessCore::normalize.quantiles.use.target(GSE41826_GLIA.qn[t1r,], target=t1r.goldstd)
GSE41826_GLIA.qn[t2,] = preprocessCore::normalize.quantiles.use.target(GSE41826_GLIA.qn[t2,], target=t2.goldstd)

GSE64495.qn = preTRAIN[,preselection.set$BATCH=="GSE64495"]
GSE64495.qn[t1g,] = preprocessCore::normalize.quantiles.use.target(GSE64495.qn[t1g,], target=t1g.goldstd)
GSE64495.qn[t1r,] = preprocessCore::normalize.quantiles.use.target(GSE64495.qn[t1r,], target=t1r.goldstd)
GSE64495.qn[t2,] = preprocessCore::normalize.quantiles.use.target(GSE64495.qn[t2,], target=t2.goldstd)

TCGA_BLCA.qn = preTRAIN[,preselection.set$BATCH=="TCGA_BLCA"]
TCGA_BLCA.qn[t1g,] = preprocessCore::normalize.quantiles.use.target(TCGA_BLCA.qn[t1g,], target=t1g.goldstd)
TCGA_BLCA.qn[t1r,] = preprocessCore::normalize.quantiles.use.target(TCGA_BLCA.qn[t1r,], target=t1r.goldstd)
TCGA_BLCA.qn[t2,] = preprocessCore::normalize.quantiles.use.target(TCGA_BLCA.qn[t2,], target=t2.goldstd)

TCGA_BRCA.qn = preTRAIN[,preselection.set$BATCH=="TCGA_BRCA"]
TCGA_BRCA.qn[t1g,] = preprocessCore::normalize.quantiles.use.target(TCGA_BRCA[t1g,], target=t1g.goldstd)
TCGA_BRCA.qn[t1r,] = preprocessCore::normalize.quantiles.use.target(TCGA_BRCA[t1r,], target=t1r.goldstd)
TCGA_BRCA.qn[t2,] = preprocessCore::normalize.quantiles.use.target(TCGA_BRCA[t2,], target=t2.goldstd)

TCGA_COAD.qn = preTRAIN[,preselection.set$BATCH=="TCGA_COAD"]
TCGA_COAD.qn[t1g,] = preprocessCore::normalize.quantiles.use.target(TCGA_COAD.qn[t1g,], target=t1g.goldstd)
TCGA_COAD.qn[t1r,] = preprocessCore::normalize.quantiles.use.target(TCGA_COAD.qn[t1r,], target=t1r.goldstd)
TCGA_COAD.qn[t2,] = preprocessCore::normalize.quantiles.use.target(TCGA_COAD.qn[t2,], target=t2.goldstd)

TCGA_HNSC.qn = preTRAIN[,preselection.set$BATCH=="TCGA_HNSC"]
TCGA_HNSC.qn[t1g,] = preprocessCore::normalize.quantiles.use.target(TCGA_HNSC.qn[t1g,], target=t1g.goldstd)
TCGA_HNSC.qn[t1r,] = preprocessCore::normalize.quantiles.use.target(TCGA_HNSC.qn[t1r,], target=t1r.goldstd)
TCGA_HNSC.qn[t2,] = preprocessCore::normalize.quantiles.use.target(TCGA_HNSC.qn[t2,], target=t2.goldstd)

TCGA_KIPAN.qn = preTRAIN[,preselection.set$BATCH=="TCGA_KIPAN"]
TCGA_KIPAN.qn[t1g,] = preprocessCore::normalize.quantiles.use.target(TCGA_KIPAN.qn[t1g,], target=t1g.goldstd)
TCGA_KIPAN.qn[t1r,] = preprocessCore::normalize.quantiles.use.target(TCGA_KIPAN.qn[t1r,], target=t1r.goldstd)
TCGA_KIPAN.qn[t2,] = preprocessCore::normalize.quantiles.use.target(TCGA_KIPAN.qn[t2,], target=t2.goldstd)

TCGA_LIHC.qn = preTRAIN[,preselection.set$BATCH=="TCGA_LIHC"]
TCGA_LIHC.qn[t1g,] = preprocessCore::normalize.quantiles.use.target(TCGA_LIHC.qn[t1g,], target=t1g.goldstd)
TCGA_LIHC.qn[t1r,] = preprocessCore::normalize.quantiles.use.target(TCGA_LIHC.qn[t1r,], target=t1r.goldstd)
TCGA_LIHC.qn[t2,] = preprocessCore::normalize.quantiles.use.target(TCGA_LIHC.qn[t2,], target=t2.goldstd)

TCGA_LUPAN.qn = preTRAIN[,preselection.set$BATCH=="TCGA_LUPAN"]
TCGA_LUPAN.qn[t1g,] = preprocessCore::normalize.quantiles.use.target(TCGA_LUPAN.qn[t1g,], target=t1g.goldstd)
TCGA_LUPAN.qn[t1r,] = preprocessCore::normalize.quantiles.use.target(TCGA_LUPAN.qn[t1r,], target=t1r.goldstd)
TCGA_LUPAN.qn[t2,] = preprocessCore::normalize.quantiles.use.target(TCGA_LUPAN.qn[t2,], target=t2.goldstd)

TCGA_PRAD.qn = preTRAIN[,preselection.set$BATCH=="TCGA_PRAD"]
TCGA_PRAD.qn[t1g,] = preprocessCore::normalize.quantiles.use.target(TCGA_PRAD.qn[t1g,], target=t1g.goldstd)
TCGA_PRAD.qn[t1r,] = preprocessCore::normalize.quantiles.use.target(TCGA_PRAD.qn[t1r,], target=t1r.goldstd)
TCGA_PRAD.qn[t2,] = preprocessCore::normalize.quantiles.use.target(TCGA_PRAD.qn[t2,], target=t2.goldstd)

TCGA_THCA.qn = preTRAIN[,preselection.set$BATCH=="TCGA_THCA"]
TCGA_THCA.qn[t1g,] = preprocessCore::normalize.quantiles.use.target(TCGA_THCA.qn[t1g,], target=t1g.goldstd)
TCGA_THCA.qn[t1r,] = preprocessCore::normalize.quantiles.use.target(TCGA_THCA.qn[t1r,], target=t1r.goldstd)
TCGA_THCA.qn[t2,] = preprocessCore::normalize.quantiles.use.target(TCGA_THCA.qn[t2,], target=t2.goldstd)

TCGA_UCEC.qn = preTRAIN[,preselection.set$BATCH=="TCGA_UCEC"]
TCGA_UCEC.qn[t1g,] = preprocessCore::normalize.quantiles.use.target(TCGA_UCEC.qn[t1g,], target=t1g.goldstd)
TCGA_UCEC.qn[t1r,] = preprocessCore::normalize.quantiles.use.target(TCGA_UCEC.qn[t1r,], target=t1r.goldstd)
TCGA_UCEC.qn[t2,] = preprocessCore::normalize.quantiles.use.target(TCGA_UCEC.qn[t2,], target=t2.goldstd)

preTRAIN.qn = cbind(BONE.qn, EMTAB1866.qn, GSE41826_GLIA.qn, NEURO.qn, GSE64495.qn, TCGA_BLCA.qn, TCGA_BRCA.qn, TCGA_COAD.qn, TCGA_HNSC.qn, TCGA_KIPAN.qn, TCGA_LIHC.qn,TCGA_LUPAN.qn, TCGA_PRAD.qn, TCGA_THCA.qn, TCGA_UCEC.qn)
preTRAIN.qn = preTRAIN.qn[, rownames(preselection.set)]

boxplot(preTRAIN.qn)

length(which(is.na(colSums(preTRAIN.qn))))

library(sva)
TRAIN <-sva::ComBat(dat=preTRAIN.qn, batch=as.factor(preselection.set$BATCH), mod=cbind(preselection.set$AGE)) #2918

plot(preTRAIN.danes_bmiq.bval["cg16867657",],preselection.set$AGE) # 0.7081259
plot(preTRAIN.qn["cg16867657",],preselection.set$AGE) # 0.7037765
plot(TRAIN["cg16867657",],preselection.set$AGE) # 0.8623394


***********************************************LHO-FSELECTION

Train.Ages<-preselection.set

train_size=dim(TRAIN)[2]
lho_size=train_size/2
count=dim(TRAIN)[1]

set.seed(1)  ; LHO_REP1  = sample(c(1:train_size),size=lho_size,replace = FALSE)
set.seed(2)  ; LHO_REP2  = sample(c(1:train_size),size=lho_size,replace = FALSE)
set.seed(3)  ; LHO_REP3  = sample(c(1:train_size),size=lho_size,replace = FALSE)
set.seed(4)  ; LHO_REP4  = sample(c(1:train_size),size=lho_size,replace = FALSE)
set.seed(5)  ; LHO_REP5  = sample(c(1:train_size),size=lho_size,replace = FALSE)
set.seed(6)  ; LHO_REP6  = sample(c(1:train_size),size=lho_size,replace = FALSE)
set.seed(7)  ; LHO_REP7  = sample(c(1:train_size),size=lho_size,replace = FALSE)
set.seed(8)  ; LHO_REP8  = sample(c(1:train_size),size=lho_size,replace = FALSE)
set.seed(9)  ; LHO_REP9  = sample(c(1:train_size),size=lho_size,replace = FALSE)
set.seed(10) ; LHO_REP10 = sample(c(1:train_size),size=lho_size,replace = FALSE)

TRAIN_R1 = TRAIN[,LHO_REP1]
TRAIN_S1 = Train.Ages[colnames(TRAIN_R1),]
SpearCorr_R1<-apply(TRAIN_R1, 1, cor.test, TRAIN_S1$AGE, method="spearman")
rho_S1<-as.double(unlist(SpearCorr_R1)[grep("estimate.rho",names(unlist(SpearCorr_R1)))])
names(rho_S1)<-gsub(pattern=".estimate.rho",x=rownames(TRAIN_R1), replacement="")
pvalue_S1<-as.double(unlist(SpearCorr_R1)[grep("p.value",names(unlist(SpearCorr_R1)))])
names(pvalue_S1)<-gsub(pattern=".p.value",x=rownames(TRAIN_R1), replacement="")

TRAIN_R2 = TRAIN[,LHO_REP2]
TRAIN_S2 = Train.Ages[colnames(TRAIN_R2),]
SpearCorr_R2<-apply(TRAIN_R2, 1, cor.test, TRAIN_S2$AGE, method="spearman")
rho_S2<-as.double(unlist(SpearCorr_R2)[grep("estimate.rho",names(unlist(SpearCorr_R2)))])
names(rho_S2)<-gsub(pattern=".estimate.rho",x=rownames(TRAIN_R2), replacement="")
pvalue_S2<-as.double(unlist(SpearCorr_R2)[grep("p.value",names(unlist(SpearCorr_R2)))])
names(pvalue_S2)<-gsub(pattern=".p.value",x=rownames(TRAIN_R2), replacement="")

TRAIN_R3 = TRAIN[,LHO_REP3]
TRAIN_S3 = Train.Ages[colnames(TRAIN_R3),]
SpearCorr_R3<-apply(TRAIN_R3, 1, cor.test, TRAIN_S3$AGE, method="spearman")
rho_S3<-as.double(unlist(SpearCorr_R3)[grep("estimate.rho",names(unlist(SpearCorr_R3)))])
names(rho_S3)<-gsub(pattern=".estimate.rho",x=rownames(TRAIN_R3), replacement="")
pvalue_S3<-as.double(unlist(SpearCorr_R3)[grep("p.value",names(unlist(SpearCorr_R3)))])
names(pvalue_S3)<-gsub(pattern=".p.value",x=rownames(TRAIN_R3), replacement="")

TRAIN_R4 = TRAIN[,LHO_REP4]
TRAIN_S4 = Train.Ages[colnames(TRAIN_R4),]
SpearCorr_R4<-apply(TRAIN_R4, 1, cor.test, TRAIN_S4$AGE, method="spearman")
rho_S4<-as.double(unlist(SpearCorr_R4)[grep("estimate.rho",names(unlist(SpearCorr_R4)))])
names(rho_S4)<-gsub(pattern=".estimate.rho",x=rownames(TRAIN_R4), replacement="")
pvalue_S4<-as.double(unlist(SpearCorr_R4)[grep("p.value",names(unlist(SpearCorr_R4)))])
names(pvalue_S4)<-gsub(pattern=".p.value",x=rownames(TRAIN_R4), replacement="")

TRAIN_R5 = TRAIN[,LHO_REP5]
TRAIN_S5 = Train.Ages[colnames(TRAIN_R5),]
SpearCorr_R5<-apply(TRAIN_R5, 1, cor.test, TRAIN_S5$AGE, method="spearman")
rho_S5<-as.double(unlist(SpearCorr_R5)[grep("estimate.rho",names(unlist(SpearCorr_R5)))])
names(rho_S5)<-gsub(pattern=".estimate.rho",x=rownames(TRAIN_R5), replacement="")
pvalue_S5<-as.double(unlist(SpearCorr_R5)[grep("p.value",names(unlist(SpearCorr_R5)))])
names(pvalue_S5)<-gsub(pattern=".p.value",x=rownames(TRAIN_R5), replacement="")

TRAIN_R6 = TRAIN[,LHO_REP6]
TRAIN_S6 = Train.Ages[colnames(TRAIN_R6),]
SpearCorr_R6<-apply(TRAIN_R6, 1, cor.test, TRAIN_S6$AGE, method="spearman")
rho_S6<-as.double(unlist(SpearCorr_R6)[grep("estimate.rho",names(unlist(SpearCorr_R6)))])
names(rho_S6)<-gsub(pattern=".estimate.rho",x=rownames(TRAIN_R6), replacement="")
pvalue_S6<-as.double(unlist(SpearCorr_R6)[grep("p.value",names(unlist(SpearCorr_R6)))])
names(pvalue_S6)<-gsub(pattern=".p.value",x=rownames(TRAIN_R6), replacement="")

TRAIN_R7 = TRAIN[,LHO_REP7]
TRAIN_S7 = Train.Ages[colnames(TRAIN_R7),]
SpearCorr_R7<-apply(TRAIN_R7, 1, cor.test, TRAIN_S7$AGE, method="spearman")
rho_S7<-as.double(unlist(SpearCorr_R7)[grep("estimate.rho",names(unlist(SpearCorr_R7)))])
names(rho_S7)<-gsub(pattern=".estimate.rho",x=rownames(TRAIN_R7), replacement="")
pvalue_S7<-as.double(unlist(SpearCorr_R7)[grep("p.value",names(unlist(SpearCorr_R7)))])
names(pvalue_S7)<-gsub(pattern=".p.value",x=rownames(TRAIN_R7), replacement="")

TRAIN_R8 = TRAIN[,LHO_REP8]
TRAIN_S8 = Train.Ages[colnames(TRAIN_R8),]
SpearCorr_R8<-apply(TRAIN_R8, 1, cor.test, TRAIN_S8$AGE, method="spearman")
rho_S8<-as.double(unlist(SpearCorr_R8)[grep("estimate.rho",names(unlist(SpearCorr_R8)))])
names(rho_S8)<-gsub(pattern=".estimate.rho",x=rownames(TRAIN_R8), replacement="")
pvalue_S8<-as.double(unlist(SpearCorr_R8)[grep("p.value",names(unlist(SpearCorr_R8)))])
names(pvalue_S8)<-gsub(pattern=".p.value",x=rownames(TRAIN_R8), replacement="")

TRAIN_R9 = TRAIN[,LHO_REP9]
TRAIN_S9 = Train.Ages[colnames(TRAIN_R9),]
SpearCorr_R9<-apply(TRAIN_R9, 1, cor.test, TRAIN_S9$AGE, method="spearman")
rho_S9<-as.double(unlist(SpearCorr_R9)[grep("estimate.rho",names(unlist(SpearCorr_R9)))])
names(rho_S9)<-gsub(pattern=".estimate.rho",x=rownames(TRAIN_R9), replacement="")
pvalue_S9<-as.double(unlist(SpearCorr_R9)[grep("p.value",names(unlist(SpearCorr_R9)))])
names(pvalue_S9)<-gsub(pattern=".p.value",x=rownames(TRAIN_R9), replacement="")

TRAIN_R10 = TRAIN[,LHO_REP10]
TRAIN_S10 = Train.Ages[colnames(TRAIN_R10),]
SpearCorr_R10<-apply(TRAIN_R10, 1, cor.test, TRAIN_S10$AGE, method="spearman")
rho_S10<-as.double(unlist(SpearCorr_R10)[grep("estimate.rho",names(unlist(SpearCorr_R10)))])
names(rho_S10)<-gsub(pattern=".estimate.rho",x=rownames(TRAIN_R10), replacement="")
pvalue_S10<-as.double(unlist(SpearCorr_R10)[grep("p.value",names(unlist(SpearCorr_R10)))])
names(pvalue_S10)<-gsub(pattern=".p.value",x=rownames(TRAIN_R10), replacement="")

save(LHO_REP10,file="LHO_REP10"); save(TRAIN_S10,file="TRAIN_S10");save(SpearCorr_R10,file="SpearCorr_R10");save(rho_S10,file="rho_S10");save(pvalue_S10,file="pvalue_S10")
save(LHO_REP9,file="LHO_REP9"); save(TRAIN_S10,file="TRAIN_S9");save(SpearCorr_R9,file="SpearCorr_R9");save(rho_S9,file="rho_S9");save(pvalue_S9,file="pvalue_S9")
save(LHO_REP8,file="LHO_REP8"); save(TRAIN_S10,file="TRAIN_S8");save(SpearCorr_R8,file="SpearCorr_R8");save(rho_S8,file="rho_S8");save(pvalue_S8,file="pvalue_S8")
save(LHO_REP7,file="LHO_REP7"); save(TRAIN_S10,file="TRAIN_S7");save(SpearCorr_R7,file="SpearCorr_R7");save(rho_S7,file="rho_S7");save(pvalue_S7,file="pvalue_S7")
save(LHO_REP6,file="LHO_REP6"); save(TRAIN_S10,file="TRAIN_S6");save(SpearCorr_R6,file="SpearCorr_R6");save(rho_S6,file="rho_S6");save(pvalue_S6,file="pvalue_S6")
save(LHO_REP5,file="LHO_REP5"); save(TRAIN_S10,file="TRAIN_S5");save(SpearCorr_R5,file="SpearCorr_R5");save(rho_S5,file="rho_S5");save(pvalue_S5,file="pvalue_S5")
save(LHO_REP4,file="LHO_REP4"); save(TRAIN_S10,file="TRAIN_S4");save(SpearCorr_R4,file="SpearCorr_R4");save(rho_S4,file="rho_S4");save(pvalue_S4,file="pvalue_S4")
save(LHO_REP3,file="LHO_REP3"); save(TRAIN_S10,file="TRAIN_S3");save(SpearCorr_R3,file="SpearCorr_R3");save(rho_S3,file="rho_S3");save(pvalue_S3,file="pvalue_S3")
save(LHO_REP2,file="LHO_REP2"); save(TRAIN_S10,file="TRAIN_S2");save(SpearCorr_R2,file="SpearCorr_R2");save(rho_S2,file="rho_S2");save(pvalue_S2,file="pvalue_S2")
save(LHO_REP1,file="LHO_REP1"); save(TRAIN_S10,file="TRAIN_S1");save(SpearCorr_R1,file="SpearCorr_R1");save(rho_S1,file="rho_S1");save(pvalue_S1,file="pvalue_S1")

rho_LHO<-cbind(rho_S1,rho_S2,rho_S3,rho_S4,rho_S5,rho_S6,rho_S7,rho_S8,rho_S9,rho_S10)
colnames(rho_LHO)<-c("REP1","REP2","REP3","REP4","REP5","REP6","REP7","REP8","REP9","REP10")

pvalue_LHO<-cbind(pvalue_S1,pvalue_S2,pvalue_S3,pvalue_S4,pvalue_S5,pvalue_S6,pvalue_S7,pvalue_S8,pvalue_S9,pvalue_S10)
colnames(pvalue_LHO)<-c("REP1","REP2","REP3","REP4","REP5","REP6","REP7","REP8","REP9","REP10")

meanrho_LHO<-rowMeans(rho_LHO)
meanpvalue_LHO<-rowMeans(pvalue_LHO)
mean_rhopval<-cbind(meanpvalue_LHO,meanrho_LHO)
colnames(mean_rhopval)<-c("pvalue","rho")

# mean_rhopval[which(mean_rhopval[,1]<0.05),] 29916
mean_rhopval = cbind(mean_rhopval, mean_rhopval[,1]*nrow(TRAIN))
colnames(mean_rhopval)[3]="padj"
mean_rhopval[which(mean_rhopval[,3]>1),3]<-1


library(h2o)
localH2O = h2o.init(nthreads = 6)

# MODELLO 1

probes_m001 = mean_rhopval[which(mean_rhopval[,1]<0.000001),]
probes_m001 = probes_m001[order(probes_m001[,1]),]

TRAIN_m001<-TRAIN[rownames(probes_m001),rownames(preselection.set)]
dim(TRAIN_m001); dim(preselection.set)

train_m001 <- rbind(preselection.set$AGE,TRAIN_m001); rownames(train_m001)[1]="AGE"
train_m001 = as.h2o(t(train_m001), destination_frame="train_m001"); dim(train_m001)

( model_m001 = h2o.glm(y = 1, x = 2:735, training_frame = train_m001,  family = "gaussian", alpha = 0, nlambda=100, nfolds=10, missing_values_handling="MeanImputation"))


# TEST EMTAB1866

load("EMTAB1866.danes.bvals")
load("EMTAB1866.Ages")

EMTAB1866.Ages<-EMTAB1866.Ages[which(!(rownames(EMTAB1866.Ages) %in% rownames(preselection.set))),]
EMTAB1866.qn = EMTAB1866.danes.bvals[rownames(TRAIN),rownames(EMTAB1866.Ages)]

EMTAB1866.qn[t1g,] = preprocessCore::normalize.quantiles.use.target(EMTAB1866.qn[t1g,], target=t1g.goldstd)
EMTAB1866.qn[t1r,] = preprocessCore::normalize.quantiles.use.target(EMTAB1866.qn[t1r,], target=t1r.goldstd)
EMTAB1866.qn[t2,] = preprocessCore::normalize.quantiles.use.target(EMTAB1866.qn[t2,], target=t2.goldstd)

EMTAB1866.test = EMTAB1866.qn[rownames(probes_m001),]
EMTAB1866.test = as.h2o(t(EMTAB1866.test), destination_frame="EMTAB1866.test"); dim(EMTAB1866.test)

EMTAB1866.preds = h2o.predict(model_m001, EMTAB1866.test)
EMTAB1866.test.preds<-as.data.frame(TEST.preds); rownames(TEST.preds)=colnames(TEST)

# The END
