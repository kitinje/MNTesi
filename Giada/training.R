setwd("~/Desktop/DANIMARCA/GE/")

# GE
load("CCLE.RAW")
GE = CCLE.RAW
head(CCLE.RAW)


# BRCA
BRCA_ids = read.table("../ge_ids.txt", sep ="\t", row.names=1, stringsAsFactors = F, header = T)
#rownames(BRCA_ids)=BRCA_ids$CLINE
matchOrder = match(rownames(BRCA_ids), colnames(GE))
BRCA = GE[,matchOrder]
colnames(BRCA)=BRCA_ids$cell_line
rownames(BRCA_ids)=BRCA_ids$cell_line

y <- DGEList(counts=BRCA)
y <-y[rownames(y)[grep("_",rownames(y), invert = TRUE)],]
keep=grep(cbind(unlist(strsplit(rownames(y),'.',fixed=TRUE))), pattern = "ENSG")
rownames(y)<-cbind(unlist(strsplit(rownames(y),'.',fixed=TRUE)))[keep]
y <- calcNormFactors(y)

# keep files genes showing more than 0.5 cpm in more than 1 sample. 
y.raw = y
exp_thr=1
samp_thr=4
keep <- rowSums(cpm(y)>exp_thr) >= samp_thr
y<- y[keep, , keep.lib.sizes=FALSE] #ricalcola anche library sizes (consigliato)
y <- calcNormFactors(y)


xx <- as.list(org.Hs.egENSEMBL2EG)
ann <- cbind(unlist(xx)); colnames(ann) = "ENTREZID"
keep = cbind(ann[(rownames(ann) %in% rownames(y)),]); colnames(keep)="ENTREZID"
keep = cbind(keep[which(!duplicated(keep[,1])),]) ; colnames(keep)="ENTREZID"
y.Hs=y[rownames(keep),]
rownames(y.Hs)=keep[,1]
y.Hs$genes$GeneID=rownames(y.Hs)
entrez_size=dim(y.Hs)[1]
# Annotating Gene-Symbol
y.Hs$genes$Symbol <- mapIds(org.Hs.eg.db, rownames(y.Hs), keytype="ENTREZID", column="SYMBOL")

y.cpm = cpm(y.Hs, normalized.lib.sizes=TRUE, log=TRUE, prior.count=1)
rownames(y.cpm) = as.character(y.Hs$genes$Symbol)
BRCA.cpm = y.cpm

BRCA.cpm = BRCA.cpm[,rownames(BRCA_ids)]

TRAIN = BRCA.cpm[,BRCA_ids$GROUP2=="TRAIN"]
TEST = BRCA.cpm[,BRCA_ids$GROUP2=="TEST"]

TRAIN.score = BRCA_ids[colnames(TRAIN),]
TEST.score =  BRCA_ids[colnames(TEST),]

TRAIN.score = -log2(TRAIN.score$AUC_9d)
TEST.score = -log2(TEST.score$AUC_9d)

names(TRAIN.score) = colnames(TRAIN)
names(TEST.score) = colnames(TEST)

###########

train_size=dim(TRAIN)[2]
lho_size=train_size/2
count=dim(TRAIN)[1]

LHO_REP1  = sample(c(1:train_size),size=lho_size,replace = FALSE)
LHO_REP2  = sample(c(1:train_size),size=lho_size,replace = FALSE)
LHO_REP3  = sample(c(1:train_size),size=lho_size,replace = FALSE)
LHO_REP4  = sample(c(1:train_size),size=lho_size,replace = FALSE)
LHO_REP5  = sample(c(1:train_size),size=lho_size,replace = FALSE)
LHO_REP6  = sample(c(1:train_size),size=lho_size,replace = FALSE)
LHO_REP7  = sample(c(1:train_size),size=lho_size,replace = FALSE)
LHO_REP8  = sample(c(1:train_size),size=lho_size,replace = FALSE)
LHO_REP9  = sample(c(1:train_size),size=lho_size,replace = FALSE)
LHO_REP10 = sample(c(1:train_size),size=lho_size,replace = FALSE)

LHO = c()
LHO$REP1 = LHO_REP1
LHO$REP2 = LHO_REP2
LHO$REP3 = LHO_REP3
LHO$REP4 = LHO_REP4
LHO$REP5 = LHO_REP5
LHO$REP6 = LHO_REP6
LHO$REP7 = LHO_REP7
LHO$REP8 = LHO_REP8
LHO$REP9 = LHO_REP9
LHO$REP10 = LHO_REP10


TRAIN_R1 = TRAIN[,LHO_REP1]
TRAIN_S1 = TRAIN.score[colnames(TRAIN_R1)]
SpearCorr_R1<-apply(TRAIN_R1, 1, cor.test, TRAIN_S1, method="spearman")
rho_S1<-as.double(unlist(SpearCorr_R1)[grep("estimate.rho",names(unlist(SpearCorr_R1)))])
names(rho_S1)<-gsub(pattern=".estimate.rho",x=rownames(TRAIN_R1), replacement="")
pvalue_S1<-as.double(unlist(SpearCorr_R1)[grep("p.value",names(unlist(SpearCorr_R1)))])
names(pvalue_S1)<-gsub(pattern=".p.value",x=rownames(TRAIN_R1), replacement="")

TRAIN_R2 = TRAIN[,LHO_REP2]
TRAIN_S2 = TRAIN.score[colnames(TRAIN_R2)]
SpearCorr_R2<-apply(TRAIN_R2, 1, cor.test, TRAIN_S2, method="spearman")
rho_S2<-as.double(unlist(SpearCorr_R2)[grep("estimate.rho",names(unlist(SpearCorr_R2)))])
names(rho_S2)<-gsub(pattern=".estimate.rho",x=rownames(TRAIN_R2), replacement="")
pvalue_S2<-as.double(unlist(SpearCorr_R2)[grep("p.value",names(unlist(SpearCorr_R2)))])
names(pvalue_S2)<-gsub(pattern=".p.value",x=rownames(TRAIN_R2), replacement="")

TRAIN_R3 = TRAIN[,LHO_REP3]
TRAIN_S3 = TRAIN.score[colnames(TRAIN_R3)]
SpearCorr_R3<-apply(TRAIN_R3, 1, cor.test, TRAIN_S3, method="spearman")
rho_S3<-as.double(unlist(SpearCorr_R3)[grep("estimate.rho",names(unlist(SpearCorr_R3)))])
names(rho_S3)<-gsub(pattern=".estimate.rho",x=rownames(TRAIN_R3), replacement="")
pvalue_S3<-as.double(unlist(SpearCorr_R3)[grep("p.value",names(unlist(SpearCorr_R3)))])
names(pvalue_S3)<-gsub(pattern=".p.value",x=rownames(TRAIN_R3), replacement="")

TRAIN_R4 = TRAIN[,LHO_REP4]
TRAIN_S4 = TRAIN.score[colnames(TRAIN_R4)]
SpearCorr_R4<-apply(TRAIN_R4, 1, cor.test, TRAIN_S4, method="spearman")
rho_S4<-as.double(unlist(SpearCorr_R4)[grep("estimate.rho",names(unlist(SpearCorr_R4)))])
names(rho_S4)<-gsub(pattern=".estimate.rho",x=rownames(TRAIN_R4), replacement="")
pvalue_S4<-as.double(unlist(SpearCorr_R4)[grep("p.value",names(unlist(SpearCorr_R4)))])
names(pvalue_S4)<-gsub(pattern=".p.value",x=rownames(TRAIN_R4), replacement="")

TRAIN_R5 = TRAIN[,LHO_REP5]
TRAIN_S5 = TRAIN.score[colnames(TRAIN_R5)]
SpearCorr_R5<-apply(TRAIN_R5, 1, cor.test, TRAIN_S5, method="spearman")
rho_S5<-as.double(unlist(SpearCorr_R5)[grep("estimate.rho",names(unlist(SpearCorr_R5)))])
names(rho_S5)<-gsub(pattern=".estimate.rho",x=rownames(TRAIN_R5), replacement="")
pvalue_S5<-as.double(unlist(SpearCorr_R5)[grep("p.value",names(unlist(SpearCorr_R5)))])
names(pvalue_S5)<-gsub(pattern=".p.value",x=rownames(TRAIN_R5), replacement="")

TRAIN_R6 = TRAIN[,LHO_REP6]
TRAIN_S6 = TRAIN.score[colnames(TRAIN_R6)]
SpearCorr_R6<-apply(TRAIN_R6, 1, cor.test, TRAIN_S6, method="spearman")
rho_S6<-as.double(unlist(SpearCorr_R6)[grep("estimate.rho",names(unlist(SpearCorr_R6)))])
names(rho_S6)<-gsub(pattern=".estimate.rho",x=rownames(TRAIN_R6), replacement="")
pvalue_S6<-as.double(unlist(SpearCorr_R6)[grep("p.value",names(unlist(SpearCorr_R6)))])
names(pvalue_S6)<-gsub(pattern=".p.value",x=rownames(TRAIN_R6), replacement="")

TRAIN_R7 = TRAIN[,LHO_REP7]
TRAIN_S7 = TRAIN.score[colnames(TRAIN_R7)]
SpearCorr_R7<-apply(TRAIN_R7, 1, cor.test, TRAIN_S7, method="spearman")
rho_S7<-as.double(unlist(SpearCorr_R7)[grep("estimate.rho",names(unlist(SpearCorr_R7)))])
names(rho_S7)<-gsub(pattern=".estimate.rho",x=rownames(TRAIN_R7), replacement="")
pvalue_S7<-as.double(unlist(SpearCorr_R7)[grep("p.value",names(unlist(SpearCorr_R7)))])
names(pvalue_S7)<-gsub(pattern=".p.value",x=rownames(TRAIN_R7), replacement="")

TRAIN_R8 = TRAIN[,LHO_REP8]
TRAIN_S8 = TRAIN.score[colnames(TRAIN_R8)]
SpearCorr_R8<-apply(TRAIN_R8, 1, cor.test, TRAIN_S8, method="spearman")
rho_S8<-as.double(unlist(SpearCorr_R8)[grep("estimate.rho",names(unlist(SpearCorr_R8)))])
names(rho_S8)<-gsub(pattern=".estimate.rho",x=rownames(TRAIN_R8), replacement="")
pvalue_S8<-as.double(unlist(SpearCorr_R8)[grep("p.value",names(unlist(SpearCorr_R8)))])
names(pvalue_S8)<-gsub(pattern=".p.value",x=rownames(TRAIN_R8), replacement="")

TRAIN_R9 = TRAIN[,LHO_REP9]
TRAIN_S9 = TRAIN.score[colnames(TRAIN_R9)]
SpearCorr_R9<-apply(TRAIN_R9, 1, cor.test, TRAIN_S9, method="spearman")
rho_S9<-as.double(unlist(SpearCorr_R9)[grep("estimate.rho",names(unlist(SpearCorr_R9)))])
names(rho_S9)<-gsub(pattern=".estimate.rho",x=rownames(TRAIN_R9), replacement="")
pvalue_S9<-as.double(unlist(SpearCorr_R9)[grep("p.value",names(unlist(SpearCorr_R9)))])
names(pvalue_S9)<-gsub(pattern=".p.value",x=rownames(TRAIN_R9), replacement="")

TRAIN_R10 = TRAIN[,LHO_REP10]
TRAIN_S10 = TRAIN.score[colnames(TRAIN_R10)]
SpearCorr_R10<-apply(TRAIN_R10, 1, cor.test, TRAIN_S10, method="spearman")
rho_S10<-as.double(unlist(SpearCorr_R10)[grep("estimate.rho",names(unlist(SpearCorr_R10)))])
names(rho_S10)<-gsub(pattern=".estimate.rho",x=rownames(TRAIN_R10), replacement="")
pvalue_S10<-as.double(unlist(SpearCorr_R10)[grep("p.value",names(unlist(SpearCorr_R10)))])
names(pvalue_S10)<-gsub(pattern=".p.value",x=rownames(TRAIN_R10), replacement="")

LHO$RHO1 = rho_S1
LHO$RHO2 = rho_S2
LHO$RHO3 = rho_S3
LHO$RHO4 = rho_S4
LHO$RHO5 = rho_S5
LHO$RHO6 = rho_S6
LHO$RHO7 = rho_S7
LHO$RHO8 = rho_S8
LHO$RHO9 = rho_S9
LHO$RHO10 = rho_S10

LHO$pValue1 = pvalue_S1
LHO$pValue2 = pvalue_S2
LHO$pValue3 = pvalue_S3
LHO$pValue4 = pvalue_S4
LHO$pValue5 = pvalue_S5
LHO$pValue6 = pvalue_S6
LHO$pValue7 = pvalue_S7
LHO$pValue8 = pvalue_S8
LHO$pValue9 = pvalue_S9
LHO$pValue10 = pvalue_S10

rho_LHO<-cbind(rho_S1,rho_S2,rho_S3,rho_S4,rho_S5,rho_S6,rho_S7,rho_S8,rho_S9,rho_S10)
colnames(rho_LHO)<-c("REP1","REP2","REP3","REP4","REP5","REP6","REP7","REP8","REP9","REP10")

pvalue_LHO<-cbind(pvalue_S1,pvalue_S2,pvalue_S3,pvalue_S4,pvalue_S5,pvalue_S6,pvalue_S7,pvalue_S8,pvalue_S9,pvalue_S10)
colnames(pvalue_LHO)<-c("REP1","REP2","REP3","REP4","REP5","REP6","REP7","REP8","REP9","REP10")

meanrho_LHO<-rowMeans(rho_LHO)
meanpvalue_LHO<-rowMeans(pvalue_LHO)

mean_rhopval<-cbind(meanpvalue_LHO,meanrho_LHO)
colnames(mean_rhopval)<-c("pval","rho")

kept = mean_rhopval[which(mean_rhopval[,1]<0.05),]
keepids = rownames(kept)

sTRAIN = TRAIN[keepids,]
sTEST = TEST[keepids,]

opt.lambda05<-glmnet(t(sTRAIN),TRAIN.score, alpha=0.5,nlambda = 100)$lambda
grid05 <- expand.grid(.lambda = opt.lambda05, .alpha = rep(0.5,1))
trc = trainControl(method="repeatedcv", number=2, repeats = 100,allowParallel = TRUE)    

#seed <- 7
#metric <- "Accuracy"
#set.seed(seed)
#mtry <- ncol(sTRAIN)/3
#tunegrid <- expand.grid(.mtry=mtry)

cl <- makeCluster(6) #onserver
registerDoParallel(cl)
pX = train(t(sTRAIN), TRAIN.score, trControl = trc, method="glmnet", tuneGrid = grid05, preProc = c("center", "scale"))
#pY = train(t(sTRAIN), TRAIN.score, trControl = trc, method="rf", tuneGrid=tunegrid, preProc = c("center", "scale"))
stopCluster(cl)

pX.coefs = cbind(as.matrix(coef(pX$finalModel,s=pX$finalModel$lambdaOpt))[!(as.matrix(coef(pX$finalModel,s=pX$finalModel$lambdaOpt))==0),])
pX.coefs[order(pX.coefs,decreasing = T),]
pX.sig = rownames(pX.coefs)
pX.sig = pX.sig[2:length(pX.sig)]
pX.sig = pX.sig[grep("-AS1", pX.sig, invert = T)]

xTRAIN = TRAIN[pX.sig,]
xTEST = TEST[pX.sig,]

opt.lambda0<-glmnet(t(sTRAIN),TRAIN.score, alpha=0,nlambda = 100)$lambda
grid0 <- expand.grid(.lambda = opt.lambda0, .alpha = rep(0,1))
trc = trainControl(method="repeatedcv", number=2, repeats = 100,allowParallel = TRUE)    

cl <- makeCluster(6) #onserver
registerDoParallel(cl)
pXX = train(t(xTRAIN), TRAIN.score, trControl = trc, method="glmnet", tuneGrid = grid0, preProc = c("center", "scale"))
stopCluster(cl)

pXX.coefs = cbind(as.matrix(coef(pXX$finalModel,s=pXX$finalModel$lambdaOpt))[!(as.matrix(coef(pXX$finalModel,s=pXX$finalModel$lambdaOpt))==0),])
pXX.coefs[order(pXX.coefs,decreasing = T),]
cbind(pXX.coefs[order(pXX.coefs,decreasing = T),][2:length(pXX.coefs[order(pXX.coefs,decreasing = T),])])

predX = predict(pX, t(sTEST))
plot(predX, TEST.score)
cor.test(predX, TEST.score)

predXX = predict(pXX, t(xTEST))
plot(predXX, TEST.score)
cor.test(predXX, TEST.score)
plot(predXX, TEST.score)

