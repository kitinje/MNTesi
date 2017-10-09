source("~/Desktop/multiplot.R")

library("edgeR")

setwd("~/Desktop/EG/NR")

AML = read.table("nrec_aml.txt", sep="\t", stringsAsFactors = T, row.names = 1, header = T)

annotations = read.table("nrec_annotations.txt", sep="\t", stringsAsFactors = T, row.names = 1, header = T)

AML.orphan = AML[,c("FAB","HIST",rownames(annotations)[annotations$TYPE=="Orphan"])]

typeof(AML.orphan)
AML.orphan = as.data.frame(AML.orphan)
typeof(AML.orphan)

library(ggplot2)

bplots = list()
for(i in 3:17) {
  print(i)
  p1 = ggplot(AML.orphan,  aes(x=FAB, y=AML.orphan[,i], fill=FAB, alpha=0.1))  + guides(fill=FALSE) +  geom_jitter( aes(color=FAB) ,size=1.5, position = position_jitter(width = .10)) + geom_boxplot(outlier.shape=NA) + theme_bw() + theme(legend.position="none")
  bplots[[i-2]] <- p1
}

p1 = ggplot(AML.orphan,  aes(x=FAB, y=AML.orphan[,3], fill=FAB, alpha=0.1))  + guides(fill=FALSE) +  geom_jitter( aes(color=FAB) ,size=1.5, position = position_jitter(width = .10)) + geom_boxplot(outlier.shape=NA) + theme_bw() + theme(legend.position="none") + ylab(label = colnames(AML.orphan)[3]) 
p2 = ggplot(AML.orphan,  aes(x=FAB, y=AML.orphan[,4], fill=FAB, alpha=0.1))  + guides(fill=FALSE) +  geom_jitter( aes(color=FAB) ,size=1.5, position = position_jitter(width = .10)) + geom_boxplot(outlier.shape=NA) + theme_bw() + theme(legend.position="none") + ylab(label = colnames(AML.orphan)[4]) 
p3 = ggplot(AML.orphan,  aes(x=FAB, y=AML.orphan[,5], fill=FAB, alpha=0.1))  + guides(fill=FALSE) +  geom_jitter( aes(color=FAB) ,size=1.5, position = position_jitter(width = .10)) + geom_boxplot(outlier.shape=NA) + theme_bw() + theme(legend.position="none") + ylab(label = colnames(AML.orphan)[5]) 
p4 = ggplot(AML.orphan,  aes(x=FAB, y=AML.orphan[,6], fill=FAB, alpha=0.1))  + guides(fill=FALSE) +  geom_jitter( aes(color=FAB) ,size=1.5, position = position_jitter(width = .10)) + geom_boxplot(outlier.shape=NA) + theme_bw() + theme(legend.position="none") + ylab(label = colnames(AML.orphan)[6]) 
p5 = ggplot(AML.orphan,  aes(x=FAB, y=AML.orphan[,7], fill=FAB, alpha=0.1))  + guides(fill=FALSE) +  geom_jitter( aes(color=FAB) ,size=1.5, position = position_jitter(width = .10)) + geom_boxplot(outlier.shape=NA) + theme_bw() + theme(legend.position="none") + ylab(label = colnames(AML.orphan)[7]) 
p6 = ggplot(AML.orphan,  aes(x=FAB, y=AML.orphan[,8], fill=FAB, alpha=0.1))  + guides(fill=FALSE) +  geom_jitter( aes(color=FAB) ,size=1.5, position = position_jitter(width = .10)) + geom_boxplot(outlier.shape=NA) + theme_bw() + theme(legend.position="none") + ylab(label = colnames(AML.orphan)[8]) 
p7 = ggplot(AML.orphan,  aes(x=FAB, y=AML.orphan[,9], fill=FAB, alpha=0.1))  + guides(fill=FALSE) +  geom_jitter( aes(color=FAB) ,size=1.5, position = position_jitter(width = .10)) + geom_boxplot(outlier.shape=NA) + theme_bw() + theme(legend.position="none") + ylab(label = colnames(AML.orphan)[9]) 
p8 = ggplot(AML.orphan,  aes(x=FAB, y=AML.orphan[,10], fill=FAB, alpha=0.1))  + guides(fill=FALSE) +  geom_jitter( aes(color=FAB) ,size=1.5, position = position_jitter(width = .10)) + geom_boxplot(outlier.shape=NA) + theme_bw() + theme(legend.position="none") + ylab(label = colnames(AML.orphan)[10]) 
p9 = ggplot(AML.orphan,  aes(x=FAB, y=AML.orphan[,11], fill=FAB, alpha=0.1))  + guides(fill=FALSE) +  geom_jitter( aes(color=FAB) ,size=1.5, position = position_jitter(width = .10)) + geom_boxplot(outlier.shape=NA) + theme_bw() + theme(legend.position="none") + ylab(label = colnames(AML.orphan)[11]) 
p10 = ggplot(AML.orphan,  aes(x=FAB, y=AML.orphan[,12], fill=FAB, alpha=0.1))  + guides(fill=FALSE) +  geom_jitter( aes(color=FAB) ,size=1.5, position = position_jitter(width = .10)) + geom_boxplot(outlier.shape=NA) + theme_bw() + theme(legend.position="none") + ylab(label = colnames(AML.orphan)[12]) 
p11 = ggplot(AML.orphan,  aes(x=FAB, y=AML.orphan[,13], fill=FAB, alpha=0.1))  + guides(fill=FALSE) +  geom_jitter( aes(color=FAB) ,size=1.5, position = position_jitter(width = .10)) + geom_boxplot(outlier.shape=NA) + theme_bw() + theme(legend.position="none") + ylab(label = colnames(AML.orphan)[13]) 
p12 = ggplot(AML.orphan,  aes(x=FAB, y=AML.orphan[,14], fill=FAB, alpha=0.1))  + guides(fill=FALSE) +  geom_jitter( aes(color=FAB) ,size=1.5, position = position_jitter(width = .10)) + geom_boxplot(outlier.shape=NA) + theme_bw() + theme(legend.position="none") + ylab(label = colnames(AML.orphan)[14]) 
p13 = ggplot(AML.orphan,  aes(x=FAB, y=AML.orphan[,15], fill=FAB, alpha=0.1))  + guides(fill=FALSE) +  geom_jitter( aes(color=FAB) ,size=1.5, position = position_jitter(width = .10)) + geom_boxplot(outlier.shape=NA) + theme_bw() + theme(legend.position="none") + ylab(label = colnames(AML.orphan)[15]) 
p14 = ggplot(AML.orphan,  aes(x=FAB, y=AML.orphan[,16], fill=FAB, alpha=0.1))  + guides(fill=FALSE) +  geom_jitter( aes(color=FAB) ,size=1.5, position = position_jitter(width = .10)) + geom_boxplot(outlier.shape=NA) + theme_bw() + theme(legend.position="none") + ylab(label = colnames(AML.orphan)[16]) 
p15 = ggplot(AML.orphan,  aes(x=FAB, y=AML.orphan[,17], fill=FAB, alpha=0.1))  + guides(fill=FALSE) +  geom_jitter( aes(color=FAB) ,size=1.5, position = position_jitter(width = .10)) + geom_boxplot(outlier.shape=NA) + theme_bw() + theme(legend.position="none") + ylab(label = colnames(AML.orphan)[17]) 


multiplot(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15, cols = 3)


