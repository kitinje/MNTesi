bk = unique(c(seq(-10,-5, length=2), seq(-5, 5, by=0.1), seq(5,10, length=2)))
col1 = rep("#4575B4", 1)
col2 = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
col3 <- rep("#D73027", 1)
colors3 <- c(col1, col2,col3)


hKDM = t(KDM[,1:7])

pheatmap(hKDM,scale="none", breaks = bk, clustering_distance_rows = "correlation", clustering_distance_cols = "euclidean", color = colors3)


pheatmap(hKDM,scale="none", clustering_distance_rows = "euclidean", clustering_distance_cols = "correlation")

pheatmap(hKDM,scale="row",clustering_distance_cols = "correlation")


ggplot(KDM,  aes(x=color, y=as.numeric(as.character(KDM4A)), fill=color, alpha=0.1))  + guides(fill=FALSE) +  geom_jitter( aes(color=color),size=2.5,position = position_jitter(width = .15)) + geom_boxplot(outlier.shape=NA) + theme_bw() + theme(legend.position="none")

ggplot(KDM,  aes(x=color, y=as.numeric(as.character(KDM4B)), fill=color, alpha=0.1))  + guides(fill=FALSE) +  geom_jitter( aes(color=color),size=2.5,position = position_jitter(width = .15)) + geom_boxplot(outlier.shape=NA) + theme_bw() + theme(legend.position="none")


KDM.exp =  read.table("~/Desktop/DANIMARCA/HMAP/kdm.txt",sep="\t", header = TRUE, stringsAsFactors = F)
KDM.exp <- as.data.frame(KDM.exp)  

ggplot(KDM.exp,  aes(x=GROUP, y=as.numeric(as.character(EXP)), fill=COLOR, alpha=0.1))  + guides(fill=FALSE) +  geom_jitter( aes(color=COLOR),size=1.5,position = position_jitter(width = .15)) + geom_boxplot(outlier.shape=NA) + theme_bw() + theme(legend.position="none") +  theme(axis.text.x = element_text(angle = 90, hjust = 1))


  

library(org.Hs.eg.db)
dim(y.tcga)
xx <- as.list(org.Hs.egENSEMBL2EG)
ann <- cbind(unlist(xx)); colnames(ann) = "ENTREZID"
keep = cbind(ann[(rownames(ann) %in% rownames(y.tcga)),]); colnames(keep)="ENTREZID"
keep = cbind(keep[which(!duplicated(keep[,1])),]) ; colnames(keep)="ENTREZID"
y.tcga.Hs=y.tcga[rownames(keep),]
rownames(y.tcga.Hs)=keep[,1]
y.tcga.Hs$genes$GeneID=rownames(y.tcga.Hs)
entrez_size=dim(y.tcga.Hs)[1]


y.tcga.Hs$genes$Symbol <- mapIds(org.Hs.eg.db, rownames(y.tcga.Hs), keytype="ENTREZID", column="SYMBOL")
rownames(y.tcga.Hs$counts) = as.character(y.tcga.Hs$genes$Symbol)
rownames(y.tcga.Hs) = y.tcga.Hs$genes$Symbol
y.tcga.log2 = cpm(y.tcga.Hs, log = TRUE, prior.count = 1)


KDM.tcga = y.tcga.log2[c("KDM4A", "KDM4B", "KDM4C", "KDM5A", "KDM5B", "KDM5C", "KDM5D"),]

KDM.tcga = t(KDM.tcga)

KDM.tcga_info = read.table("~/Desktop/DANIMARCA/HMAP/ihc_tcga.txt",sep="\t", header = T, row.names=1, stringsAsFactors = F)

KDM.tcga = KDM.tcga[rownames(KDM.tcga_info),]
KDM.tcga = as.data.frame(KDM.tcga)

KDM.tcga$GROUP = KDM.tcga_info$GROUP

hKDM.tcga = matrix(nrow = nrow(KDM.tcga)*7, ncol = 2)

espressione = c()
gruppo = c()
for(i in 1:7) {
  espressione = c(espressione, KDM.tcga[,i])
  temp = KDM.tcga$GROUP
  currGene = colnames(KDM.tcga)[i]
  gruppo = c(gruppo, paste(temp, currGene ,sep = "_") )
}



hKDM.tcga[,1] = espressione;
hKDM.tcga[,2] = gruppo;

colnames(hKDM.tcga) = c("EXP", "GROUP")
hKDM.tcga = as.data.frame(hKDM.tcga)
hKDM.tcga[,1] = as.numeric(as.character(hKDM.tcga[,1]))

CLR = c( rep("KDM4A",923), rep("KDM4B",923), rep("KDM4C",923), rep("KDM5A",923), rep("KDM5B",923), rep("KDM5C",923), rep("KDM5D",923))
hKDM.tcga$COLOR = CLR

ggplot(hKDM.tcga,  aes(x=GROUP, y=as.numeric(as.character(EXP)), fill=as.factor(COLOR), alpha=0.1))  + guides(fill=FALSE) +  geom_jitter( aes(color=COLOR),size=0.05,position = position_jitter(width = .15)) + geom_boxplot(outlier.shape=NA) + theme_bw() + theme(legend.position="none") +  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(KDM.exp,  aes(x=GROUP, y=as.numeric(as.character(EXP)), fill=as.factor(COLOR), alpha=0.1))  + guides(fill=FALSE) +  geom_jitter( aes(color=COLOR),size=0.05,position = position_jitter(width = .15)) + geom_boxplot(outlier.shape=NA) + theme_bw() + theme(legend.position="none") +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
