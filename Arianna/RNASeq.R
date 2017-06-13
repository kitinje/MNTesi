#### FASE 1 - Importazione dei dati

setwd("/Users/age/Desktop/ARIANNA/ATRA/")
dir()

targets = dir(pattern = "*.tab")

temp = read.table(targets[1], sep="\t", skip = 4, row.names = 1, stringsAsFactors = F) #skip 4 perchè nelle prime 4 righe ci sono commenti
head(temp)

ngenes = nrow(temp)
nsamples = length(targets)

RAWCOUNTS = matrix(nrow = ngenes, ncol = nsamples, data = 0)
head(RAWCOUNTS)
colnames(RAWCOUNTS) = targets
rownames(RAWCOUNTS) = rownames(temp)
head(RAWCOUNTS)

for(i in 1:length(targets))
{
  print(i)
  RAWCOUNTS[,i] = read.table(targets[i],sep="\t",row.names = 1, skip = 4, stringsAsFactors = FALSE, header = FALSE)[,3]
}


pb <- txtProgressBar(min = 0, max = length(targets), style = 3)
for(i in 1:length(targets))
{
  setTxtProgressBar(pb, i, title = "Working")
  RAWCOUNTS[,i] = read.table(targets[i],sep="\t",row.names = 1, skip = 4, stringsAsFactors = FALSE, header = FALSE)[,3]
}

head(RAWCOUNTS)

colnames(RAWCOUNTS) = gsub( colnames(RAWCOUNTS), pattern = "_ReadsPerGene.out.tab", replacement = "" )

head(RAWCOUNTS)
dim(RAWCOUNTS)
typeof(RAWCOUNTS)

# FASE 2 - Visualizzazione e normalizzazione 
library(edgeR)
y.raw <- DGEList(counts=RAWCOUNTS)

names(y.raw)

head(y.raw$counts)
y.raw$samples

# remove everything after "." character
keep=grep(cbind(unlist(strsplit(rownames(y.raw),'.',fixed=TRUE))), pattern = "ENSG")
rownames(y.raw)<-cbind(unlist(strsplit(rownames(y.raw),'.',fixed=TRUE)))[keep]

# keep files genes showing more than 1 cpm in at least 3 samples. 
exp_thr=1 #expression threshold
samp_thr=3 #sample threshold
keep <- rowSums(cpm(y.raw)>exp_thr) >= samp_thr
y <- y.raw[keep, , keep.lib.sizes=FALSE] #ricalcola anche library sizes (consigliato)

# Normalization and TMM correction
y <- calcNormFactors(y)

y.raw$samples
y$samples

y$samples$group = factor( c( rep("A",3), rep("C", 3) )  )

y$samples

#library sizes
barplot(y$samples$lib.size*1e-6, names=as.character(rownames(y$samples)), ylab="Library size (millions)", xlab = "included samples") 

library(plotly)
librarysize_df <- data.frame("mapped_counts" = y$samples$lib.size*1e-6, "sample" = colnames(y$counts), "cline" = y$samples$group)
library_plot <- plot_ly(
  data=librarysize_df,
  x = seq(1:dim(y$counts)[2]),
  y = ~mapped_counts, 
  type = 'bar',
  color = ~cline,
  hoverinfo = "text",
  text = ~paste(sample, ' / ', mapped_counts, ' millions reads')
)  %>%
  layout(legend = list(x = 0.84, y = 1))
library_plot


# PCA
y.cpm.log2 = cpm(y, log = TRUE, prior.count = 1)
hist(y.cpm.log2, breaks = 20)
plot(density(y.cpm.log2))


pcaY <- prcomp(t(y.cpm.log2), scale. = TRUE)
plot(pcaY, type = "l") # Prinicipal components explain variance:

pcaDf <- data.frame(pcaY$x, "group" = y$samples$group)


library(ggplot2)
library(ggrepel)

# 2D STATIC PCA-PLOT
plot.pca_2d <- ggplot(pcaDf,aes(x=PC1, y=PC2)) +
  theme_classic() +
  geom_hline(yintercept = 0, color = "gray70") +
  geom_vline(xintercept = 0, color = "gray70") +
  geom_point(aes(color = group), alpha = 0.55, size = 3) +
  xlab("PC1") +
  ylab("PC2") + 
  ggtitle("PCA from cell-line data")
plot.pca_2d_numbered = plot.pca_2d + geom_text_repel(aes(y = PC2 + 0.25, label = c(seq(1:dim(y)[2]) )))
plot.pca_2d_numbered

library(plotly)
# 2D INTERACTIVE PCA-PLOT
plot.pca_2d_int <- plot_ly(pcaDf, x = pcaDf[,1] , y = pcaDf[,2],  text = pcaDf$group, mode = "markers", color = pcaDf$group,  marker = list(size = 11), type = "scatter" ) 
plot.pca_2d_int <- layout(plot.pca_2d_int,title = "Principal Component Analysis", 
                          xaxis = list(title = "PC 1"),
                          yaxis = list(title = "PC 2"))
plot.pca_2d_int


# 3D INTERACTIVE PCA-PLOT
plot.pca_3d_int <- plot_ly(pcaDf, x = pcaDf[,1] , y = pcaDf[,2], z = pcaDf[,3],  text = pcaDf$group,
                           mode = "markers", color = pcaDf$group, marker = list(size = 7), type="scatter3d") 
plot.pca_3d_int <- layout(plot.pca_3d_int,  title = "Principal Component Analysis")
plot.pca_3d_int


#hierarchical clustering
library(ape)
data <- scale(t(y.cpm.log2))
hc <- hclust( dist(data, method = "euclidean"), method = "ward.D2" )
hcd <- as.dendrogram(hc)
colTips= c( rep("red",3), rep("blue",3) )

plot(as.phylo(hc), type = "phylogram", cex = 0.8, edge.color = "black", edge.width = 2, edge.lty = 1, tip.color = colTips)
plot(as.phylo(hc), type = "cladogram", cex = 0.8, edge.color = "black", edge.width = 2, edge.lty = 1, tip.color = colTips)
plot(as.phylo(hc), type = "unrooted", cex = 0.7, edge.color = "black", edge.width = 1, edge.lty = 1, tip.color = colTips, lab4ut='axial', no.margin=TRUE)
plot(as.phylo(hc), type = "fan", cex = 0.8, edge.color = "black", edge.width = 2, edge.lty = 1, tip.color = colTips)

A <- factor(c(1,1,1,0,0,0))
design = model.matrix(~0 + A)
colnames(design)=c("C","A")

y.de <- estimateDisp(y, design, robust=TRUE)
plotBCV(y.de) #plot NB dispersion
fit_int <- glmQLFit(y.de, design, robust=TRUE)
AvsC=makeContrasts(A-C, levels=design)
y.de.AvsC <- glmQLFTest(fit_int,  contrast = AvsC)


topTags(y.de.AvsC) # primi maggiormente DE

y.de.AvsC.res <- topTags(y.de.AvsC, n=Inf, sort.by = "none")$table
y.de.AvsC.res <- topTags(y.de.AvsC, n=Inf, sort.by = "p.value")$table

library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- rownames(y.cpm.log2)
symbols <- genes
symbols[] = NA
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "entrezgene", "hgnc_symbol", "description"),values=genes,mart= mart)
G_list = G_list[!(duplicated(G_list$ensembl_gene_id)), ]
rownames(G_list)=G_list$ensembl_gene_id
G_list=G_list[rownames(y.cpm.log2),]
annotations = G_list[,2:3]

y.de.AvsC.res$Symbol=annotations[rownames(y.de.AvsC.res),2]
y.de.AvsC.res$Entrez=annotations[rownames(y.de.AvsC.res),1]

head(y.de.AvsC.res, n=100)

results.de = y.de.AvsC.res[y.de.AvsC.res$FDR<0.05,]

write.table(results.de, file="results.txt", sep="\t", col.names=NA)

# PATHWAY

#load PATHWAYS
# pathway da MSIGDB da broad institute

#PATHWAYS
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_H_v5p1.rdata"))
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c1_v5p1.rdata"))
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c2_v5p1.rdata"))
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c3_v5p1.rdata"))
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c4_v5p1.rdata"))
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c5_v5p1.rdata"))
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c6_v5p1.rdata"))
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c7_v5p1.rdata"))




# Hs.H
typeof(Hs.H)
names(Hs.H)
Hs.H$HALLMARK_NOTCH_SIGNALING

library(org.Hs.eg.db)
# Gene-set analysis
xx <- as.list(org.Hs.egENSEMBL2EG)
ann <- cbind(unlist(xx)); colnames(ann) = "ENTREZID"
keep = cbind(ann[(rownames(ann) %in% rownames(y)),]); colnames(keep)="ENTREZID"
keep = cbind(keep[which(!duplicated(keep[,1])),]) ; colnames(keep)="ENTREZID"
y.Hs=y[rownames(keep),]
rownames(y.Hs)=keep[,1]

y.Hs <- estimateDisp(y.Hs, design, robust=TRUE)
y.idx.H <- ids2indices(Hs.H,id=rownames(y.Hs)) #molto importante senno sfasa
y.cam.H.AvsC <- camera(y.Hs, y.idx.H, design,inter.gene.cor=0.01, contrast=AvsC, sort = T)
colnames(y.cam.H.AvsC)=c("GENESET_SIZE", "AvsC_dir", "AvsC_pvalue",  "AvsC_fdr")
y.cam.H.AvsC


#selezioni un pathway
# HALLMARK_MYC_TARGETS_V2 

selgeni = Hs.H$HALLMARK_MYC_TARGETS_V2 

y.Hs.mycv2  = y.Hs$counts[which(rownames(y.Hs$counts) %in% selgeni), ]
y.Hs.mycv2.log2 = cpm(y.Hs.mycv2, log = TRUE, prior.count = 1)

library(pheatmap)
pheatmap(y.Hs.mycv2.log2,scale="row",  clustering_distance_rows = "correlation", clustering_distance_cols = "correlation")


controlli = y.Hs.mycv2.log2[,4:6]
avg.controlli = cbind(rowMeans(controlli))

a1 = y.Hs.mycv2.log2[,1]-avg.controlli
a2 = y.Hs.mycv2.log2[,2]-avg.controlli
a3 = y.Hs.mycv2.log2[,3]-avg.controlli

atraFC = cbind(a1,a2,a3)
colnames(atraFC) = c("A1", "A2", "A3")
pheatmap(atraFC, scale="none",  clustering_distance_rows = "euclidean", clustering_distance_cols = "correlation")

annot1 = annotations[which(annotations$entrezgene %in% rownames(atraFC)),] # non sono necessariamente in ordine
rownames(annot1)= annot1$entrezgene
annot1 = annot1[rownames(atraFC),]
rownames(atraFC) = annot1$hgnc_symbol
pheatmap(atraFC, scale="none",  clustering_distance_rows = "euclidean", clustering_distance_cols = "correlation")

