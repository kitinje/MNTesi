# HCC1187 TAB files in /pico/scratch/userexternal/mfratell/MARTINA/RNASeq/
# HS578T TAB files in /pico/scratch/userexternal/mfratell/ENRICO/RNASeq/TAB
# HCC202 TAB files in /pico/scratch/userexternal/mfratell/HCC202/TAB
# le altre /pico/scratch/userexternal/mfratell/JOB/TAB
# scp mfratell@login.pico.cineca.it:/pico/scratch/userexternal/mfratell/MARTINA/RNASeq_R/TAB/ ./

library(edgeR)
library(plotly)
library(ggplot2)
library(ggrepel)
library(ape)
library(biomaRt)
library(org.Hs.eg.db)

setwd("/Users/atra/Desktop/NUOVE/")

targets = dir(pattern = "*.tab")

temp = read.table(targets[1], sep="\t", skip = 4, row.names = 1, stringsAsFactors = F)
ngenes = nrow(temp)
nsamples = length(targets)

RAWCOUNTS = matrix(nrow = ngenes, ncol = nsamples, data = 0)
colnames(RAWCOUNTS) = targets
rownames(RAWCOUNTS) = rownames(temp)

for(i in 1:length(targets))
{
  RAWCOUNTS[,i] = read.table(targets[i],sep="\t",row.names = 1, skip = 4, stringsAsFactors = FALSE, header = FALSE)[,3]
}

colnames(RAWCOUNTS) = gsub( colnames(RAWCOUNTS), pattern = "_001_ReadsPerGene.out.tab", replacement = "" )

head(RAWCOUNTS)
dim(RAWCOUNTS)
typeof(RAWCOUNTS)

y.raw <- DGEList(counts=RAWCOUNTS)
y = y.raw

# discard EnsemblID version number and remove not valid IDS such as ENSG00000001.1_PAR_Y
y <-y[rownames(y)[grep("_",rownames(y), invert = TRUE)],]

# remove everything after "." character
keep=grep(cbind(unlist(strsplit(rownames(y),'.',fixed=TRUE))), pattern = "ENSG")
rownames(y)<-cbind(unlist(strsplit(rownames(y),'.',fixed=TRUE)))[keep]

# keep files genes showing more than 1 cpm in at least 3 samples. 
exp_thr=1
samp_thr=3
keep <- rowSums(cpm(y)>exp_thr) >= samp_thr
y<- y[keep, , keep.lib.sizes=FALSE] #ricalcola anche library sizes (consigliato)

# Normalization and TMM correction
y <- calcNormFactors(y)

y$samples$group = factor( rep(c(rep("A",3), rep("C", 3)),3) )

#library sizes
barplot(y$samples$lib.size*1e-6, names=as.character(rownames(y$samples)), ylab="Library size (millions)", xlab = "included samples") 

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
hist(y.cpm.log2, breaks = 100)
plot(density(y.cpm.log2))

y$samples$cell_line <- factor( c(rep("HCC1187",6), rep("HCC202", 6), rep("HS578T", 6)))
y$samples$treatment = factor( rep(c(rep("A",3), rep("C", 3)),3) )

pcaY <- prcomp(t(y.cpm.log2), scale. = TRUE)

pcaDf <- data.frame(pcaY$x, "cell_line" = y$samples$cell_line, "treatment" = y$samples$treatment, "group" = y$samples$group)
plot(pcaY, type = "l") # Prinicipal components explain variance:

# 2D STATIC PCA-PLOT
plot.pca_2d <- ggplot(pcaDf,aes(x=PC1, y=PC2)) +
  theme_classic() +
  geom_hline(yintercept = 0, color = "gray70") +
  geom_vline(xintercept = 0, color = "gray70") +
  geom_point(aes(color = cell_line, shape=treatment), alpha = 0.55, size = 3) +
  xlab("PC1") +
  ylab("PC2") + 
  ggtitle("PCA from cell-line data")
plot.pca_2d_numbered = plot.pca_2d + geom_text_repel(aes(y = PC2 + 0.25, label = c(seq(1:dim(y)[2]) )))
plot.pca_2d_numbered

# 2D INTERACTIVE PCA-PLOT
plot.pca_2d_int <- plot_ly(pcaDf, x = pcaDf[,1] , y = pcaDf[,2],  text = pcaDf$group, mode = "markers", color = pcaDf$cell_line, symbol = y$samples$treatment,  marker = list(size = 11), type = "scatter" ) 
plot.pca_2d_int <- layout(plot.pca_2d_int,title = "Principal Component Analysis", 
                          xaxis = list(title = "PC 1"),
                          yaxis = list(title = "PC 2"))
plot.pca_2d_int

# 3D INTERACTIVE PCA-PLOT
plot.pca_3d_int <- plot_ly(pcaDf, x = pcaDf[,1] , y = pcaDf[,2], z = pcaDf[,3],  text = pcaDf$group,
                           mode = "markers", color = pcaDf$cell_line, marker = list(size = 7), type="scatter3d", symbol = y$samples$treatment) 
plot.pca_3d_int <- layout(plot.pca_3d_int,  title = "Principal Component Analysis")
plot.pca_3d_int

#hierarchical clustering
data <- scale(t(y.cpm.log2))
hc <- hclust( dist(data, method = "euclidean"), method = "ward.D2" )
hcd <- as.dendrogram(hc)
colTips= c( rep("forestgreen",6), rep("magenta",12), rep("dodgerblue4",6), rep("red",6), rep("darkmagenta",6), rep("orange4",6))

plot(as.phylo(hc), type = "phylogram", cex = 0.8, edge.color = "black", edge.width = 2, edge.lty = 1, tip.color = colTips)
plot(as.phylo(hc), type = "cladogram", cex = 0.8, edge.color = "black", edge.width = 2, edge.lty = 1, tip.color = colTips)
plot(as.phylo(hc), type = "unrooted", cex = 0.7, edge.color = "black", edge.width = 1, edge.lty = 1, tip.color = colTips, lab4ut='axial', no.margin=TRUE)
plot(as.phylo(hc), type = "fan", cex = 0.8, edge.color = "black", edge.width = 2, edge.lty = 1, tip.color = colTips)


# CELL-LINE HCC1187
hcc1187<- y[,grep("^HCC1187", colnames(y))]
dim(hcc1187)

#create design matrix

A <- factor(c(1,1,1,0,0,0))
design = model.matrix(~0 + A)
colnames(design)=c("C","A")

hcc1187 <- estimateDisp(hcc1187, design, robust=TRUE)
plotBCV(hcc1187) #plot NB dispersion
fit_int <- glmQLFit(hcc1187, design, robust=TRUE)

AvsC=makeContrasts(A-C, levels=design)

hcc1187.qlf.AvsC <- glmQLFTest(fit_int,  contrast = AvsC)
topTags(hcc1187.qlf.AvsC)
hcc1187.AvsC <- topTags(hcc1187.qlf.AvsC, n=Inf, sort.by = "p.value")$table

dim(hcc1187.AvsC)

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- rownames(y.cpm.log2)
symbols <- genes
symbols[] = NA
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "entrezgene", "hgnc_symbol", "description"),values=genes,mart= mart)
G_list = G_list[!(duplicated(G_list$ensembl_gene_id)), ]
rownames(G_list)=G_list$ensembl_gene_id
G_list=G_list[rownames(y.cpm.log2),]
annotations = G_list[,2:4]

hcc1187.AvsC$Symbol=annotations[rownames(hcc1187.AvsC),2]
hcc1187.AvsC$Entrez=annotations[rownames(hcc1187.AvsC),1]

hcc1187.AvsC.fdr_05 = hcc1187.AvsC[which(hcc1187.AvsC$FDR<0.05),]
hcc1187.AvsC.fdr_05_up = hcc1187.AvsC.fdr_05[which(hcc1187.AvsC.fdr_05$logFC>0),]
hcc1187.AvsC.fdr_05_dn = hcc1187.AvsC.fdr_05[which(hcc1187.AvsC.fdr_05$logFC<0),]
hcc1187.AvsC.fdr_05[which(hcc1187.AvsC.fdr_05$Symbol=="MYC"),]

# CELL-LINE HCC202
hcc202<- y[,grep("^HCC202", colnames(y))]
dim(hcc202)

#create design matrix

A <- factor(c(1,1,1,0,0,0))
design = model.matrix(~0 + A)
colnames(design)=c("C","A")

hcc202 <- estimateDisp(hcc202, design, robust=TRUE)
plotBCV(hcc202) #plot NB dispersion
fit_int <- glmQLFit(hcc202, design, robust=TRUE)

AvsC=makeContrasts(A-C, levels=design)

hcc202.qlf.AvsC <- glmQLFTest(fit_int,  contrast = AvsC)
topTags(hcc202.qlf.AvsC)
hcc202.AvsC <- topTags(hcc202.qlf.AvsC, n=Inf, sort.by = "p.value")$table

dim(hcc202.AvsC)

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- rownames(y.cpm.log2)
symbols <- genes
symbols[] = NA
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "entrezgene", "hgnc_symbol", "description"),values=genes,mart= mart)
G_list = G_list[!(duplicated(G_list$ensembl_gene_id)), ]
rownames(G_list)=G_list$ensembl_gene_id
G_list=G_list[rownames(y.cpm.log2),]
annotations = G_list[,2:4]

hcc202.AvsC$Symbol=annotations[rownames(hcc202.AvsC),2]
hcc202.AvsC$Entrez=annotations[rownames(hcc202.AvsC),1]

hcc202.AvsC.fdr_05 = hcc202.AvsC[which(hcc202.AvsC$FDR<0.05),]
hcc202.AvsC.fdr_05_up = hcc202.AvsC.fdr_05[which(hcc202.AvsC.fdr_05$logFC>0),]
hcc202.AvsC.fdr_05_dn = hcc202.AvsC.fdr_05[which(hcc202.AvsC.fdr_05$logFC<0),]
hcc202.AvsC.fdr_05[which(hcc202.AvsC.fdr_05$Symbol=="MYC"),]

# CELL-LINE HS578T
hs578t<- y[,grep("^HS578T", colnames(y))]
dim(hs578t)

#create design matrix

A <- factor(c(1,1,1,0,0,0))
design = model.matrix(~0 + A)
colnames(design)=c("C","A")

hs578t <- estimateDisp(hs578t, design, robust=TRUE)
plotBCV(hs578t) #plot NB dispersion
fit_int <- glmQLFit(hs578t, design, robust=TRUE)

AvsC=makeContrasts(A-C, levels=design)

hs578t.qlf.AvsC <- glmQLFTest(fit_int,  contrast = AvsC)
topTags(hs578t.qlf.AvsC)
hs578t.AvsC <- topTags(hs578t.qlf.AvsC, n=Inf, sort.by = "p.value")$table

dim(hs578t.AvsC)

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- rownames(y.cpm.log2)
symbols <- genes
symbols[] = NA
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "entrezgene", "hgnc_symbol", "description"),values=genes,mart= mart)
G_list = G_list[!(duplicated(G_list$ensembl_gene_id)), ]
rownames(G_list)=G_list$ensembl_gene_id
G_list=G_list[rownames(y.cpm.log2),]
annotations = G_list[,2:4]

hs578t.AvsC$Symbol=annotations[rownames(hs578t.AvsC),2]
hs578t.AvsC$Entrez=annotations[rownames(hs578t.AvsC),1]

hs578t.AvsC.fdr_05 = hs578t.AvsC[which(hs578t.AvsC$FDR<0.05),]

hs578t.AvsC.fdr_05[which(hs578t.AvsC.fdr_05$Symbol=="RARB"),]

# Gene-set analysis
xx <- as.list(org.Hs.egENSEMBL2EG)
ann <- cbind(unlist(xx)); colnames(ann) = "ENTREZID"
keep = cbind(ann[(rownames(ann) %in% rownames(y)),]); colnames(keep)="ENTREZID"
keep = cbind(keep[which(!duplicated(keep[,1])),]) ; colnames(keep)="ENTREZID"

y.Hs=y[rownames(keep),]
rownames(y.Hs)=keep[,1]
y.Hs$genes$GeneID=rownames(y.Hs) #boh
entrez_size=dim(y.Hs)[1]
# Annotating Gene-Symbol
y.Hs$genes$Symbol <- mapIds(org.Hs.eg.db, rownames(y.Hs), keytype="ENTREZID", column="SYMBOL")



# HCC1187
hcc1187.Hs<- y.Hs[,grep("^HCC1187", colnames(y))]
A <- factor(c(1,1,1,0,0,0))
design = model.matrix(~0 + A )
colnames(design)=c("C","A")
hcc1187.Hs <- estimateDisp(hcc1187.Hs , design, robust=TRUE)
plotBCV(hcc1187) #plot NB dispersion
fit_int <- glmQLFit(hcc1187.Hs , design, robust=TRUE)
AvsC=makeContrasts(A-C, levels=design)

#load PATHWAYS
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_H_v5p1.rdata"))

# H
hcc1187.idx.H <- ids2indices(Hs.H,id=rownames(hcc1187.Hs))
hcc1187.cam.H.AvsC <- camera(hcc1187.Hs, hcc1187.idx.H, design,inter.gene.cor=0.01, contrast=AvsC, sort = T)
colnames(hcc1187.cam.H.AvsC)=c("GENESET_SIZE", "HCC1187_AvsC_dir", "HCC1187_AvsC_fdr")

# C2
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c2_v5p1.rdata"))

hcc1187.idx.c2 <- ids2indices(Hs.c2,id=rownames(hcc1187.Hs))
hcc1187.cam.c2.AvsC <- camera(hcc1187.Hs, hcc1187.idx.c2, design,inter.gene.cor=0.01, contrast=AvsC, sort = T)
colnames(hcc1187.cam.c2.AvsC)=c("GENESET_SIZE", "HCC1187_AvsC_dir", "HCC1187_AvsC_fdr")

# HCC202
hcc202.Hs<- y.Hs[,grep("^HCC202", colnames(y))]
A <- factor(c(1,1,1,0,0,0))
design = model.matrix(~0 + A )
colnames(design)=c("C","A")
hcc202.Hs <- estimateDisp(hcc202.Hs , design, robust=TRUE)
plotBCV(hcc202) #plot NB dispersion
fit_int <- glmQLFit(hcc202.Hs , design, robust=TRUE)
AvsC=makeContrasts(A-C, levels=design)

#load PATHWAYS
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_H_v5p1.rdata"))

# H
hcc202.idx.H <- ids2indices(Hs.H,id=rownames(hcc202.Hs))
hcc202.cam.H.AvsC <- camera(hcc202.Hs, hcc202.idx.H, design,inter.gene.cor=0.01, contrast=AvsC, sort = T)
colnames(hcc202.cam.H.AvsC)=c("GENESET_SIZE", "HCC202_AvsC_dir", "HCC202_AvsC_fdr")

# C2
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c2_v5p1.rdata"))

hcc202.idx.c2 <- ids2indices(Hs.c2,id=rownames(hcc202.Hs))
hcc202.cam.c2.AvsC <- camera(hcc202.Hs, hcc202.idx.c2, design,inter.gene.cor=0.01, contrast=AvsC, sort = T)
colnames(hcc202.cam.c2.AvsC)=c("GENESET_SIZE", "HCC202_AvsC_dir", "HCC202_AvsC_fdr")


# HS578T
hs578t.Hs<- y.Hs[,grep("^HS578T", colnames(y))]
A <- factor(c(1,1,1,0,0,0))
design = model.matrix(~0 + A )
colnames(design)=c("C","A")
hs578t.Hs <- estimateDisp(hs578t.Hs , design, robust=TRUE)
plotBCV(hs578t) #plot NB dispersion
fit_int <- glmQLFit(hs578t.Hs , design, robust=TRUE)
AvsC=makeContrasts(A-C, levels=design)

#load PATHWAYS
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_H_v5p1.rdata"))

# H
hs578t.idx.H <- ids2indices(Hs.H,id=rownames(hs578t.Hs))
hs578t.cam.H.AvsC <- camera(hs578t.Hs, hs578t.idx.H, design,inter.gene.cor=0.01, contrast=AvsC, sort = T)
colnames(hs578t.cam.H.AvsC)=c("GENESET_SIZE", "HCC1187_AvsC_dir", "HCC1187_AvsC_fdr")

# C2
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c2_v5p1.rdata"))

hs578t.idx.c2 <- ids2indices(Hs.c2,id=rownames(hs578t.Hs))
hs578t.cam.c2.AvsC <- camera(hs578t.Hs, hs578t.idx.c2, design,inter.gene.cor=0.01, contrast=AvsC, sort = T)
colnames(hs578t.cam.c2.AvsC)=c("GENESET_SIZE", "HCC1187_AvsC_dir", "HCC1187_AvsC_fdr")
    

write.table(variabile, sep="\t", col.names=NA, file="nomefile.txt")


