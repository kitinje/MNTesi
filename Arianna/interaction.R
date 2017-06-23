setwd("/Users/age/Desktop/ARIANNA/RAW/")
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

y$samples$cell_line = c( rep("HCC1599",12), rep("MB157",12) )
y$samples$treatment = factor( rep(c(rep("A",3), rep("AD", 3), rep("C", 3), rep("D", 3)  ),2) )

y.cpm.log2 = cpm(y, log = TRUE, prior.count = 1)


pcaY <- prcomp(t(y.cpm.log2), scale. = TRUE)

pcaDf <- data.frame(pcaY$x, "cell_line" = y$samples$cell_line, "treatment" = y$samples$treatment)
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


mb157 <- y[,grep(x = colnames(y),  pattern = "MB157")]

y$samples

#create design matrix
A <- factor(c(1,1,1,1,1,1,0,0,0,0,0,0))
D <- factor(c(0,0,0,1,1,1,0,0,0,1,1,1))


design <- model.matrix(~ A + D + A:D, data=mb157$samples)

mb157 <- estimateDisp(mb157, design, robust=TRUE)
plotBCV(mb157) #plot NB dispersion

fit_int <- glmQLFit(mb157, design, robust=TRUE)

mb157.qlf.AvsC <- glmQLFTest(fit_int,  coef = 2)
resultsA = topTags(mb157.qlf.AvsC, n=Inf)$table
dim(resultsA[resultsA$FDR<0.05,])

mb157.qlf.DvsC <- glmQLFTest(fit_int,  coef = 3)
resultsD = topTags(mb157.qlf.DvsC, n=Inf)$table
dim(resultsD[resultsD$FDR<0.05,])

mb157.qlf.ADvsC <- glmQLFTest(fit_int,  coef = 4)
resultsAD = topTags(mb157.qlf.ADvsC, n=Inf)$table
dim(resultsAD[resultsAD$FDR<0.05,])

