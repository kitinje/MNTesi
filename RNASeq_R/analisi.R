# TAB files in /pico/scratch/userexternal/mfratell/MARTINA/RNASeq/TAB
# scp mfratell@login.pico.cineca.it:/pico/scratch/userexternal/mfratell/MARTINA/RNASeq/TAB/*.tab ./

setwd("/Users/age/Desktop/NUOVE/")

targets = dir(pattern = "*.tab")

temp = read.table("HCC1187_24H_A1_S4_001_ReadsPerGene.out.tab", sep="\t", skip = 4, row.names = 1, stringsAsFactors = F)
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


library(edgeR)
y <- DGEList(counts=RAWCOUNTS)

# remove everything after "." character
keep=grep(cbind(unlist(strsplit(rownames(y),'.',fixed=TRUE))), pattern = "ENSG")
rownames(y)<-cbind(unlist(strsplit(rownames(y),'.',fixed=TRUE)))[keep]
dim(y)

y <- calcNormFactors(y)

# keep files genes showing more than 0.5 cpm in more than 1 sample. 
exp_thr=0.5
samp_thr=3

keep <- rowSums(cpm(y)>exp_thr) >= samp_thr
y <- y[keep, , keep.lib.sizes=FALSE] #ricalcola anche library sizes (consigliato)

dim(y)

y <- calcNormFactors(y)



y$samples$group = factor( c( rep("ATRA",3), rep("CTRL", 3) ) )

# create design matrix
design <- model.matrix(~0 + y.Hs$samples$group)
colnames(design) <- c("A", "C")

y <- estimateDisp(y, design, robust=TRUE)
plotBCV(y)
fit <- glmQLFit(y, design, robust=TRUE) 
plotQLDisp(fit)

AvsC <- makeContrasts(A-C, levels=design)
y.AvsC <- glmQLFTest(fit, contrast=AvsC)

topTags(y.AvsC, n=10)$table
plot(y$counts["ENSG00000168447",])

y.de <- topTags(y.AvsC, n=Inf)$table
y.de_0.05 = y.de[which(y.de$FDR<0.05) ,]
y.de_up = y.de_0.05[which(y.de_0.05$logFC>0), ]
y.de_dn = y.de_0.05[which(y.de_0.05$logFC<0), ]


y.cpm = cpm(y)
y.cpm.log2 = cpm(y, log = TRUE, prior.count = 1)
hist(y.cpm.log2)
plot(density(y.cpm.log2))

# se vuoi distribuzione normale
y.voom <-voom(y, plot=TRUE)
hist(y.voom$E)
