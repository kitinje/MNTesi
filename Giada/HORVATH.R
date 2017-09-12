setwd("~/Desktop/MRC/temp/E-GEOD-84003.processed.1/")
targets = dir()

temp = c()
for(i in 1:length(targets)){
  temp = cbind(temp, read.table(file=targets[i], sep="\t", header = T, row.names = 1, stringsAsFactors = F)[,1])
}
colnames(temp) = targets
rownames(temp) = rownames(read.table(file=targets[i], sep="\t", header = T, row.names = 1, stringsAsFactors = F))
colnames(temp) = gsub(x = colnames(temp), pattern = "_sample_table.txt", replacement = "")

write.table(temp, sep="\t", file = "EGEOD84003.txt", col.names = NA)

EGEOD84003 = temp

source("https://bioconductor.org/biocLite.R")
biocLite("impute")
biocLite("preprocessCore")
biocLite("GO.db")
biocLite("WGCNA")
biocLite("sqldf")
biocLite("RPMM")

library(WGCNA)
library(sqldf)
library(impute)

source("~/Desktop/MRC/temp/NORMALIZATION.R")

trafo= function(x,adult.age=20) { 
  x=(x+1)/(1+adult.age) 
  y=ifelse(x<=1, log( x),x-1)
  y }

anti.trafo= function(x,adult.age=20) { 
  ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) 
}

setwd("~/Desktop/MRC/temp/")
probeAnnotation21kdatMethUsed=read.csv("probeAnnotation21kdatMethUsed.csv")
probeAnnotation27k=read.csv("datMiniAnnotation27k.csv")

datClock=read.csv("AdditionalFile3.csv")


#Read in the DNA methylation data (beta values)

# For a small file, e.g. measured on the 27k platform you could just use read.csv.
# But for large files, e.g. those measured on the 450K platform, I recommend you use read.csv.sql.




#dat0=read.csv("MethylationDataExample55.csv") ;
dat0 = read.table(file = "EGEOD84003.txt", sep="\t", header = T)
nSamples=dim(dat0)[[2]]-1
nProbes= dim(dat0)[[1]]

XchromosomalCpGs=as.character(probeAnnotation27k$Name[probeAnnotation27k$Chr=="X"])
selectXchromosome=is.element(dat0[,1], XchromosomalCpGs )
selectXchromosome[is.na(selectXchromosome)]=FALSE
meanXchromosome=rep(NA, dim(dat0)[[2]]-1)
if ( sum(selectXchromosome) >=500 ) {
  meanXchromosome= as.numeric(apply( as.matrix(dat0[selectXchromosome,-1]),2,mean,na.rm=TRUE)) }
if ( sum(is.na(meanXchromosome)) >0 ) { cat(paste( "\n \n Comment: There are lots of missing values for X chromosomal
                                                   probes for some of the samples. This is not a problem when it comes to estimating age but I cannot predict the gender of these
                                                   samples.\n " ),file="LogFile.txt",append=TRUE) } 

#STEP 2: Restrict the data to 21k probes and ensure they are numeric
match1=match(probeAnnotation21kdatMethUsed$Name , dat0[,1])
if ( sum( is.na(match1))>0 ) stop(paste(sum( is.na(match1)), "CpG probes cannot be matched"))
dat1= dat0[match1,]
asnumeric1=function(x) {as.numeric(as.character(x))}
dat1[,-1]=apply(as.matrix(dat1[,-1]),2,asnumeric1)

set.seed(1)
# Do you want to normalize the data (recommended)?
normalizeData=TRUE
source("StepwiseAnalysis.txt")

datout

normalizeData=FALSE
source("StepwiseAnalysis.txt")

datout




