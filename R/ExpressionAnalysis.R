library(edgeR)
library(reshape)
library(ggplot2)
setwd("/Users/johanreimegard/Vetenskap/Data/Wierup/Project_GBP_Pig")


exp.data <- read.table("Gene.Stranded.htseq.count.table.txt", header=TRUE, row.names=1, sep="\t")
metaInfo <- read.table("MetaData.tab.table.txt", sep = "\t", header = TRUE);

# If colnames start with a number
sampleNames = substring(colnames(exp.data),2)
colnames(exp.data) <- sampleNames

#temp becuse X405_ileum_TGACCA_L001 did not work
exp.data.matrix = as.matrix(exp.data[1:((dim(exp.data)[1])-5),])
nrOfSamples = length(sampleNames)

d = DGEList(counts=exp.data.matrix)
d <- calcNormFactors(d, method="TMM"  )
d$samples

cpmN  <- cpm(d, normalized.lib.sizes=TRUE ,log = TRUE )

head(cpmN)


RsquareValue = matrix(0,ncol = nrOfSamples,nrow = nrOfSamples)
colnames(RsquareValue) <- sampleNames
rownames(RsquareValue) <- sampleNames
for(i in 1:(nrOfSamples-1)){
  for(j in (i+1):nrOfSamples){
    lmod = lm(cpmN[,i]~cpmN[,j])
    R2 = summary(lmod)$r.squared
    RsquareValue[i,j]= R2
    RsquareValue[j,i] = R2
  }
}

max(RsquareValue)
dftest = melt(RsquareValue)
colnames(dftest) <- c("sampleX","sampley","z")

(dftest) <- sampleNames
lmod = lm(cpmN[,2]~cpmN[,6])
summary(lmod)$r.squared
plot(colSums(cpmN))
colSums(cpmN)


class(cpmN)
head(cpm)
apply(cpm ,2,sum)
norm.cpm

norm.cpm <- t(t(exp.data.matrix)/(scale.factors*lib.size))
hist(rowSums(log(norm.cpm)),breaks = 50)
norm.data <- log(norm.data)
hist(rowSums(norm.data),breaks = 50)

norm.data.cleared <- norm.data[rowSums(norm.data) >-250 , ]
hist(rowSums(norm.data.cleared),breaks = 50)





