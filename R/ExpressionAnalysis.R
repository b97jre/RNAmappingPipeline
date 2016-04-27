source("http://bioconductor.org/workflows.R")
workflowInstall("rnaseqGene")

library(rnaseqGene)
library(DESeq2)
library(reshape)
library(ggplot2)
library("pheatmap")
library("RColorBrewer")

source("/Users/johanreimegard/git/GBPanalysis/R/ExpressionAnalysisFunctions.R")

wd = "/Users/johanreimegard/Vetenskap/Data/jakobsson/nobackup/mRNA/PD/"



colNamesInCountData = "sampleName"
var1 = "Days"
var2 = "AAV.a.syn"
intGroups = c(var1,var2)

minRowSums = 10



setwd(wd)

#Load data (counts from HTseq)
exp.data <- read.table("Gene.Stranded.htseq.count.table.txt", header=TRUE, row.names=1, sep="\t",comment.char="")
#Load metaInfo
metaInfo <- read.table("MetaData.table.txt", sep = "\t", header = TRUE,comment.char="");


head(exp.data)
head(metaInfo)
# Make sure that the sample names are correct and in the same order in the meta data and the count data
sampleNames = substring(colnames(exp.data),2)
colnames(exp.data) <- sampleNames
colNamesInCountData = sampleNames

exp.data = exp.data[,order(colnames(exp.data))]
metaInfo = metaInfo[order(metaInfo[colNamesInCountData]),]
rownames(metaInfo) = metaInfo[[colNamesInCountData]]
colnames(exp.data) = metaInfo[[colNamesInCountData]]

#Remove last five columns that are summaries in HTseq
exp.data = exp.data[1:((dim(exp.data)[1])-5),]


# Check only
mhTT.metaInfo = metaInfo[grep("mHTT", rownames(metaInfo)),]
yn.metaInfo = metaInfo[grep("yn", rownames(metaInfo)),]


mhTT.exp.data = exp.data[,grep("mHTT", colnames(exp.data))]
head(mhTT.exp.data)
write.table(mhTT.exp.data, file = "mHTT/Gene.Stranded.htseq.count.table.txt",quote = FALSE, row.names = TRUE,col.names = TRUE , sep = "\t")

yn.exp.data = exp.data[,grep("yn", colnames(exp.data))]
head(yn.exp.data)
write.table(yn.exp.data, file = "PD/Gene.Stranded.htseq.count.table.txt",quote = FALSE, row.names = TRUE,col.names = TRUE , sep = "\t")

# add two variables in the end of metadata to make it applicable always
metaInfo['var1'] = metaInfo[var1]
metaInfo['var2'] = metaInfo[var2]





head(exp.data)

############################################################################
###Creating DEseq2 object and normalising it for visualisation
###########################################################################
(dds <- DESeqDataSetFromMatrix(countData = exp.data,
                                  colData = metaInfo,
                                  design = ~ var1 + var2))




#Remove rows with low counts and normalise samples for visualisation
dds <- dds[ rowSums(counts(dds)) > minRowSums, ]
rld <- rlog(dds, blind=FALSE)

dds <- estimateSizeFactors(dds)
countDataMatrix <-counts(dds, normalized=TRUE)


#########################################################################
##Visualising the data
###########################################################################

plotSample2SampleDistance(assay(rld))



values = assay(rld)
metaDataTable = metaInfo
n.comp = 5

data <- plotPCA(rld, intgroup = intGroups)

(data <- plotPCA(rld, intgroup = intGroups, returnData=TRUE))
percentVar <- round(100 * attr(data, "percentVar"))

ggplot(data, aes(PC1, PC2, color=lane, shape=as.factor(flowcell))) + geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))



mdsData <- data.frame(cmdscale(sampleDistMatrix))
mds <- cbind(mdsData, as.data.frame(colData(rld)))
ggplot(mds, aes(X1,X2,color=var1,shape=var2,label = nr)) + geom_point(size=3) + geom_text()

par( mfrow = c( 1, 2 ) )
dds <- estimateSizeFactors(dds)
plot(log2(counts(dds, normalized=TRUE)[,1:2] + 1),
     pch=16, cex=0.3)
plot(assay(rld)[,1:2],
     pch=16, cex=0.3)


par( mfrow = c( 1, 2 ) )
dds <- estimateSizeFactors(dds)
plot(log2(counts(dds, normalized=TRUE)[,1:2] + 1),
     pch=16, cex=0.3)
plot(assay(rld)[,1:2],
     pch=16, cex=0.3)


d = DGEList(counts=exp.data.matrix)
d <- calcNormFactors(d, method="TMM"  )
d$samples

cpmN  <- cpm(d, normalized.lib.sizes=TRUE ,log = TRUE )

head(cpmN)


RsquareValue = matrix(0,ncol = nrOfSamples,nrow = nrOfSamples)
colnames(RsquareValue) <- sampleNames
rownames(RsquareValue) <- sampleNames
for(i in 1:(nrOfSamples-1)){
  RsquareValue[i,i] = 1
  for(j in (i+1):nrOfSamples){
    lmod = lm(cpmN[,i]~cpmN[,j])
    R2 = summary(lmod)$r.squared
    RsquareValue[i,j]= R2
    RsquareValue[j,i] = R2
  }
}
RsquareValue[i+1,i+1] = 1

max(RsquareValue)
dftest = melt(RsquareValue)
colnames(dftest) <- c("sampleX","sampley","R2")
Rsquare = ggplot(dftest, aes(sampleX, sampley, fill = R2)) +
  geom_raster() +
  ggtitle("R-square values") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

Rsquare
qplot(sampleX, sampley, data = dftest, fill = z, geom = "raster", title = "test")




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





