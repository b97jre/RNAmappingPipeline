source("http://bioconductor.org/workflows.R")
workflowInstall("rnaseqGene")

library(rnaseqGene)
library(DESeq2)
library(reshape)
library(ggplot2)
library("pheatmap")
library("RColorBrewer")

setwd("/Users/johanreimegard/Vetenskap/Data/Wierup/Project_GBP_Pig")


colNamesInCountData = "Mapping_Name"
var1 = "Surgery"
var2 = "Location"
intGroups = c(var1,var2,"nr","Operation_date", "weight")

minRowSums = 10




#Load data (counts from HTseq)
exp.data <- read.table("Gene.Stranded.htseq.count.table.txt", header=TRUE, row.names=1, sep="\t")
#Load metaInfo
metaInfo <- read.table("MetaData.tab.table.txt", sep = "\t", header = TRUE);


# Make sure that the sample names are correct and in the same order in the meta data and the count data
sampleNames = substring(colnames(exp.data),2)
colnames(exp.data) <- sampleNames

exp.data = exp.data[,order(colnames(exp.data))]
metaInfo = metaInfo[order(metaInfo[colNamesInCountData]),]
rownames(metaInfo) = metaInfo$Sample
colnames(exp.data) = metaInfo$Sample


#Remove last five columns that are summaries in HTseq
exp.data = exp.data[1:((dim(exp.data)[1])-5),]


# add two variables in the end of metadata to make it applicable always
metaInfo['var1'] = metaInfo[var1]
metaInfo['var2'] = metaInfo[var2]


exp.data = exp.data[ ,-which(colnames(exp.data) == "Sample_821_ileum")]
metaInfo = metaInfo[colnames(exp.data), ]

############################################################################
###Creating DEseq2 object and normalising it for visualisation
###########################################################################
(dds <- DESeqDataSetFromMatrix(countData = exp.data,
                                  colData = metaInfo,
                                  design = ~ var1 + var2))




#Remove rows with low counts and normalise samples for visualisation
dds <- dds[ rowSums(counts(dds)) > minRowSums, ]
rld <- rlog(dds, blind=FALSE)


#########################################################################
##Visualising the data
###########################################################################

plotSample2SampleDistance(assay(rld))



# remove sample_821_ileum


values = assay(rld)
metaDataTable = metaInfo
n.comp = 5

(data <- plotPCA(rld, intgroup = intGroups, returnData=TRUE))
percentVar <- round(100 * attr(data, "percentVar"))

ggplot(data, aes(PC1, PC2, color=Surgery, shape=as.factor(Operation_date))) + geom_point(size=3) +
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





