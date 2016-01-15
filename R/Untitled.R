library(edgeR)
setwd("/Users/johanreimegard/Vetenskap/Data/Wierup/Project_GBP_Pig")

exp.data <- read.table("Gene.Stranded.htseq.count.table.txt", header=TRUE, row.names=1, sep="\t")
metaInfo <- read.table("MetaData.tab.table.txt", sep = "\t", header = TRUE);

substr(colnames(exp.data),start = 2,last = length(colnames(exp.data)))

#temp becuse X405_ileum_TGACCA_L001 did not work
exp.data = exp.data[, 1:18]

exp.data <- exp.data + 1
lib.size <- apply(exp.data,2,sum)
scale.factors <- calcNormFactors(exp.data, method="TMM")
norm.data <- t(t(exp.data)/(scale.factors*lib.size))
norm.data <- log(norm.data)
hist(rowSums(norm.data),breaks = 50)

norm.data.cleared <- norm.data[rowSums(norm.data) >-250 , ]
hist(rowSums(norm.data.cleared),breaks = 50)





