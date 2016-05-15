##########################################
## Plot counts of reads in different biotypes (miRNA, tRNA etc.),
## and plot reads lengths.


## Arguments
args <- commandArgs(trailingOnly = TRUE)
input.dir <- args[1]
plot.header <- args[2]
out.plot.file <- args[3]
out.tab.file <- args[4]

## Hard coded values
suffix <- "_lengths.txt"
files =  list.files(path = input.dir, full.names = FALSE)
biotypesFiles = files[grep(pattern = suffix, x = files)]
biotypes  = gsub(pattern = suffix,replacement = "", x = biotypesFiles)
lwd=3
xlim=c(15,48)
show.lengths <- as.character(xlim[1]:xlim[2])


##################
## Read data, allow missing values (which will be set to 0)
input.files <- list.files(path=input.dir, pattern=paste("*", suffix, sep=""), full.names=TRUE)


length.count.mat <- matrix(0, ncol=length(biotypes), nrow=length(show.lengths))
rownames(length.count.mat) <- show.lengths
colnames(length.count.mat) <- gsub(suffix, "", sapply(input.files, basename))
for (input.file in input.files){
	print(file=stderr(), input.file)
	tmp.data <- read.table(input.file, as.is=TRUE)
	new.lengths <- tmp.data[,1]
	names(new.lengths) <- tmp.data[,2]
	new.lengths <- new.lengths[intersect(show.lengths, names(new.lengths))]
	
	new.biotype <- gsub(suffix, "", basename(input.file))
	length.count.mat[names(new.lengths),new.biotype] <- new.lengths
}

## Sort matrix
if(ncol(length.count.mat) > 1){
  length.count.mat <- length.count.mat[, biotypes]
}


## Print table with combined data
write.table(length.count.mat, file=out.tab.file, sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)


################
## Make plots
pdf(out.plot.file, width=7, height=7)
par(mfrow=c(2,2), oma = c(0, 0, 2, 0))

## 1) Plot overall length distribution
if(ncol(length.count.mat) > 1){
  length.distr <- apply(length.count.mat,1,sum)
}else{
  length.distr <- length.count.mat[,1]
}
length.distr <- length.distr[show.lengths]
plot(x=names(length.distr), y=length.distr, type="l", lwd=lwd, col="gray",
		 xlab="nucleotides", ylab="nr reads")
mtext(plot.header, outer = TRUE, cex = 1.5)

## 2) Plot biotype distribution
require(grDevices)
biotypes <- apply(length.count.mat,2, sum)
cols <- rainbow(ncol(length.count.mat))

barplot(biotypes, las=3, col = cols, ylab="nr reads")


## 3) Length distribution within each biotype
cols <- rainbow(ncol(length.count.mat))
legend.cex=0.6
if(ncol(length.count.mat) > 1){
  plot.data <- length.count.mat[as.character(xlim[1]:xlim[2]),]
}else{
  plot.data <- length.count.mat
}



plot(x=rownames(plot.data), y=plot.data[,1], type="l", col=cols[1],
     lwd=lwd, xlim=xlim, ylim=c(0, max(plot.data)), xlab="nucleotides", ylab="nr reads")
if(ncol(length.count.mat) > 1){
  for (i in 2:ncol(plot.data)){
    lines(x=rownames(plot.data), y=plot.data[,i], col=cols[i], lwd=lwd)
  }
}
legend(x="topright", fill=cols[1:ncol(plot.data)], legend=colnames(plot.data), cex=legend.cex)


## 4) Length distribution within each biotype
plot.data <- t(t(plot.data)/apply(plot.data, 2, sum)) ## Normalize for each biotype

plot(x=rownames(plot.data), y=plot.data[,1], type="l", col=cols[1],
     lwd=lwd, xlim=xlim, ylim=c(0, max(plot.data)), xlab="nucleotides", ylab="fraction reads")
if(ncol(length.count.mat) > 1){
  for (i in 2:ncol(plot.data)){
    lines(x=rownames(plot.data), y=plot.data[,i], col=cols[i], lwd=lwd)
  }
}
legend(x="topright", fill=cols[1:ncol(plot.data)], legend=colnames(plot.data), cex=legend.cex)

dev.off()


