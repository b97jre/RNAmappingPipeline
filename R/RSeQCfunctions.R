
library(reshape2)
library(ggplot2)


LinetypeValues = rep(c(rep("solid", 9), rep("dashed", 9),rep("dotted", 9), rep("dotdash", 9),rep("longdash", 9), rep("twodash", 9)), 40)
colorValues =rep( rep(RColorBrewer::brewer.pal(9,"Set1"),6), 40)

plotValues = data.frame(linetype = LinetypeValues, color = colorValues)


##################################################################################main function


getAllPlots <- function(RSeQCdir, metaInfo){

  # Getting read distribution info
  Groups = c("TES_down_10kb","TES_down_5kb","TES_down_1kb","TSS_up_10kb","TSS_up_5kb","TSS_up_1kb","Introns","3'UTR_Exons","5'UTR_Exons","CDS_Exons")
  RD <- data.frame(Group= factor(levels = Groups ) ,Total_bases = integer(0),Tag_count=integer(0) ,Tags_Kb = numeric(0), TPKM = numeric(0),sampleName = character(0) )
  for(i in 1: length(metaInfo$sampleName)){
    fileName = paste(metaInfo$sampleName[i],"read_distribution.txt", sep = '.')
    RD <-rbind(RD,read_distributionFile(RSeQCdir,fileName,sampleName = metaInfo$sampleName[i]))
  }



  # Getting TIN info
  TIN <- data.frame(Group= character(0) ,chrom = character(0),tx_start=integer(0) ,tx_end = integer(0), TIN = numeric(0),sampleName = character(0))
  for(i in 1: length(metaInfo$sampleName)){
    fileName = paste(metaInfo$sampleName[i],"tin.xls", sep = '.')
    TIN <-rbind(TIN,readTINFile(RSeQCdir,fileName,sampleName = metaInfo$sampleName[i]))
  }


  CP <- data.frame(position = integer(0), R1  = integer(0), R1percentage = numeric(0) ,R2 = integer(0), R2percentage = numeric(0),sampleName = character(0))
  for(i in 1: length(metaInfo$sampleName)){
    fileName = paste(metaInfo$sampleName[i],"clipping_profile.r", sep = '.')
    CP <-rbind(CP,readClippingProfile(RSeQCdir,fileName,sampleName = metaInfo$sampleName[i]))
  }


  ID <- data.frame(from = integer(0), position  = integer(0),count = integer(0), Percent = numeric(0),sampleName = character(0))
  for(i in 1: length(metaInfo$sampleName)){
    fileName = paste(metaInfo$sampleName[i],"inner_distance_freq.txt", sep = '.')
    ID <-rbind(ID,readInnerDistance(RSeQCdir,fileName,sampleName = metaInfo$sampleName[i]))
  }



  GBC <- data.frame(position  = integer(0), count = integer(0), Percent = numeric(0),sampleName = character(0))
  for(i in 1: length(metaInfo$sampleName)){
    fileName = paste(metaInfo$sampleName[i],"geneBodyCoverage.txt", sep = '.')
    GBC <-rbind(GBC,readGeneBodyCoverage(RSeQCdir,fileName,sampleName = metaInfo$sampleName[i]))
  }

  JS <- data.frame(position  = integer(0), knownJunctions = integer(0), AllJunctions = integer(0), NovelJunctions = integer(0), knownJuncionsFraction = numeric(0),NovelJunctionsFraction = numeric(0),sampleName = character(0))
  for(i in 1: length(metaInfo$sampleName)){
    fileName = paste(metaInfo$sampleName[i],"junctionSaturation_plot.r", sep = '.')
    JS <-rbind(JS,read_junctionSaturation(RSeQCdir,fileName,sampleName = metaInfo$sampleName[i]))
  }


  TPKMplot <- plotTPKMInfo(tagDistribution =  RD)
  TPKMplots = list(TPKMplot = TPKMplot)

  TINplot <- plotTINInfo(TIN)
  TINplots = list(TINplot = TINplot)

  ClippingInfoPlot <- plotClippingInfo(CP)
  ClippingInfoPlot2 <- plotClippingInfo2(CP)
  ClippingInfoPlots = list(ClippingInfoPlot = ClippingInfoPlot)
  ClippingInfoPlots2 = list(ClippingInfoPlot2 = ClippingInfoPlot2)

  InnerDistancePlot <- plotInnerDistance(ID)
  InnerDistancePlot2 <- plotInnerDistance2(ID)

  InnerDistancePlots1 = list(InnerDistancePlot = InnerDistancePlot)
  InnerDistancePlots2 = list(InnerDistancePlot2 = InnerDistancePlot2)
  GBCplot <- plotGeneBodyCoverage(GBC)
  GBCplot2 <- plotGeneBodyCoverage2(GBC)
  GBCplots = list(GBCplot = GBCplot, GBCplot2 = GBCplot2)
  junctionSaturationPlot <- plotjunctionSaturation(JS)
  junctionSaturationPlot2 <- plotjunctionSaturation2(JS)
  junctionSaturationPlots = list(junctionSaturationPlot = junctionSaturationPlot, junctionSaturationPlot2 = junctionSaturationPlot2)

  Allplots = list(list(TPKMplot), list(TINplot), list(ClippingInfoPlot) ,
                  lisjat(ClippingInfoPlot2), list(InnerDistancePlot),
                  list(InnerDistancePlot2),list(GBCplot), list(GBCplot2),
                  list(junctionSaturationPlot) , list(junctionSaturationPlot2))

  plotAllPlotsInOnePDf(Allplots, "test.pdf")

}

getTPKMPlots <- function(RSeQCdir, metaInfo){

  # Getting read distribution info
  Groups = c("TES_down_10kb","TES_down_5kb","TES_down_1kb","TSS_up_10kb","TSS_up_5kb","TSS_up_1kb","Introns","3'UTR_Exons","5'UTR_Exons","CDS_Exons")
  RD <- data.frame(Group= factor(levels = Groups ) ,Total_bases = integer(0),Tag_count=integer(0) ,Tags_Kb = numeric(0), TPKM = numeric(0),sampleName = character(0) )
  for(i in 1: length(metaInfo$sampleName)){
    fileName = paste(metaInfo$sampleName[i],"read_distribution.txt", sep = '.')
    RD <-rbind(RD,read_distributionFile(RSeQCdir,fileName,sampleName = metaInfo$sampleName[i]))
  }


  TPKMplot <- plotTPKMInfo(tagDistribution =  RD)
}



getTINPlots <- function(RSeQCdir, metaInfo){

  # Getting TIN info
  TIN <- data.frame(Group= character(0) ,chrom = character(0),tx_start=integer(0) ,tx_end = integer(0), TIN = numeric(0),sampleName = character(0))
  for(i in 1: length(metaInfo$sampleName)){
    fileName = paste(metaInfo$sampleName[i],"tin.xls", sep = '.')
    TIN <-rbind(TIN,readTINFile(RSeQCdir,fileName,sampleName = metaInfo$sampleName[i]))
  }

    TINplot <- plotTINInfo(TIN)
}


getClippingProfilePlots <- function(RSeQCdir, metaInfo){
 CP <- data.frame(position = integer(0), R1  = integer(0), R1percentage = numeric(0) ,R2 = integer(0), R2percentage = numeric(0),sampleName = character(0))
 for(i in 1: length(metaInfo$sampleName)){
   fileName = paste(metaInfo$sampleName[i],"clipping_profile.r", sep = '.')
   CP <-rbind(CP,readClippingProfile(RSeQCdir,fileName,sampleName = metaInfo$sampleName[i]))
 }
 ClippingInfoPlot <- plotClippingInfo(CP)
 ClippingInfoPlot2 <- plotClippingInfo2(CP)
 ClippingInfoPlots = list(ClippingInfoPlot = ClippingInfoPlot,ClippingInfoPlot2 = ClippingInfoPlot2)
}




getInnerDistancePlot<- function(RSeQCdir, metaInfo){
  ID <- data.frame(from = integer(0), position  = integer(0),count = integer(0), Percent = numeric(0),sampleName = character(0))
  for(i in 1: length(metaInfo$sampleName)){
    fileName = paste(metaInfo$sampleName[i],"inner_distance_freq.txt", sep = '.')

    ID <-rbind(ID,readInnerDistance(RSeQCdir,fileName,sampleName = metaInfo$sampleName[i]))
  }
  InnerDistancePlot <- plotInnerDistance(ID)
  InnerDistancePlot2 <- plotInnerDistance2(ID)

  InnerDistancePlots1 = list(InnerDistancePlot = InnerDistancePlot)
  InnerDistancePlots2 = list(InnerDistancePlot2 = InnerDistancePlot2)
}

getGeneBodyCoveragePlots<- function(RSeQCdir, metaInfo){

  GBC <- data.frame(position  = integer(0), count = integer(0), Percent = numeric(0),sampleName = character(0))
  for(i in 1: length(metaInfo$sampleName)){
    fileName = paste(metaInfo$sampleName[i],"geneBodyCoverage.txt", sep = '.')
    GBC <-rbind(GBC,readGeneBodyCoverage(RSeQCdir,fileName,sampleName = metaInfo$sampleName[i]))
  }

  GBCplot <- plotGeneBodyCoverage(GBC)
  GBCplot2 <- plotGeneBodyCoverage2(GBC)
  GBCplots = list(GBCplot = GBCplot, GBCplot2 = GBCplot2)
}



getJunctionSaturationPlots<- function(RSeQCdir, metaInfo){

  JS <- data.frame(position  = integer(0), knownJunctions = integer(0), AllJunctions = integer(0), NovelJunctions = integer(0), knownJuncionsFraction = numeric(0),NovelJunctionsFraction = numeric(0),sampleName = character(0))
  for(i in 1: length(metaInfo$sampleName)){
    fileName = paste(metaInfo$sampleName[i],"junctionSaturation_plot.r", sep = '.')
    JS <-rbind(JS,read_junctionSaturation(RSeQCdir,fileName,sampleName = metaInfo$sampleName[i]))
  }

  junctionSaturationPlot <- plotjunctionSaturation(JS)
  junctionSaturationPlot2 <- plotjunctionSaturation2(JS)
  junctionSaturationPlots = list(junctionSaturationPlot = junctionSaturationPlot, junctionSaturationPlot2 = junctionSaturationPlot2)

}



#################################################### Read functions ############################################

read_distributionFile <- function(dir, fileName, sampleName){
  readDistributionInfo = readLines(paste(dir,fileName, sep="/"))
  nrOfTags = as.numeric(strsplit(x = readDistributionInfo[3] , split = "           ")[[1]][2])
  tagDistribution = read.table(text = readDistributionInfo[5:15],header = TRUE)
  tagDistribution['TPKM'] = tagDistribution$Tags.Kb/(nrOfTags/1000000)
  tagDistribution['sampleName'] = sampleName
  return (tagDistribution)
}

readTINFile <- function(dir, fileName, sampleName){
  TINinfo = read.table(paste(dir,fileName, sep="/"), sep = "\t", header = TRUE)
  TINinfo['sampleName'] = sampleName
  return(TINinfo)
}



readClippingProfile <- function(dir, fileName, sampleName){
  clipping.profile = read.table(paste(dir,fileName, sep="/"), sep = "\t", header = FALSE)

  clippingPositions =regmatches(clipping.profile[2,1], gregexpr("(?<=\\().*?(?=\\))", clipping.profile[2,1], perl=T))[[1]]


  clippingInfo = as.data.frame(as.numeric (strsplit (clippingPositions,"," )[[1]]) )
  junctionR1= regmatches(clipping.profile[3,1], gregexpr("(?<=\\().*?(?=\\))", clipping.profile[3,1], perl=T))[[1]]
  clippingInfo$R1 = as.numeric(strsplit (junctionR1,"," )[[1]])

  TotalCount = strsplit(as.character(clipping.profile[4,1]), "=")[[1]][[2]]
  TotalCount = as.numeric(strsplit(TotalCount,"-")[[1]][[1]])
  clippingInfo$R1percentage = clippingInfo$R1/TotalCount

  #If paired end
  if(nrow(clipping.profile) > 8){

    junctionR2= regmatches(clipping.profile[9,1], gregexpr("(?<=\\().*?(?=\\))", clipping.profile[9,1], perl=T))[[1]]
    clippingInfo$R2 = as.numeric(strsplit (junctionR2,"," )[[1]])

    TotalCount = strsplit(as.character(clipping.profile[10,1]), "=")[[1]][[2]]
    TotalCount = as.numeric(strsplit(TotalCount,"-")[[1]][[1]])
    clippingInfo$R2percentage = clippingInfo$R2/TotalCount
  } else{ # create info for consistency
    clippingInfo$R2 = 0
    clippingInfo$R2percentage = 0
  }


  colnames(clippingInfo) <-  c("position","R1","R1percentage","R2","R2percentage")
  clippingInfo['sampleName'] = sampleName

  return (clippingInfo)
}



readInnerDistance <- function(dir, fileName, sampleName){
  inner.distance = read.table(paste(dir,fileName, sep="/"), sep = "\t", header = FALSE)
  colnames(inner.distance) <- c("from","position","count")
  inner.distance['Percent'] = inner.distance$count /max(inner.distance$count)


  inner.distance['sampleName'] = sampleName
  return (inner.distance)
}

readGeneBodyCoverage <- function(dir, fileName, sampleName){
  GBC =t(data.matrix(read.table(paste(dir,fileName, sep="/"), sep = "\t", header = FALSE)))
  GBCinfo = data.frame(position = GBC[2:nrow(GBC),1], count = GBC[2:nrow(GBC),2])
  GBCinfo['Percent'] = GBCinfo$count /max(GBCinfo$count)
  GBCinfo['sampleName'] = sampleName

  return (GBCinfo)
}



read_junctionSaturation <- function(dir, fileName, sampleName){
  JSlines = readLines(paste(dir,fileName, sep="/"))
  AP = regmatches(JSlines[2], gregexpr("(?<=\\().*?(?=\\))", JSlines[2], perl=T))[[1]]
  position = as.numeric (strsplit (AP,"," )[[1]])
  JSinfo = data.frame(position)
  KJ = regmatches(JSlines[3], gregexpr("(?<=\\().*?(?=\\))", JSlines[3], perl=T))[[1]]
  JSinfo['knownJunctions'] = as.numeric (strsplit (KJ,"," )[[1]])
  KJ = regmatches(JSlines[4], gregexpr("(?<=\\().*?(?=\\))", JSlines[4], perl=T))[[1]]
  JSinfo['AllJunctions'] = as.numeric (strsplit (KJ,"," )[[1]])
  KJ = regmatches(JSlines[5], gregexpr("(?<=\\().*?(?=\\))", JSlines[5], perl=T))[[1]]
  JSinfo['NovelJunctions'] = as.numeric (strsplit (KJ,"," )[[1]])
  JSinfo['knownJuncionsFraction'] =JSinfo['knownJunctions'] / JSinfo['AllJunctions']
  JSinfo['NovelJunctionsFraction'] =JSinfo['NovelJunctions'] / JSinfo['AllJunctions']
  JSinfo['sampleName'] = sampleName
  return (JSinfo)
}




#########################################################Plot functions##################################

plotTPKMInfo <- function(tagDistribution){
  p <- ggplot(tagDistribution, aes(x=factor(sampleName),y = TPKM, fill=factor(Group))) +
    geom_bar(position="fill",stat="identity") +
    labs(title = "Tags Per 1000 bases per Million tags (TPKM)" , x = "Sample", y = "Fraction of TPKM") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+ coord_flip()

  return (p)
}




plotTINInfo <- function(TINinfo){
  plotTINInfoData = TINinfo[TINinfo$TIN >1, ]
  g <- ggplot(plotTINInfoData, aes(as.factor(sampleName),TIN ,colour = sampleName))+
    geom_boxplot()+
    labs(title = "Distribution of TIN values across samples" , x = "Sample", y = "TIN value") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="none")+ coord_flip()


}

plotTINInfo <- function(TINinfo){
  plotTINInfoData = TINinfo[TINinfo$TIN >1, ]
  g <- ggplot(plotTINInfoData, aes(as.factor(sampleName),TIN ,colour = sampleName))+
    geom_boxplot()+
    labs(title = "Distribution of TIN values across samples" , x = "Sample", y = "TIN value") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="none")+ coord_flip()


}



stat_sum_df <- function(fun, geom="crossbar",colour="black", ...) {
  stat_summary_bin(fun.data=fun, colour=colour, geom=geom, width = 0.90 , ...)
}




plotClippingInfo <- function(clippingInfo){
  clippingInfoPlot = clippingInfo[, c("position","sampleName","R1percentage","R2percentage")]
  colnames(clippingInfoPlot) <-c("position","sampleName","Forward read","Reverse read")
  clippingInfoPlot = melt(clippingInfoPlot,id.vars = c("position","sampleName"))
  nrOfSamples = length(unique(clippingInfoPlot$sampleName))

  plot <- ggplot(clippingInfoPlot, aes(x = position, y = value))
  plot = plot +  geom_line(aes(colour = sampleName))+ theme(legend.position="left" ) + guides(colour = guide_legend(reverse=TRUE))
  # To complicated and did not work properly maybe use later
  #  plot = plot +  geom_line(aes(colour = sampleName, linetype = sampleName))
  #  plot = plot + scale_linetype_manual(values = LinetypeValues[1:nrOfSamples])
  #  plot = plot + scale_color_manual(values = colorValues[1:nrOfSamples])

  #plot =  plot + stat_sum_df("mean_cl_normal", geom = "smooth")
  plot = plot + labs(title = "Distribution of clipped reads across samples" , x = "Position on read", y = "Fraction clipped reads") +facet_grid(variable ~ . )
  return (plot)
}

plotClippingInfo2 <- function(clippingInfo){
  clippingInfoPlot = clippingInfo[, c("position","sampleName","R1percentage","R2percentage")]
  colnames(clippingInfoPlot) <-c("position","sampleName","Forward read","Reverse read")
  clippingInfoPlot = melt(clippingInfoPlot,id.vars = c("position","sampleName"))
  nrOfSamples = length(unique(clippingInfoPlot$sampleName))

  plot <- ggplot(clippingInfoPlot, aes(x = position, fill = value, y = sampleName))
  plot = plot + geom_tile()+ scale_fill_gradient(low="white", high="blue")
  plot = plot + labs(title = "Distribution of clipped reads across samples" , x = "Position on read", y = "Sample") +facet_grid(variable ~ . )
  return (plot)
}

plotClippingInfo2Single <- function(clippingInfo){
  clippingInfoPlot = clippingInfo[, c("position","sampleName","R1percentage")]
  colnames(clippingInfoPlot) <-c("position","sampleName","Forward read")
  clippingInfoPlot = melt(clippingInfoPlot,id.vars = c("position","sampleName"))
  nrOfSamples = length(unique(clippingInfoPlot$sampleName))

  plot <- ggplot(clippingInfoPlot, aes(x = position, fill = value, y = sampleName))
  plot = plot + geom_tile()+ scale_fill_gradient(low="white", high="blue")
  plot = plot + labs(title = "Distribution of clipped reads across samples" , x = "Position on read", y = "Sample")
  return (plot)
}





plotInnerDistance <- function(inner.distance){
  plot = ggplot(inner.distance, aes(x = position, y = Percent, colour = sampleName))+geom_line()
  plot = plot + labs(title = "Distance between reads" , x = "Distance", y ="Fraction")+ theme(legend.position="left" ) + guides(colour = guide_legend(reverse=TRUE))
  return (plot)
}

plotInnerDistance2 <- function(inner.distance){
  plot = ggplot(inner.distance, aes(x = position, fill = Percent, y = sampleName))
  plot = plot + geom_tile()+ scale_fill_gradient(low="white", high="blue")
  plot = plot + labs(title = "Distance between reads" , x = "Distance", y ="Sample")

  return (plot)
}




plotGeneBodyCoverage <- function(GBCinfo){
  plot = ggplot(GBCinfo, aes(x = position, y = Percent, colour = sampleName))
  plot =  plot  +geom_line() + theme(legend.position="left" ) + guides(colour = guide_legend(reverse=TRUE))
  #plot =  plot + stat_sum_df("mean_cl_normal", geom = "smooth")
  plot = plot + labs(title = "Gene body coverage" , x = "5' -> 3'", y ="Relative coverage")
  return (plot)
}

plotGeneBodyCoverage2 <- function(GBCinfo){
  plot = ggplot(GBCinfo, aes(x = position, y = sampleName, fill = Percent))
  plot = plot + geom_tile() + scale_fill_gradient(low="white", high="blue")
  plot = plot + labs(title = "Gene body coverage" , x = "5' -> 3'", y ="Sample name")

  return (plot)
}



plotjunctionSaturation <- function(junctionSaturation){
  junctionSaturationPlot = junctionSaturation[, c("position","sampleName","knownJuncionsFraction")]
  colnames(junctionSaturationPlot) <-c("PercentOfReads","sampleName","Fraction known")
  junctionSaturationPlot = melt(junctionSaturationPlot,id.vars = c("PercentOfReads","sampleName"))

  plot <- ggplot(junctionSaturationPlot, aes(x = PercentOfReads, y = value))
  plot = plot +  geom_line(aes(colour = sampleName)) + theme(legend.position="left" ) + guides(colour = guide_legend(reverse=TRUE))
  #plot =  plot + stat_sum_df("mean_cl_normal", geom = "smooth")
  plot = plot + labs(title = "Distribution of known junctions depending on percent of reads used" , x = "Percent of reads", y = "Fraction") +facet_grid(variable ~ . )
  return (plot)
}

plotjunctionSaturation2 <- function(junctionSaturation){
  junctionSaturationPlot = junctionSaturation[, c("position","sampleName","knownJuncionsFraction")]
  colnames(junctionSaturationPlot) <-c("PercentOfReads","sampleName","Fraction known")
  junctionSaturationPlot = melt(junctionSaturationPlot,id.vars = c("PercentOfReads","sampleName"))

  plot <- ggplot(junctionSaturationPlot, aes(x = PercentOfReads, y = sampleName, fill = value))
  plot = plot + geom_tile()+ scale_fill_gradient(low="white", high="blue")
  plot = plot + labs(title = "Distribution of known junctions depending on percent of reads used" , x = "Percent of reads", y = "Fraction")

  return (plot)
}


plotAllPlotsInOnePDf <- function(plots, pdfFileName){
  pdf(pdfFileName, onefile = TRUE)
  for (i in seq(length(plots))) {
    do.call("grid.arrange", plots[[i]])
  }
  dev.off()
}






