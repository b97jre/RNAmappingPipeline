source("scripts/")
RSeQCdir = 'RSeQC'
metaInfo <- read.table("metaData.table.txt", sep = "\t", header = TRUE,comment.char="");

getAllPlots(dir = RSeQCdir,metaInfo)



renderMyDocument <- function(RmdScriptLocation,tableLocation) {
  rmarkdown::render(RmdScriptLocation, params = list(
    region = region,
    start = start
  ))
}
