RSeQCdir = 'RSeQC'
metaInfo <- read.table("metaData.table.txt", sep = "\t", header = TRUE);

getAllPlots(dir = RSeQCdir,metaInfo)



renderMyDocument <- function(RmdScriptLocation,tableLocation) {
  rmarkdown::render(RmdScriptLocation, params = list(
    region = region,
    start = start
  ))
}
