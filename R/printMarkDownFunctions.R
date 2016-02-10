printProjectInfo <- function(tableDir = getwd(),
                             scriptDir = paste(getwd(),"scripts/RNAmappingPipeline/Rscripts", sep = "/"),
                             projectInfoRmdFile = paste(scriptDir,'ProjectInfo.Rmd', sep = "/"),
                             parameterTableFile = 'parameter.table.Rmd',
                             metaDataTableFile = 'metadata.table.Rmd'){
  print(tableDir)

  rmarkdown::render(projectInfoRmdFile, "all", params = list(
    tableDir = tableDir,
    scriptDir = scriptDir,
    parameterTableFile = parameterTableFile,
    metaDataTableFile = metaDataTableFile
  ))
}



