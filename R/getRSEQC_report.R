#install.packages("optparse")
library('rmarkdown')
library("optparse")

option_list = list(
  make_option(c("-d", "--dir"), type="character", default=NULL, 
              help="Working directory that contains all parameter file", metavar="character"),
  make_option(c("-p", "--parameterFile"), type="character", default="parameters.table.txt", 
              help="The name of the parameter file that contains all the info", metavar="character"),
  make_option(c("-o", "--outputDir"), type="character", default="Working directory", 
              help="output directory were the report will end up ", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$dir)){
  print_help(opt_parser)
  stop("A working directory has to be called (--dir) \n", call.=FALSE)
}

parametersTable = read.table("parameters.table.txt", sep = "\t", header = TRUE)
workingDirectory = as.character(parametersTable$Value[which(parametersTable$Key == "ProjectWD")])
RfunctionsDirectory = as.character(parametersTable$Value[which(parametersTable$Key == "RscriptsDir")])
resultsDir = paste(workingDirectory, "RSeQC", sep = "/")
metaDataTableTxt= "metaData.table.txt"
sampleNameColumn= "sampleName" 
metaDataTableRmd = 'metadata.table.Rmd'


paramsList = list(
  workingDirectory = workingDirectory, 
  RfunctionsDirectory = RfunctionsDirectory,
  resultsDir = resultsDir,
  metaDataTableTxt = metaDataTableTxt,
  sampleNameColumn = sampleNameColumn,
  metaDataTableRmd = metaDataTableRmd)


rmarkdown::render(paste(RfunctionsDirectory,'RSeQCdocument.Rmd', sep = "/"),
                  output_dir = workingDirectory, 
                  params = paramsList
                  )
