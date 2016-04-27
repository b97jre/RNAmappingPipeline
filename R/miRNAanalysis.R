
#gtfFileName = "Mus_musculus.GRCm38.77.fix.gtf"

get_BiotypePerGeneID_from_Ensemble_GTFfile <- function(gtfFileName, id = "gene_id",biotype = "gene_biotype" ){
  gtfFile <- read.table(file = gtfFileName, sep = "\t")
  gtfFileGene  <- gtfFile[gtfFile$V3 == "gene",]
  gtfFileGeneInfo <- read.table(text = as.character(gtfFileGene$V9), sep = ";",colClasses = "character")

  elems <- unlist( strsplit( as.character(gtfFileGeneInfo) , ";" ) )
  gene_ids = elems[grep(id,elems)]
  gene_biotypes = elems[grep(biotype,elems)]

  if(length(gene_ids) != length(gene_ids)){
    print(paste(
      paste("There is not a one to one ratio  between ",id, "and",biotype,"in the gtf file",gtfFileName, sep =" ")
      ,".", sep = ""))
  }

  Geneinfo = data.frame(gene_ids, gene_biotypes)
  Geneinfo$gene_ids = trim(gsub("gene_id ", "", Geneinfo$gene_ids))
  Geneinfo$gene_biotypes = trim(gsub(" gene_biotype ", "", as.character(Geneinfo$gene_biotypes)))


  head(Geneinfo$gene_biotypes)
  gtfFileGene$V9

  levels(as.factor(Geneinfo$gene_biotypes))
  head(gtfFileGene)
  gtfFileGene[grep (pattern = "tRNA", gtfFileGene$V9),]

}
