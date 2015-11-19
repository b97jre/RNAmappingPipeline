
library("ggplot2")
library("plyr")

setwd("/Users/johanreimegard/Vetenskap/Data/Wierup/Project_GBP_Pig")

mappingInfo <- read.table("Log.summary.ver2",sep="\t")
metaDataInfo <- read.table("MetaData.tab.table.txt", sep = "\t", header = TRUE);

head(mappingInfo)

colNames <- c("Started_job_on","Started_mapping_on","Finished_on","Mapping_speed_,_Million_of_reads_per_hour","Number_of_input_reads","Average_input_read_length","Uniquely_mapped_reads","Uniquely_mapped_reads_Percent","Average_mapped_length","Number_of_splices_Total","Number_of_splices_Annotated_(sjdb)","Number_of_splices_GT/AG","Number_of_splices_GC/AG","Number_of_splices_AT/AC","Number_of_splices_Non-canonical","MismatchRate","Deletion_rate_per_base","Deletion_average_length","Insertion_rate_per_base","Insertion_average_length","Number_of_reads_mapped_to_multiple_loci","Percent_of_reads_mapped_to_multiple_loci","Number_of_reads_mapped_to_too_many_loci","Percent_of_reads_mapped_to_too_many_loci","Percent_of_reads_unmapped_too_many_mismatches","Percent_of_reads_unmapped_too_short","Percent_of_reads_unmapped_other","fileName")

colnames(mappingInfo)<-colNames

test <-strsplit(as.character(mappingInfo$fileName),split="\\.")
df <- ldply(test)
colnames(df) <- c("Mapping_Name","log")
mappingInfo$Mapping_Name = df$Mapping_Name

Final = merge(mappingInfo,metaDataInfo, by = "Mapping_Name")



ER <- strsplit(as.character(Final$MismatchRate),split="%")
df <- ldply(ER)
colnames(df) <- c("ErrorRate")
Final$MismatchRate = as.numeric(df$ErrorRate)




ER <- strsplit(as.character(Final$Percent_of_reads_unmapped_too_short),split="%")
df <- ldply(ER)
colnames(df) <- c("ErrorRate")
Final$Percent_of_reads_unmapped_too_short = as.numeric(df$ErrorRate)



Final$fractionMappedReads =Final$Uniquely_mapped_reads /  Final$Number_of_input_reads
Final$fractionMappedReadsMultiple =Final$Number_of_reads_mapped_to_multiple_loci /  Final$Number_of_input_reads
Final$fractionMappedReadsMultiple2 =Final$Number_of_reads_mapped_to_too_many_loci /  Final$Number_of_input_reads

Final$totalMapped =(Final$Number_of_reads_mapped_to_too_many_loci +  Final$Uniquely_mapped_reads + Final$Number_of_reads_mapped_to_multiple_loci)/  Final$Number_of_input_reads


ggplot(Final, aes(x=factor(nr),fill=factor(Location),weight=MismatchRate))+
  geom_bar(position="dodge")+ theme(text = element_text(size=23,colour = "black"),axis.text.x = element_text(angle = 90, hjust = 0,vjust = 0.5,colour = "black"))+ylab("Percent error per base")+xlab("Individual")
ggsave(filename="errorRate.pdf")


ggplot(Final, aes(x=factor(nr),fill=factor(Location),weight=Uniquely_mapped_reads))+
  geom_bar(position="dodge")+ theme(text = element_text(size=23,colour = "black"),axis.text.x = element_text(angle = 90, hjust = 0,vjust = 0.5,colour = "black"))+ylab("Uniquely mapped reads")+xlab("Individual")
ggsave(filename="NrOfMappedReads.pdf")


ggplot(Final, aes(x=factor(nr),fill=factor(Location),weight=fractionMappedReads))+
  geom_bar(position="dodge")+ theme(text = element_text(size=23,colour = "black"),axis.text.x = element_text(angle = 90, hjust = 0,vjust = 0.5,colour = "black"))+ylab("Fraction Uniquely Mapped Reads")+xlab("Individual")
ggsave(filename="FractionUniquleyMappedReads.pdf")

ggplot(Final, aes(x=factor(nr),fill=factor(Location),weight=fractionMappedReadsMultiple))+
  geom_bar(position="dodge")+ theme(text = element_text(size=23,colour = "black"),axis.text.x = element_text(angle = 90, hjust = 0,vjust = 0.5,colour = "black"))+ylab("Fraction Mapped Reads to multiple locations (<=3)")+xlab("Individual")
ggsave(filename="FractionMappedReadsMultiple.pdf")

ggplot(Final, aes(x=factor(nr),fill=factor(Location),weight=fractionMappedReadsMultiple2))+
  geom_bar(position="dodge")+ theme(text = element_text(size=23,colour = "black"),axis.text.x = element_text(angle = 90, hjust = 0,vjust = 0.5,colour = "black"))+ylab("Fraction Mapped Reads to multiple locations (>3)")+xlab("Individual")
ggsave(filename="FractionMappedReadsToomany.pdf")

ggplot(Final, aes(x=factor(nr),fill=factor(Location),weight=totalMapped))+
  geom_bar(position="dodge")+ theme(text = element_text(size=23,colour = "black"),axis.text.x = element_text(angle = 90, hjust = 0,vjust = 0.5,colour = "black"))+ylab("Reads mapped (total)")+xlab("Individual")
ggsave(filename="FractionMappedReads.pdf")




ggplot(Final, aes(x=factor(nr),fill=factor(Location),weight=Number_of_input_reads))+
  geom_bar(position="dodge")+ theme(text = element_text(size=23,colour = "black"),axis.text.x = element_text(angle = 90, hjust = 0,vjust = 0.5,colour = "black"))+ylab("Number of input Reads")+xlab("Sample")

ggsave(filename="NrOfReads.pdf")




ggplot(Final, aes(x=factor(nr),fill=factor(Location),weight=Percent_of_reads_unmapped_too_short))+
  geom_bar(position="dodge")+ theme(text = element_text(size=23,colour = "black"),axis.text.x = element_text(angle = 90, hjust = 0,vjust = 0.5,colour = "black"))+ylab("Percent")+xlab("Sample")
ggsave(filename="PercentToShort.pdf")



