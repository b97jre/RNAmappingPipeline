# RNA-seq pipeline setup

This is the instruction for how to should do to run this RNAseq pipeline on your data. It requires some manual steps but not to much. The program will run two pipelines. One that will map all the reads to a reference and then create a table with all the count reads in one table.Then it will run a pipeline that will generate different QC measurments and generate a html file and a pdf file of all the quality measurments. 


## Pre-pipeline run
To run the program you will have to do some manual steps before.

1.  Creating two tables that contain information that will be used in the analysis. (This can be done by using the SnakemakeSetup.jar package)   
2.  Fill in the table that contains all the metadata for each sample and the sampleNames. The more things that will be included in this the better. This table is highly dependent on the study. (This step is done manually)   
3.  Create a table that contains all the information regarding where the files are and the paths to the programs. If the program is run on uppmax this step can be ignored. (This step is done manually). 
4.  Setup the path structure so that the snakemake file can be used. (This can be done by using the SnakemakeSetup.jar file).    
5.  Start the pipeline  

### 1. Creating a reads table location and initiate the meta data table  
  The important thing here is to in the end have a tab delimeted table with all reads that you want to include in your final analysis. This can be done manually by just writing it. A more efficient way to start it is to use the java package that you can find here. Download it and install it. If you are working on __UPPMAX__ it is already installed on __/proj/b2013006/sw/modules/SnakemakeSetup.jar__ . To start the search you type the following.

    java -jar /glob/johanr/bin/SnakemakeSetup.jar -p findReads -i <localPath> -suffix <endOfFile> -sep <endOfForwardFile endOfReverseFile> 
     
where: 

*  -p findReads tells the program to find reads and report them.   
*  -i  <folderName> gives the folder where the program will find all the reads that ends with the suffix specified below.  
*  -suffix <endOfreadFileName> will determine will files that will be reported.     
*  -sep <endOfForwardFile endOfReverseFile> will tell you what will be reported.    
     
     
The output in the files will look like this. 

**readTableFile**

|SampleName|forward|reverse|
|----------|-------|-----|
|SampleName1     
|SampleName2     |
|SampleName3     |

**metadataTableFile**

|SampleName     |Name  |Prep_date|Seq_Date|Seq_type|Sep1|Sep2|Sep3|forward|reverse|	
|--|--|--|--|--|--|--|--|--|--|--|--|	
|SampleName1     |SN1 |1|1|Ilumina|1|1|1|/path/to/forwardFile1 |/path/to/reverseFile1|
|SampleName3     |SN2 |1|2|Solid|2|2|1|/path/to/forwardFile2 |/path/to/reverseFile2|
|SampleName3     |SN2 |2|3|PacBio|3|1|2|/path/to/forwardFile3 |/path/to/reverseFile3|

**ParametersFile**

#### Filling the previously generated metadataTableFile. 

This step is important for making sense of the data and should be filled in by the user before the analysis start.The Sep columns should be filled in with meta data that is different between the samples. Things that often is used is __[sick | healthy]__ ,__[male | female]__, __[weight]__,__[Individual]__, __[TissueType]__. This list can be filled in with as many columns as the user wants. It is quite important that it is filled in as much as possible.   This can be done in excel but the end result should be a tab delimeted file (__UNIX__ formated).  





#### Filling the previously generated parameters file. 

This will be generated in the first step aswell. The values that needs to be changed and are specific for each project is the three parameters in the  **#ReferenceInfo** section. The others are set to work on **UPPMAX** otherwise they need to be changed to fit the path.  



|Key|Value|Comment|
|---|-----|-------|
|#Project info|||
|ProjectWD|/path/to/projectWD|Give the full path to were the project should be initiated 
|#Reads info|||
|readsPath|/path/to/reads|End of reads file in below reads directory 
|ReadsSuffix|.fastq.gz|End of reads file in below reads directory 
|readsTable|/path/to/reads_ReadInfo_table.tab.txt|Table with reads info   |
|parametersTable|/path/to/reads_Parameters_table.tab.txt|Table with path info   |
|#ReferenceInfo|||  
|ReferenceFasta|/path/to/fasta/reference|Reference sequence file|
|ReferenceSuffix|.fa|End of reference file (most likely .fa or .fasta)|
|ReferenceGTF|/path/to/fasta/reference|Reference annotation file|
|#ReferenceInfo|||  
|ReferenceFasta|/proj/b2015262/private/references/pig/fastaSequences/Sus_scrofa.Sscrofa10.2.dna_sm.toplevel.fa|Reference sequence file|
|#Program and scriptInfo||Information regarding where all the scripts and programs will be located|  
|snakemakeFile|/proj/b2015262/private/scripts/GBP_scripts/snakemake|Location where all the local snake make files are located|
|Rscripts|/proj/b2015262/private/scripts/GBP_scripts/Snakemake|Location where all the local snake make files are located|
|STAR|star |write out full path if default path is not applicable|
|fastqc|star |write out full path if default path is not applicable|
|multiqc|multiqc |write out full path if default path is not applicable|
|samtools|samtools |write out full path if default path is not applicable|
|RSeQC|/Path/To/RSeQC/scripts/|write out full path if default path is not applicable|



