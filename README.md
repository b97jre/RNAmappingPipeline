# Analysing HTS data to identify differentially expressed genes.

#Introduction
This is a description of different part in analysing data that in the end will identify pathways that are important for GBPs

All method will be written so that it can be used in any kind of RNAseq data analysis where there is a genome and annotation available. It will also contain information how to campare different datasets and how to identify other parameters that might affect which genes that are differentially expressed that are not dependent on the question that is asked. 


## Mapping reads to a known reference.

This is the first steps that will analyse the reads and how well they map to the referene and if they differ between the samples.

This is done in some steps. 

1.  Creating a table that contains all the reads that will be used in the analysis. (This can be done by using the SnakemakeSetup.jar package)
2.  Creating a table that contains all the metadata for each sample. The more things that will be included in this the better. This table is highly dependent on the study.(This step is done manually) 
3.  Create a table that contains all the informatoin regarding where the files are and the paths to the programs. If the programs are in the path only the program name can be used as the path. (This step is done manually)
4. Setup the path structure so that the snakemake file can be used. (This can be done by using the SnakemakeSetup.jar file)
5. Run the snakemake file that will run the pipleine that will map the reads and count the number of reads that map to each gene.



