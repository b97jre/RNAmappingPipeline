## Script for analyzing small RNA data with repsecto to biotype and reads length.
## By Jakub Orzechowski Westholm, 2016-02-22
## Requires
## - bowtie
## - samtools
## - perl
## - R
## - bowtie indexes for each biotype. These are currently hard coded:
##   miRNA, tRNA, snRNA, snoRNA, rRNA, protein_coding, lincRNA, pseudogene and repeat


## Input parameters
IN_FASTQGZ=$1
BIOTYPE_DIR=$2
OUT_PLOT_FILE=$3
TMP_DIR=$4


## Hard coded
OUT_TAB_FILE=$(dirname $OUT_PLOT_FILE)/$(basename $OUT_PLOT_FILE .pdf).txt

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" ## directory with this and other helper scripts
PLOT_SCRIPT=${SCRIPT_DIR}/plot_rna_lib_summary.R               ## R script for plotting biotypes and reads lengths 
LABEL=$(basename $IN_FASTQGZ .fastq.gz)                        ## basename of out files
TMP_FASTQ=${TMP_DIR}/tmp_${LABEL}.fastq                        ## temporary fastq file 1
TMP2_FASTQ=${TMP_DIR}/tmp2_${LABEL}.fastq                      ## temporary fastq file 2

## Creat working directory
if [ ! -d "$TMP_DIR" ]; then
  mkdir $TMP_DIR
fi

zcat $IN_FASTQGZ > $TMP_FASTQ


## Look through each biotype and 
## 1) Gap with bowtie (only 1 hit per read), 
## 2) Get hits that map on + strand
## 3) Get read sequences
## 4) Get read lengths
## 5) Make "histogram" table of read lengths    
for BIOTYPE in miRNA tRNA snRNA snoRNA rRNA protein_coding lincRNA pseudogene repeat
do
    echo ${BIOTYPE} 1>&2
    bowtie -M 1 -n 0 -l 17 -S --best --un $TMP2_FASTQ ${BIOTYPE_DIR}/${BIOTYPE} $TMP_FASTQ  | samtools view -S -F 0x14 -| cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort | uniq -c > ${TMP_DIR}/${BIOTYPE}_lengths.txt
    mv $TMP2_FASTQ $TMP_FASTQ
done

## Plot in R: read length distribution and biotype distribution
Rscript $PLOT_SCRIPT $TMP_DIR $LABEL $OUT_PLOT_FILE $OUT_TAB_FILE

## Clean up
rm -r ${TMP_DIR}



