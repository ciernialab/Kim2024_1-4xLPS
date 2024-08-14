#!/bin/bash


#######################################################################################
#RNA QuantSeq Analysis: https://www.lexogen.com/quantseq-data-analysis/

#######################################################################################
#intall packages if don't have them:
#install fastqc with brew or conda
#install multiqc with: pip install multiqc
#From Lexogen: The programs umi2index and collapse_UMI_bam are distributed as binaries compiled on Ubuntu 16.04, CentOS 6.6, CentOS 7, and MACos. The HTS library from samtools (http://www.htslib.org/download/) is required by collapse_UMI_bam. After installation, the HTS library should be added to the standard library search path. Alternatively, environment variable LD_LIBRARY_PATH can be set to point to the location of the HTS library.
#from bbmap:https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/installation-guide/ download bbduk.sh
#######################################################################################
#make folders
mkdir testdata/fastqc_pretrim
mkdir testdata/fastqc_posttrim
mkdir testdata/trimmed_sequences
#run trim.sh srcipt:
#./trim.sh

#it contains this:
for sample in `cat testdata/samples.txt`
do
R1=${sample}
echo ${R1} "pre qc"

#fastqc on pre-trimmed file
fastqc testdata/${R1}*.fastq.gz --outdir testdata/fastqc_pretrim

#index UMIs > add UMI sequence for each read to the read identifier
#The program umi2index pre-processes fastq files by extracting the UMI from a read and storing it as a read index. This step is performed right after demultiplexing to ensure that artificial UMI sequences do not influence subsequent analysis of endogenous read sequences. The program is invoked as follows:
#umi2index <fastq_gz_file_in> <fastq_gz_file_out>
echo ${R1} "umi2index"
./umi2index testdata/${R1}*.fastq.gz testdata/${R1}.UMI.fastq.gz

### remove the adapter contamination, polyA read through, and low quality tails
##https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/
#In ktrim=r mode, once a reference kmer is matched in a read, that kmer and all the bases to the right will be trimmed, leaving only the bases to the left; this is the normal mode for adapter trimming
## quality-trim to Q10 using the Phred algorithm,
echo ${R1} "trimming"

/Users/aciernia/Desktop/programs/bbmap/bbduk.sh in=testdata/${R1}.UMI.fastq.gz out=testdata/${R1}.trim.fastq.gz ref=polyA.fa.gz,truseq_rna.fa.gz k=13 ktrim=r forcetrimleft=11 useshortkmers=t mink=5 qtrim=t trimq=10 minlength=20


#fastqc on post-trimmed file
echo ${R1} "post qc"
fastqc testdata/${R1}.trim.fastq.gz --outdir testdata/fastqc_posttrim


done

#######################################################################################
#collect all qc together with:

multiqc testdata/fastqc_pretrim/ --filename testdata/PreTrim_multiqc_report.html --ignore-samples Undetermined* --interactive

multiqc testdata/fastqc_posttrim/ --filename testdata/PostTrim_multiqc_report.html --ignore-samples Undetermined* --interactive


#move to trimmed sequences folder
mv testdata/*.trim.fastq.gz testdata/trimmed_sequences/
#######################################################################################
