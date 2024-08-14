#!/bin/bash


#######################################################################################
#RNA QuantSeq Analysis: https://www.lexogen.com/quantseq-data-analysis/

#######################################################################################
#From Lexogen: The programs umi2index and collapse_UMI_bam are distributed as binaries compiled on Ubuntu 16.04, CentOS 6.6, CentOS 7, and MACos. The HTS library from samtools (http://www.htslib.org/download/) is required by collapse_UMI_bam. After installation, the HTS library should be added to the standard library search path. Alternatively, environment variable LD_LIBRARY_PATH can be set to point to the location of the HTS library.

#######################################################################################
#make folders

#it contains this:
for sample in `cat testdata/samples.txt`
do
R1=${sample}
echo ${R1} "index bam"

#Indexed bam files are necessary for many visualization and downstream analysis tools

samtools index testdata/star_out/${R1}Aligned.sortedByCoord.out.bam testdata/star_out/${R1}Aligned.sortedByCoord.out.bai

echo ${R1} "flagstat on unfiltered bam"
#flagstat on reads
samtools flagstat testdata/star_out/${R1}Aligned.sortedByCoord.out.bam > testdata/star_out/${R1}.preUMIcollapse.flagstat.txt

#collapse UMIs
#Program collapse_UMI_bam takes a bam file as input which is generated by aligning a fastq file to a reference genome. This fastq file is derived from <fastq_gz_file_out>, usually by adapter and quality trimming. The output bam file is used for downstream analysis. This version of collapse_UMI_bam has been tested with bam files generated by the star aligner
#collapse_UMI_bam <bam_file_in> <bam_file_out>

echo ${R1} "umi2index"
./collapse_UMI_bam testdata/star_out/${R1}Aligned.sortedByCoord.out.bam testdata/star_out/${R1}.collapseUMI.bam


echo ${R1} "coordinate sort and index"
samtools sort  testdata/star_out/${R1}.collapseUMI.bam -o testdata/star_out/${R1}.collapseUMI.sort.bam

samtools index testdata/star_out/${R1}.collapseUMI.sort.bam testdata/star_out/${R1}.collapseUMI.sort.bai

#remove unmapped reads
#samtools view -h -F 4 -b ${R1}.collapseUMI.sort.bam > ${R1}.collapseUMI.mapped.bam

echo ${R1} "flagstat on filtered bam"
#flagstat on reads
samtools flagstat testdata/star_out/${R1}.collapseUMI.sort.bam > testdata/star_out/${R1}.collapseUMI.flagstat.txt


done

#######################################################################################
#collect all qc flagstats together with:

multiqc testdata/star_out/*.flagstat.txt --filename testdata/STARAlignment.html --interactive

#######################################################################################
