#!/bin/bash

#######################################################################################
#RNA QuantSeq Analysis: https://www.lexogen.com/quantseq-data-analysis/

#######################################################################################
#install if don't have:
#http://bioinf.wehi.edu.au/featureCounts/
#download binary and add to path
#echo 'export PATH=/Users/aciernia/Desktop/programs/subread-1.6.5-MacOSX-x86_64/bin:$PATH' >>~/.bash_profile
#source ~/.bash_profile
#######################################################################################

#run trim.sh srcipt:
#./STAR_alignmm10.sh

mkdir testdata/counts


#it contains this:
for sample in `cat testdata/samples.txt`
do
R1=${sample}


#quant seq kit is FWD stranded

echo ${R1} "count started"
#0 (unstranded), 1 (stranded) and 2 (reversely stranded)
# Number of CPU threads
#set for gene_id > ensembl id mm10
#exon level counts

featureCounts -T 12 -s 1 -t exon -g gene_id -a /Users/aciernia/Desktop/programs/STAR_libs/mm10/annotation/Mus_musculus.GRCm38.92.gtf -o testdata/counts/${R1}.preUMI.counts.txt testdata/star_out/${R1}Aligned.sortedByCoord.out.bam

featureCounts -T 12 -s 1 -t exon -g gene_id -a /Users/aciernia/Desktop/programs/STAR_libs/mm10/annotation/Mus_musculus.GRCm38.92.gtf -o testdata/counts/${R1}.collapsedUMI.counts.txt testdata/star_out/${R1}.collapseUMI.sort.bam

echo ${R1} "count completed"
done

#######################################################################################
#collect all qc together with:
multiqc testdata/counts/*.preUMI.counts.txt.summary --filename testdata/Counts.preUMI --interactive
multiqc testdata/counts/*.collapsedUMI.counts.txt.summary --filename testdata/Counts.collapsedUMI --interactive
#######################################################################################
#######################################################################################
