#!/bin/bash


#######################################################################################
#RNA QuantSeq Analysis: https://www.lexogen.com/quantseq-data-analysis/

#######################################################################################
#intall packages if don't have them:
#install qualimap
#add the location of the Qualimap tool to our PATH variable
#echo 'export PATH=/Users/aciernia/Desktop/programs/qualimap_v2.2.1/:$PATH' >>~/.bash_profile
#source ~/.bash_profile
#qualimap rnaseq --help
#######################################################################################
#make folders
mkdir -p testdata/qualimapQC_preUMI
mkdir -p testdata/qualimapQC_postUMI

#By default, Qualimap will try to open a GUI to run Qualimap, so we need to run the unset DISPLAY command
unset DISPLAY



#it contains this:
for sample in `cat testdata/samples.txt`
do
R1=${sample}
echo ${R1} "qualimapQC RNAseq"

##-outdir: output directory for html report
#-a: Counting algorithm - uniquely-mapped-reads(default) or proportional (each multi-mapped read is weighted according to the number of mapped locations)
#-bam: path/to/bam/file(s)
#-p: Sequencing library protocol - strand-specific-forward, strand-specific-reverse or non-strand-specific (default)
#-gtf: path/to/gtf/file - needs to match the genome build and GTF used in alignment
#--java-mem-size=: set Java memory

#Note that Qualimap must be run with the -outdir option as well as -outformat HTML (which is on by default). MultiQC uses files found within the raw_data_qualimapReport folder (as well as genome_results.txt).

qualimap rnaseq -bam testdata/star_out/${R1}Aligned.sortedByCoord.out.bam -gtf /Users/aciernia/Desktop/programs/STAR_libs/mm10/annotation/Mus_musculus.GRCm38.92.gtf -outdir testdata/qualimapQC_preUMI/${R1}.preUMI --java-mem-size=8G -p strand-specific-forward

qualimap rnaseq -bam testdata/star_out/${R1}.collapseUMI.sort.bam -gtf /Users/aciernia/Desktop/programs/STAR_libs/mm10/annotation/Mus_musculus.GRCm38.92.gtf -outdir testdata/qualimapQC_postUMI/${R1}.postUMI --java-mem-size=8G -p strand-specific-forward


done

#######################################################################################
#collect all qc together with:

multiqc testdata/qualimapQC_preUMI/*.preUMI/ --filename testdata/PreUMI_multiQualimapRNAseq.html --interactive
multiqc testdata/qualimapQC_postUMI/*.postUMI/ --filename testdata/PostUMI_multiQualimapRNAseq.html --interactive
