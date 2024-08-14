#!/bin/bash


#######################################################################################
#RNA QuantSeq Analysis: https://www.lexogen.com/quantseq-data-analysis/

#######################################################################################
#install STAR
#https://github.com/alexdobin/STAR
#git clone https://github.com/alexdobin/STAR.git
#brew install gcc
#Build STAR:
# note that the path to c++ executable has to be adjusted to its current version > 9
#cd source
#make STARforMacStatic CXX=/usr/local/Cellar/gcc/9.1.0/bin/g++-9
#run with ./STAR -h
#or add to path: export PATH=$PATH:/the_dir_with_STAR
#export PATH=$PATH:/Users/aciernia/Desktop/programs/STAR/source
#echo $PATH
#echo 'export PATH=/Users/aciernia/Desktop/programs/STAR/source:$PATH' >>~/.bash_profile
#source .bash_profile
# build indicies with mm10STARbuild.sh
#######################################################################################

#run trim.sh srcipt:
#./STAR_alignmm10.sh
mkdir testdata/star_out

#it contains this:
for sample in `cat testdata/samples.txt`
do
R1=${sample}

echo ${R1} "unzipping"

#gunzip testdata/trimmed_sequences/${R1}.trim.fastq.gz

echo ${R1} "mapping started"

#allows 10000 files to be open at once:
ulimit -n 10000

STAR --runThreadN 12 --genomeDir /Users/aciernia/Desktop/programs/STAR_libs/mm10/star_indices/ --readFilesIn  testdata/trimmed_sequences/${R1}.trim.fastq --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMattributes NH HI nM AS MD --outSAMtype BAM SortedByCoordinate --outFileNamePrefix testdata/star_out/${R1}


echo ${R1} "mapping completed"
done


#######################################################################################
#mapping qc
#takes Log.final.out files and compiles into HTML

multiqc testdata/star_out/*Log.final.out --filename testdata/STARAlignmentLogs.html --interactive
