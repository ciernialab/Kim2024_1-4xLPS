### Script to find TF motif enrichments in promoters of DEGs
# JK September 16, 2019
# JK February 29, 2024

## gist of findMotifs.pl command:
# looking across top genes of interest for shared (enriched) motifs!
# findMotifs.pl will take your top genes and compare them to a background set of genes to make suremthat the enriched motifs identified are not by chance or artifact

### COMMANDS AND PARAMETERS
# findMotifs.pl <inputfile.txt> <promoter set> <output directory> [options]
# bg <background file> (ids to use as background, default: all genes)
# -p <#> (Number of processors to use, default: 1)
# -start <#> (offset from TSS, default=-300) [max=based on Promoter Set]
# -end <#> (offset from TSS, default=50) [max=based on Promoter Set]
## for -start and -end parameters, negative value usually indicates upstream of TSS

# -bg all of the genes in our differential gene expression analysis (16,7888 genes)

# -p 2 max on laptop, -p 12 if working on computing cluster or high-power desktop
# -start -1000 (indicates 1000bp upstream, in 5' --> 3' direction)
# -end 200 (indicates 200bp downstream, in 3' --> 5' direction)

# ultimate output that you're interested in is the knownresults.html file 
# make a sheet that compares motif enrichment results across conditions and brain regions
# R package -- pull out directly from Homer results to compile comparisons


# prepare files to load into findMotifs.pl 
# a. input files of geneIDs (per cluster: 2xLPS, 4xLPS, LPSDownregulated)
# notes for pairwise comparisons:
#### HC 3xLPS only has 1 DEG so not included
#### STR 3xLPS only has 5 DEGs so not included

# b. -bg files of all geneIDs in the list: make a .txt file of all geneIDs, regardless of significance
#### -bg <background file> ids to use as background, default: all genes included in DE analysis

# loading in mm10 promoter set
# mouse promoter set ready in homer suite is mm9, need to load in mm10 version yourself
# create custom promoter set for mm10 (can't use mm9)
loadPromoters.pl -name mm10_promoters -org mouse -id refseq -genome mm10 -tss mm10.tss

### findMotifs.pl command with new promoter set:

# 2xLPS-sensitive cluster
findMotifs.pl 2xLPS_cluster_genes.txt mm10_promoters 2xLPS_cluster_genes_output/ -bg HOMER_Background_AllGenes.txt -start -1000 -end 200 -p 12

# 4xLPS-sensitive cluster
findMotifs.pl 4xLPS_cluster_genes.txt mm10_promoters 4xLPS_cluster_genes_output/ -bg HOMER_Background_AllGenes.txt -start -1000 -end 200 -p 12

# LPSDownregulated cluster
findMotifs.pl LPSDownregulated_cluster_genes.txt mm10_promoters LPSDownregulated_cluster_genes_output/ -bg HOMER_Background_AllGenes.txt -start -1000 -end 200 -p 12

## only 2xLPS-sensitive cluster had significant TF enrichments

#### Generate annotation of homer motifs in individual promoters, just for 2xLPS cluster

# 2xLPS-sensitive cluster
# concatenate all known motifs into one file
cat 2xLPS_cluster_genes_output/knownResults/known*.motif > 2xLPS_cluster_genes_output/knownResults/all.known.motif

# generate annotations
findMotifs.pl 2xLPS_cluster_genes.txt mm10_promoters 2xLPS_cluster_genes_output/ -find 2xLPS_cluster_genes_output/knownResults/all.known.motif > 2xLPS_cluster_genes_output/2xLPS_cluster_genes_outputmotifs.txt -bg HOMER_Background_AllGenes.txt -start -1000 -end 200 -p 12