#!/usr/bin/bash

# Annotation description and genomic feature of 5-cytosine methylation detection
# Gene annotation used in this study was in-house curated by another Master’s student and I utilized only annotation of 
# scaffolds 1 through 9 as these contain the most important gene and information that drives phenotypic expression in the table beet plant. 
# The same goes for the repeat annotation which was used in this study. For further DNA methylation level and pattern analysis, 
# I utilized the gene, exon, intron, CDS, transcription start and termination site, UTRs, and the repeat information from the annotation files. 
# Then the number of methylated positions in each of these regions was assessed by overlapping the 5mC position in all three contexts (CG, CHG, CHH) 
# with methylation probability of  ≥ 20%  as truly methylated cytosine called from the different approaches to the regions of the genomic features stated above. 
# The overlapping was done using the line of command line below in bash

bedtools intersect  -f -u 1.0 -a annotation.gff3 -b 5mC_positions.bed | wc -l

# Validating and analyzing DNA methylation patterns in characteristic table beet genes (MYB1, CYP76AD1, and DODA1)
# To determine the precise coding regions of the MYB1, CYP76AD1, and DODA1 genes within the table beet reference genome under investigation, 
# relevant information pertaining to these genes and their sequences was extracted from UniProtKB. The primary accession names assigned to these 
# genes according to UniProtKB were M1ETK3, I3PFJ5.1, and I3PFJ9.1, respectively. Subsequently, the annotated coding regions were acquired in-house 
# through the following procedure: 

gffread annotation.gff3 -g TB-assembly.posscaf.fa -o output.gff -w exon.fa -x CDS.fa -y tr_CDS.fa

# Next, a BLASTp search was conducted using the translated coding sequences by aligning them against the SwissProt database. 
# The following line of code was used for this purpose:

blastp -query tr_CDS.fa -db swissprot/swissprot -out gene_last.txt -outfmt 6 -num_threads 20

# The resulting alignment output was then subjected to sorting based on percentage identity, e-value, and bit score using the following lines of code:

export LC_ALL=C LC_LANG=C;sort -k1,1 -k12,12gr -k11,11g -k3,3gr gene_last.txt > glast_sort.out;

# Collection of the top 1 BLASTp result for each gene

for i in $(cut -f1 glast_sort.out|sort -u); do grep -w -m 1 "$i" gene_last.txt; done > glast_top1.txt.

# Then leveraging the accession name information obtained from UniProtKB, the corresponding gene names in the gene annotation file were identified. 
# These aligned gene names to UniProtKB accession identifiers were used to search the annotation file containing coding sequence (CDS) information. 
# This enabled the identification of the exact CDS sequences and their corresponding information for the target genes. Subsequently, these specific CDS 
# were extracted and subjected to a BLAST alignment against the sequences derived from UniProtKB using default settings, allowing for the determination 
# of sequence similarity and identity values. Rz2 gene sequence information was obtained in-house (from other students).

# To investigate the DNA methylation patterns across all sequence contexts within the identified coding regions of the genes, 
# I employed a comprehensive methodology. Firstly, I collected the reference genome-aligned BAM files derived from two sequencing approaches, 
# along with their associated called methylation probability data, specifically obtained from the Primrose-Pb_CpG_tool and Bsbolt pipeline. 
# The focus was primarily directed toward the CDS regions of genes of interest, determined by their precise start and end coordinates. 

# Next, leveraging the available genome assembly information, I input the obtained methylation data with their corresponding annotation details 
# into the IGV browser to observe and compare the DNA methylation patterns in these genes’ coding exons and regions upstream of the transcription start sites. 
