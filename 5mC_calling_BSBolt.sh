#!/usr/bin/bash

# Methylation Calling and Analysis from the EM_seq dataset
# To analyze the EM-seq data, I utilized the Bisulfite Bolt (Bsbolt) pipeline. 
# This pipeline is a Python-based program and was designed to utilize forked versions of Burrow Wheelers Alignment (BWA) 
# for read alignment, and it engages the following steps of alignment, sorting, and indexing, and the actual methylation calling (Farrell et al., 2021). 
# The full description of all the necessary steps is described in Figure 8 and the lines of codes below;

# Quality check was first performed on the paired-end read using the Fastqc function in bash, then adaptors were removed using trim_galore in the line of code as; 

trim_galore -j 30 -e 0.1 -q 20 -O 1 -a -illumina ../Em_seq/Em_seq_R1 ../Em_seq/Em_seq_R2

# Indexing the reference genome

python3 -m bsbolt Index -G /path/to/ TB-assembly.posscaf.fa -DB /path/to/output

# Read alignment of the obtained trimmed Em_seq read fastq files to the indexed reference

python3 -m bsbolt Align -DB path/to/DB_output -F1 input.fastq -F1 input.fastq -O ./output/Em_meth -t 20

# Note: After filtering and alignment of the reads to the reference genome, the mapped reads accounted for approximately 99x coverage of the genome.

# Pre-processing of the aligned output file using default settings of samtools with options fixmate, 
# sort, markdup, and index to prepare the reads for duplicate removal. The reads were sorted by coordinate, 
# duplicates were removed and reads were indexed respectively. 

# Then the methylation probabilities of cytosine position were called using the code below

python3 -m bsbolt CallMethylation -I Em_dup.bam -O wgbs_out -DB ../path/to/DB_output -t 20 -verbose > Em_meth_stat
