#!/usr/bin/bash

# Pacbio sequences and Enzymatic methyl sequences were obtained from ICB, BOKU.
# After these files where obtained, the Pacbio Seq was firts analysed using the procedure below:

# Methylation analysis of the PacBio Sequence 
# Raw Pacbio Hifi reads in BAM format were first filtered to collect only double-stranded reads 
# for downstream processes. Double-starnded reads were collected because at the time of this analysis, 
# the packages available couldn't work with single-stranded reads. This was done using the line of code below; 
   
samtools view -d "fi" -b TableBeet.Kinetics.CCS.hifi.bam > TableBeet.Kinetics.CCS.hifi.ds.bam

# Next, I applied the Primrose software (Version 1.3.0) to the retrieved double-stranded Hifi read data 
to predict the 5mC probabilities using their respective kinetic information with the line of code below;

primrose --min-passes 3 TKC_hifi.ds.bam TKC.hifi.ds.out.bam -j 20

# Next, the obtained HiFi reads containing 5mC tags (provided in an unaligned BAM format) were aligned to a 
tablebeet reference using pbmm2 version 1.12.0 (Pacific Bioscience, 2023) Using the line of code below;

pbmm2 align TB-assembly.posscaf.fa TKC.hifi.ds.out.bam TKC.hifi.ds.out.align.bam -j 30

# Then the indexed bam files were used to obtain Pileup scores for 5mC across CpG sites using the Pb_CpG_tool 
(version 2.3.1) with the line of code below;

aligned_bam_to_cpg_scores -b hifi.sort.bam -f TB-assembly.posscaf.fa -o pb -p count -m reference -t 30

# Then, I collected methylation information from the first 9 scaffolds, which is known to contain the most important 
sequences of the table beet genome, with coverage of ≥ 4 and methylation probability ≥20% using the line of code below;

awk '$1 ~ /^scaffold_[1-9]$/' pb.combined.reference.mincov4.bed | awk '$5>=20 && $6>=4' | awk	'{OFS="\t"} {print $1"\t"$2"\t"$3"\t"$4"\t"$6}' > mCG_pb1-9.bed

