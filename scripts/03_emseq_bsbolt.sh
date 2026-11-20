#!/usr/bin/env bash
# =============================================================================
# 03 | EM-seq methylation calling with BSBolt (truth set for benchmarking)
# -----------------------------------------------------------------------------
# Input : paired-end EM-seq reads (Em_seq_R1/R2)
# Output: wgbs_out.CGmap.gz  (per-cytosine methylation in CG, CHG and CHH)
#
# BSBolt (Farrell et al., 2021) uses a forked BWA-MEM for bisulfite/enzymatic
# conversion-aware alignment. After filtering and alignment the mapped reads
# gave ~99x genome coverage.
# =============================================================================
set -euo pipefail

REF="TB-assembly.posscaf.fa"
DB="bsbolt_index"
R1="Em_seq/Em_seq_R1.fastq.gz"
R2="Em_seq/Em_seq_R2.fastq.gz"
THREADS=20

# 1. QC and adapter / quality trimming
fastqc -t "${THREADS}" "${R1}" "${R2}"
trim_galore --paired --illumina -j 30 -e 0.1 -q 20 -O 1 "${R1}" "${R2}"
T1="Em_seq_R1_val_1.fq.gz"
T2="Em_seq_R2_val_2.fq.gz"

# 2. Index the reference
python3 -m bsbolt Index -G "${REF}" -DB "${DB}"

# 3. Align trimmed pairs
python3 -m bsbolt Align -DB "${DB}" -F1 "${T1}" -F2 "${T2}" -O Em_meth -t "${THREADS}"

# 4. Fix mates, coordinate-sort, remove duplicates, index
samtools fixmate -@ "${THREADS}" -p -m Em_meth.bam Em_fixmate.bam
samtools sort    -@ "${THREADS}" -o Em_sorted.bam Em_fixmate.bam
samtools markdup -@ "${THREADS}" -r Em_sorted.bam Em_dup.bam
samtools index   -@ "${THREADS}" Em_dup.bam

# 5. Call methylation in all contexts
python3 -m bsbolt CallMethylation -I Em_dup.bam -O wgbs_out -DB "${DB}" -t "${THREADS}" -verbose > Em_meth_stat.txt
