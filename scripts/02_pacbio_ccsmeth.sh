#!/usr/bin/env bash
# =============================================================================
# 02 | PacBio HiFi 5mCpG calling with ccsmeth (deep learning)
# -----------------------------------------------------------------------------
# Input : HiFi reads aligned with pbmm2 (hifi.pbmm2.bam)
# Output: mCG_ccs.bed  (strand-collapsed CpG sites, coverage >= 4)
#
# ccsmeth predicts 5mCpG from CCS kinetics with an attention-based bidirectional
# GRU network (Ni et al., 2023, Nat. Commun. 14:4054).
# =============================================================================
set -euo pipefail

REF="TB-assembly.posscaf.fa"
ALIGNED_BAM="hifi.pbmm2.bam"
MODEL="model_ccsmeth_5mCpG_call_mods_attbigru2s_b21.v1.ckpt"
OUTDIR="ccsmeth_out"
mkdir -p "${OUTDIR}"

# 1. Per-read modification probabilities.
#    The per-read-site TSV is kept (no --rm_per_readsite) because step 2 reads it.
ccsmeth call_mods \
    --input "${ALIGNED_BAM}" \
    --ref "${REF}" \
    --model_file "${MODEL}" \
    --model_type attbigru2s \
    --output "${OUTDIR}/cmod" \
    --threads 50 --threads_call 10 \
    --mode align \
    --modbam modbam \
    --motif CG

# 2. Aggregate to per-site methylation frequency (BED output)
ccsmeth call_freqt \
    --input_path "${OUTDIR}/cmod.per_readsite.tsv" \
    --result_file "${OUTDIR}/CG_ccsmeth.bed" \
    --bed --sort \
    --prob_cf 0.0 \
    --refsites_only \
    --motifs CG --mod_loc 0 \
    --ref "${REF}" \
    --threads 30

# 3. Reformat: chrom, start, end, strand, % methylated, coverage.
#    Minus-strand CpGs are shifted by -1 so both strands collapse onto the
#    same CpG coordinate; then deduplicate and keep coverage >= 4.
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $6, $11, $10}' "${OUTDIR}/CG_ccsmeth.bed" \
  | awk 'BEGIN{OFS="\t"} $4=="+"{print $1, $2, $3, $4, $5, $6}
                         $4=="-"{print $1, $2-1, $3-1, $4, $5, $6}' \
  | sort -u -k1,1 -k2,2n \
  | awk '$6 >= 4' > mCG_ccs.bed
