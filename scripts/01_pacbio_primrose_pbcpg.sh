#!/usr/bin/env bash
# =============================================================================
# 01 | PacBio HiFi 5mCpG calling with Primrose + pb-CpG-tools
# -----------------------------------------------------------------------------
# Input : HiFi kinetics BAM (CCS reads with ipd/pw tags) from ICB, BOKU
# Output: mCG_pb1-9.bed  (CpG sites on scaffolds 1-9, coverage >= 4, 5mC >= 20 %)
#
# Steps
#   1. Keep only double-stranded (by-strand "fi" kinetics) reads. At the time of
#      analysis the available callers could not use single-stranded reads.
#   2. Primrose v1.3.0 predicts per-read 5mC probability from polymerase kinetics.
#   3. pbmm2 v1.12.0 aligns the 5mC-tagged reads to the table beet assembly.
#   4. pb-CpG-tools v2.3.1 turns per-read calls into per-site pileup scores.
#   5. Keep scaffolds 1-9 (the nine pseudo-chromosomes), coverage >= 4, score >= 20.
#
# Note: Primrose has since been replaced by PacBio "jasmine"; the commands below
#       reproduce the versions used in this study.
# =============================================================================
set -euo pipefail

REF="TB-assembly.posscaf.fa"                    # table beet reference assembly
RAW_BAM="TableBeet.Kinetics.CCS.hifi.bam"       # HiFi BAM with kinetics
PREFIX="TableBeet.hifi"
THREADS=30

# 1. Double-stranded reads only (-d "fi" keeps reads carrying the forward-strand IPD tag)
samtools view -@ "${THREADS}" -d "fi" -b "${RAW_BAM}" > "${PREFIX}.ds.bam"

# 2. 5mC prediction from kinetics (minimum 3 passes)
primrose --min-passes 3 -j 20 "${PREFIX}.ds.bam" "${PREFIX}.ds.5mC.bam"

# 3. Align to reference; --sort writes a coordinate-sorted, indexed BAM
pbmm2 align --sort -j "${THREADS}" "${REF}" "${PREFIX}.ds.5mC.bam" "${PREFIX}.ds.5mC.aligned.bam"

# 4. Site-level CpG scores (reference-anchored, count-based pileup)
aligned_bam_to_cpg_scores \
    -b "${PREFIX}.ds.5mC.aligned.bam" \
    -f "${REF}" \
    -o pb \
    -p count \
    -m reference \
    -t "${THREADS}"

# 5. Filter. pb-CpG-tools v2 BED columns:
#    1 chrom | 2 start | 3 end | 4 5mC score (0-100) | 5 haplotype | 6 coverage | ...
awk 'BEGIN{OFS="\t"} $1 ~ /^scaffold_[1-9]$/ && $4 >= 20 && $6 >= 4 {print $1, $2, $3, $4, $6}' \
    pb.combined.reference.mincov4.bed > mCG_pb1-9.bed
