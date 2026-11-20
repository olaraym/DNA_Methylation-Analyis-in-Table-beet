#!/usr/bin/env bash
# =============================================================================
# 04 | 5mC distribution across genomic features + candidate-gene methylation
# -----------------------------------------------------------------------------
# Part A: count methylated cytosines (CG, CHG, CHH; probability >= 20 %) that
#         fall inside genes, exons, introns, CDS, UTRs, TSS/TTS and repeats on
#         scaffolds 1-9 (in-house gene and repeat annotation).
# Part B: locate the coding regions of betalain-pathway genes BvMYB1, CYP76AD1
#         and DODA1 (plus the rhizomania-resistance gene Rz2) and inspect their
#         methylation in IGV.
# =============================================================================
set -euo pipefail

REF="TB-assembly.posscaf.fa"
GFF="annotation.gff3"              # gene + repeat annotation, scaffolds 1-9
SITES="5mC_positions.bed"          # methylated cytosines from scripts 01-03
THREADS=20

# ---------------------------------------------------------------------------
# Part A. Methylated positions per feature type
#   -u     : report each methylated site once if it overlaps the feature
#   -f 1.0 : the whole site must lie inside the feature
# ---------------------------------------------------------------------------
for feature in $(awk '!/^#/ {print $3}' "${GFF}" | sort -u); do
    n=$(awk -v f="${feature}" 'BEGIN{OFS="\t"} $3==f' "${GFF}" \
        | bedtools intersect -u -f 1.0 -a "${SITES}" -b stdin | wc -l)
    printf "%s\t%s\n" "${feature}" "${n}"
done > 5mC_per_feature.tsv

# ---------------------------------------------------------------------------
# Part B. Candidate genes
# UniProtKB accessions: BvMYB1 = M1ETK3, CYP76AD1 = I3PFJ5.1, DODA1 = I3PFJ9.1
# ---------------------------------------------------------------------------

# B1. Extract transcripts, CDS and translated CDS from the annotation
gffread "${GFF}" -g "${REF}" -o output.gff -w exon.fa -x CDS.fa -y tr_CDS.fa

# B2. BLASTp translated CDS against SwissProt
blastp -query tr_CDS.fa -db swissprot/swissprot -out gene_blast.txt -outfmt 6 -num_threads "${THREADS}"

# B3. Sort hits: query, bitscore (desc), e-value (asc), % identity (desc)
export LC_ALL=C
sort -k1,1 -k12,12gr -k11,11g -k3,3gr gene_blast.txt > gene_blast.sorted.txt

# B4. Keep the best hit per query (first line per query in the sorted file)
awk '!seen[$1]++' gene_blast.sorted.txt > gene_blast.top1.txt

# B5. Look up the UniProt accessions of the target genes
grep -E "M1ETK3|I3PFJ5|I3PFJ9" gene_blast.top1.txt > target_genes_hits.txt

# B6. (Manual) Extract the CDS of the matched gene models, BLAST them back
#     against the UniProtKB sequences to confirm identity, then load the
#     aligned BAMs + methylation tracks (Primrose/pb-CpG-tools and BSBolt)
#     with the annotation into IGV to compare methylation over coding exons
#     and the region upstream of the TSS. Rz2 coordinates were provided in-house.
