# Benchmarking PacBio HiFi 5mCpG callers against EM-seq

Two PacBio HiFi methylation pipelines were compared against EM-seq data:

| Pipeline | Principle | Script |
|---|---|---|
| **Primrose + pb-CpG-tools** | Kinetics-based model (PacBio) | [`01_pacbio_primrose_pbcpg.sh`](../scripts/01_pacbio_primrose_pbcpg.sh) |
| **ccsmeth** | Attention-based BiGRU deep learning | [`02_pacbio_ccsmeth.sh`](../scripts/02_pacbio_ccsmeth.sh) |
| **BSBolt (EM-seq)**, the truth set | Enzymatic conversion + short-read alignment, ~99× coverage | [`03_emseq_bsbolt.sh`](../scripts/03_emseq_bsbolt.sh) |

## Truth-set definitions

Only CpG sites with **coverage ≥ 4** on both platforms were compared. Each site was then put into one of three classes:

| Class | 5mC probability | Rationale |
|---|---|---|
| Fully methylated, high confidence | ≥ 99 % | Strict positive set |
| Methylated, relaxed | ≥ 50 % | Accounts for the lower HiFi coverage |
| Unmethylated | ≤ 1 % | Negative set |

## Metrics

With EM-seq as the reference, each HiFi call is a true positive (TP), false positive (FP), true negative (TN) or false negative (FN).

| Metric | Formula | What it tells you |
|---|---|---|
| Precision | $\dfrac{TP}{TP + FP}$ | Share of predicted 5mC sites that are truly methylated |
| Recall | $\dfrac{TP}{TP + FN}$ | Share of truly methylated sites that were recovered |
| F1 score | $\dfrac{2 \times Precision \times Recall}{Precision + Recall}$ | Balance of precision and recall |
| Accuracy | $\dfrac{TP + TN}{TP + TN + FP + FN}$ | Overall correctness |
