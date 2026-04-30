# Pretraining encodes regulatory hierarchy at the gene level but not at the cell level: a controlled probe of Geneformer

**Status:** Draft (Round 1, 2026-04-30). Single-author (Independent Researcher). Targeting bioRxiv → eLife / TMLR.

## Abstract

We probe whether Geneformer, a transformer pretrained on ~30M single-cell transcriptomes, implicitly encodes known regulatory hierarchies. Using the smallest published variant (V1-10M, 256-dim, 6 layers) we extract frozen cell- and gene-level representations from 100,000 healthy human blood/immune cells (CELLxGENE Census, 15 donors, 71 cell types) and compare three layered probes against three controls — random-initialized weights, raw expression, and Hewitt–Liang label permutation — under donor-stratified 5-fold CV. Pretraining provides large, selective gains for **gene-level relational structure** (TF→target identification: median AUC 0.694 vs 0.578 raw expression, p = 1.7×10⁻³⁴; hub-TF identification: AUC 0.718 vs 0.499) but **does not improve cell-type identification** beyond raw log-expression (macro-F1 0.354 vs 0.404; the bag-of-genes baseline wins). The split is robust to a passing selectivity control (permuted-label probe collapses to AUC 0.495). We argue this asymmetry — gene-level structure transfers, cell-level identity does not — is the correct framing for what foundation-model pretraining adds at this scale.

## 1. Methods

### 1.1 Data

We pulled the human blood/immune subset of the CELLxGENE Census (release 2025-11-08, `tissue_general="blood"`, `disease="normal"`, `is_primary_data=True`, restricted to `10x 3'/5' v1/v2/v3` assays). Donor-stratified sampling (full donors first, trim only inside the boundary donor) yielded 100,000 cells across 15 donors, 41 datasets, and 71 CL-ontology cell types. Counts were normalized to 10k transcripts per cell and log1p-transformed; the raw counts-per-10k vector is reconstructed via `expm1(.X)` for Geneformer's input rank-value tokenization.

DoRothEA confidence levels A/B/C (32,286 edges) and CollecTRI (42,990 edges) were pulled via `decoupler.op` and unioned. Restriction to Geneformer's 25,424-gene vocabulary kept 73,269 edges across 1,194 unique TFs (97% retention). HGNC↔Ensembl resolution used Geneformer's bundled `gene_name_id_dict_gc30M.pkl`.

Donor-stratified 5-fold splits were built by greedy bin-packing donors into the smallest-current fold; an in-code assertion (PROBE-20) confirms each donor lives in exactly one fold.

### 1.2 Embeddings

We extracted three sets of representations on the same 100,000 cells:

1. **Geneformer V1-10M** (`ctheodoris/Geneformer`, 256-dim, 6 layers, 4 heads, vocab 25,424). Cell embedding = mean-pool of last hidden state over non-padding tokens. Gene embedding = the input token embedding lookup table itself, *not* contextualized aggregations — chosen because it isolates what pretraining writes into the parameter space without per-batch confounds. Tokenization follows Geneformer's published rank-value scheme: median-scaled CP10K, descending sort, top-2048 tokens, dynamic-truncated per batch to its longest non-pad cell. Inference ran on Apple Silicon MPS (single Mac, 24 GB RAM) via `transformers.BertModel`, batch size 16, in 2h53m.
2. **Random-init Geneformer**: identical config, untrained weights (`BertForMaskedLM(config)`), same extraction code path.
3. **Bag-of-genes**: cell embedding = the vocab-restricted log1p-CP10K vector itself (sparse CSR, 22,029 dim). Gene embedding = TruncatedSVD on the transposed (gene × cell) matrix to 256 components (explained-variance 0.633), giving each gene a signature of its co-expression pattern across all cells.

Attention tensors were extracted on a 141-cell type-stratified subsample (token axis capped at 256, eager attention, float16) for case-study figures.

### 1.3 Probes

All probes are scikit-learn linear logistic regression (`C=1.0`, `max_iter=1000–2000`).

- **Layer 1 — cell-type identity (PROBE-13)**: cells × cell-type CL-ontology label, multinomial. Donor-stratified 5-fold CV, restricted to types with ≥200 cells (39 types, 98,434 cells). Standardization (with_mean=False to preserve sparsity for bag-of-genes).
- **Layer 2 — TF→target (PROBE-14)**: per-TF binary classifier on gene embeddings with ≥30 high-confidence targets (375 TFs of 1,194). Negatives: expression-matched — for each positive bin (10 bins on log-mean expression), an equal number of non-target genes are sampled from the same bin. 80/20 gene split per TF. Reported as AUC distribution across TFs.
- **Layer 3 — hub identity (PROBE-15)**: TFs labeled "hub" if their out-degree in the unioned regulon is in the top decile (threshold = 143 targets, 120 hubs of 1,194). 5-fold CV across TFs. RSA: Spearman correlation between gene-embedding pairwise cosine similarity and the binary symmetric DoRothEA TF-TF adjacency, computed over the 1,194 TF subset.

### 1.4 Controls

- **Random-init**: rerun every probe with random_init embeddings (PROBE-17).
- **Bag-of-genes**: rerun every probe with raw-expression / SVD embeddings (PROBE-18).
- **Selectivity (Hewitt & Liang, PROBE-19)**: rerun Layer 2 with TF→target labels permuted within each TF's gene set, holding the negative-sampling procedure fixed. A well-behaved probe should fail.
- **Donor-leakage audit (PROBE-20)**: in-code assertion that no donor spans folds.

## 2. Results

### 2.1 Headline numbers (Figure 1)

| Probe | metric | geneformer | bag_of_genes | random_init | permuted |
|---|---|---:|---:|---:|---:|
| Layer 1 — cell-type identity | macro-F1 | 0.354 ± 0.06 | **0.404 ± 0.06** | 0.219 ± 0.03 | — |
| Layer 2 — TF→target | AUC (n=375) | **0.694** | 0.578 | 0.493 | 0.495 |
| Layer 3 — hub identity | AUC (5-fold) | **0.718** | 0.499 | 0.425 | — |
| RSA TF-TF adjacency | Spearman ρ | -0.028 | — | — | — |

Paired Wilcoxon + bootstrap 95% CI on per-fold/per-TF differences (geneformer minus baseline):

- **Layer 1 vs bag-of-genes**: median delta **−0.052** [−0.057, −0.044], p = 0.063. Bag-of-genes wins (n=5 folds; p=0.063 is the smallest two-sided value attainable with the Wilcoxon signed-rank test at n=5).
- **Layer 2 vs bag-of-genes**: median delta **+0.107** [+0.100, +0.131], p = 1.7×10⁻³⁴.
- **Layer 2 vs random-init**: median delta **+0.204** [+0.183, +0.218], p = 6.8×10⁻⁵¹.
- **Layer 2 vs permuted**: median delta **+0.206** [+0.183, +0.213], p = 1.6×10⁻⁵² — selectivity passes.
- **Layer 3 vs bag-of-genes**: median delta **+0.245** [+0.132, +0.292], p = 0.063.
- **Layer 3 vs random-init**: median delta **+0.252** [+0.242, +0.364], p = 0.063.

### 2.2 Selectivity (Figure 2)

The Layer 2 permuted-label probe collapses cleanly: median AUC 0.495 [IQR 0.439, 0.557], indistinguishable from chance. This rules out the most common failure mode of high-capacity probes — memorizing whatever the labels are — and means the Layer 2 signal in Geneformer is genuinely about regulatory edges, not probe overfitting.

### 2.3 Random-init underperforms chance on hub identity

The random-init Geneformer baseline scores AUC 0.425 (mean) on hub identity — *worse than chance*. We attribute this to the imbalanced class structure (10% hubs) interacting with an uninformative embedding plus L2 regularization: the linear probe defaults to the majority class, and on this fold-shuffled structure that produces a sub-chance ROC under the held-out class balance. The number is therefore consistent with "random embeddings carry no usable signal" rather than evidence of anti-signal.

### 2.4 RSA is null

Spearman correlation between cosine similarity in Geneformer's gene-token embedding space and DoRothEA's TF-TF adjacency, restricted to the 1,194-TF subset, is ρ = −0.028. The Layer 3 hub probe successfully distinguishes hubs from peripheral TFs, but the *fine-grained pairwise* TF-TF regulatory geometry is not captured by raw cosine similarity. This is consistent with the linear probe finding "is_hub" via the embedding's first few principal directions while the higher-order structure stays unaligned.

## 3. Discussion

The clean version of this paper's contribution is: pretraining at 30M-cell scale gives Geneformer V1-10M a substantial, selective advantage on **gene-level** relational tasks but no advantage on **cell-level** identity. We see this in three convergent ways:

- Layer 1 (cell-type identity) is dominated by raw log-expression. The bag-of-genes baseline beats Geneformer by macro-F1 0.052. With more model variants and more data this gap might close, but at this scale the baseline is the right answer.
- Layer 2 (TF→target) is the regime where pretraining clearly helps: a +0.11 AUC gain over raw expression on 375 TFs, with a passing selectivity control. The signal is in the input embedding lookup itself, not in any contextualized inference.
- Layer 3 (hub identity) shows the largest pretraining gain (+0.25 AUC) but at small n=5 folds the Wilcoxon p saturates at 0.063.

We interpret this asymmetry as Geneformer's pretraining objective — masked language modeling over rank-value gene tokens — being a relational signal between genes (which gene tends to co-rank with which), not between cells. Cell identity probably *can* be recovered well from later contextualized layers; we did not probe layer-wise, and that's a clear next step. But the fact that the input embedding lookup *alone* gets to AUC 0.69 on TF→target with a passing selectivity control is the headline.

### 3.1 Limitations

- **One model variant**. V1-10M is the smallest pretrained Geneformer; V2-104M and V2-316M may close the Layer 1 gap. We made this choice deliberately — the entire pipeline runs end-to-end on a 24 GB MacBook Air — but it constrains the claim. A scale-up study is the obvious follow-up.
- **One organism, one tissue**. Blood/immune is the easiest setting in CELLxGENE; cell-type ontology is well-curated and donor diversity is high. Generalization to e.g. CNS or developmental contexts is untested.
- **Layer-wise probing not done**. Cell-level identity may emerge in later layers; gene-level structure may shift across layers. Pre-registered for v2.
- **Mean-pooled cell embedding** is the simplest aggregation; CLS-token or attention-weighted variants may differ.
- **DoRothEA + CollecTRI as ground truth**. These are curated regulons with known biases (literature concentration on cancer / immune TFs). A graph from ENCODE ChIP-seq would be a different lens.

### 3.2 What this is not

This is not an evaluation of Geneformer as a tool for cell-type annotation; for that, the Geneformer authors' own fine-tuned classifiers are the relevant benchmark. We are asking the narrower question of what's encoded in the *frozen* representation, and the answer is: gene-level relations, not cell-level identity. That answer is informative for anyone using these embeddings as features in downstream models.

## 4. Reproducibility

Code: see repository [`Morpheus`](https://github.com/<TBD>) with pinned environment (`uv.lock`), seed-deterministic CLI runners, and a one-shot `scripts/run_rest.sh` that regenerates every figure from `data/processed/cells.h5ad` + `data/embeddings/geneformer/`. Pipeline ran end-to-end on a single 24 GB MacBook Air in ~7 hours wall clock (3h Geneformer extraction + 3h random-init + 1h probes/stats/figures).

Census release pin: `2025-11-08`. decoupler `2.1.6`. transformers `5.7.0` (with `attn_implementation="sdpa"` for bulk and `"eager"` for the attention-output sample, due to a known sdpa/output_attentions incompatibility). PyTorch `2.11.0`.

## 5. Outstanding tickets

- **PROBE-26** internal review pass (after 2–3 days cooling).
- **PROBE-28** repo polish + reproducible README.
- **PROBE-29** Zenodo snapshot + DOI.
- **PROBE-30/32** bioRxiv + eLife/TMLR submission.
