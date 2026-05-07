# What scale buys for single-cell foundation models: a controlled probe of Geneformer V1, V2-104M, and V2-316M

**Status:** Draft v2 (2026-05-05). Single-author (Independent Researcher). Targeting bioRxiv → eLife / TMLR / PLOS Comp Bio.

> **Major reframe vs. v1:** v1 was structured around "different biological properties scale differently in Geneformer." After computing the within-V2 paired effect-size CIs, the headline finding shifted: **pure parameter scaling within Geneformer-V2 (104M → 316M, 3× parameters at fixed pretraining corpus) produces no measurable improvement on any of our three probes.** The improvements visible across V1→V2 are confounded with corpus expansion (~30M → ~104M cells) and vocabulary change. v2 leads with the cleanest positive finding (hub topology is encoded at the smallest variant and survives every control we threw at it) and is more careful about what the within-V2 nulls actually rule out.

## Abstract

We probe whether Geneformer's frozen single-cell representations encode known regulatory structure, using donor-stratified linear probes on 100,000 healthy human blood/immune cells (CELLxGENE Census) across three pretrained variants: V1-10M (Genecorpus-30M), V2-104M and V2-316M (both pretrained on ~104M cells). Three probes: cell-type identity (L1), TF→target structure (L2, expression-matched negatives), and hub-TF identity (L3, top-decile out-degree in DoRothEA ∪ CollecTRI). All comparisons include a label-permutation selectivity control (Hewitt & Liang 2019) and a bag-of-genes baseline. The within-V2 family forms our only clean parameter-scaling axis.

**Three findings.** (1) **Hub-TF identity is encoded at the smallest variant** (V1-10M, AUC 0.72) and exceeds the bag-of-genes baseline by +0.20 AUC at every scale (CI [+0.13, +0.29] across all three scales) — a gap that survives a random-initialization control at matched embedding dimensionality (random-init AUC 0.42), ruling out the alternative that the signal is geometric room rather than learned structure. (2) **TF→target probe gains are large and selectivity-passing across V1→V2 (+0.054 AUC, p=10⁻¹⁶, n=375 TFs) but do not improve under pure parameter scaling within V2** (V2-316M − V2-104M Δ=−0.003 [−0.014, +0.008], p=0.45). The V1→V2 gain is therefore not separable from pretraining-corpus expansion. (3) **Cell-type identity (L1) shows no measurable improvement under pure parameter scaling** (V2-316M − V2-104M Δ=−0.007 [−0.034, +0.017]); V2-316M reaches a statistical tie with the bag-of-genes baseline (Δ=+0.019 [−0.001, +0.043], p=0.31; we report this as "fails to exceed," not "matches," because we did not perform formal equivalence testing).

The contribution is methodological as much as empirical: scaling claims for single-cell foundation models that compare across model-family transitions confound parameter count with pretraining-corpus and vocabulary changes, and the only clean within-family parameter-scaling axis available to us shows no measurable gain on any probe.

## 1. Methods

### 1.1 Data

We pulled the human blood/immune subset of the CZ CELLxGENE Census (release `2025-11-08`; CZI Cell Science Program et al., 2025), restricted to `tissue_general="blood"`, `disease="normal"`, `is_primary_data=True`, and 10x 3'/5' v1/v2/v3 assays. Donor-stratified sampling (full donors first, trimming only inside the boundary donor) yielded 100,000 cells across 15 donors, 41 datasets, and 71 CL-ontology cell types (Diehl et al. 2016; Osumi-Sutherland et al. 2021). Counts were log1p-normalized to 10,000 per cell (the simple-and-effective transformation per Ahlmann-Eltze & Huber 2023); raw CP10K is reconstructed via `expm1(.X)` for Geneformer's rank-value tokenization.

DoRothEA confidence levels A/B/C (Garcia-Alonso et al. 2019; 1,402 human TFs total) and CollecTRI (Müller-Dott et al. 2023; 1,186 TFs, 43,175 signed interactions) were pulled via `decoupler` (Badia-i-Mompel et al. 2022) and unioned. After resolving endpoints to Ensembl IDs, we obtain 73,269 directed edges across 1,194 unique TFs. We do *not* apply a model-vocabulary filter to this edge file; per-model vocab restriction happens at probe time via each model's gene-index parquet, so V1 and V2 probes draw from the same canonical edge set.

**Donor-stratified 5-fold splits.** Following Zimmerman, Espeland & Langefeld (2021), who show that intra-donor cell correlations exceed inter-donor correlations and treating cells from the same donor as independent inflates type-1 error, we partition cells by donor: greedy bin-packing donors (largest-first) into the smallest-current fold. An assertion in `morpheus/data/build_splits.py` confirms each donor lives in exactly one fold; the assertion passes (15 donors, 5 folds). The same audit was absent from the most-cited cell-type classification benchmark in scRNA-seq (Abdelaal et al. 2019), which uses cell-type-stratified-only CV — meaning published scFM accuracy numbers under that protocol may be inflated by donor leakage. Our protocol is conservative.

### 1.2 Models and embeddings

Three pretrained Geneformer variants:

- **V1-10M** (Theodoris et al. 2023): 10M parameters, 6 layers, 256-dim, 4 heads, vocab 25,424 (gc30M dictionary). Pretrained on Genecorpus-30M (~30M cells).
- **V2-104M** (Geneformer model card, December 2024): 104M parameters, 12 layers, 768-dim, 12 heads, vocab 20,275 (gc104M dictionary). Pretrained on ~104M cells.
- **V2-316M**: 316M parameters, 18 layers, 1152-dim, 18 heads, same vocab and pretraining corpus as V2-104M.

**The V1→V2 transition is not a clean parameter-scaling step.** It changes parameter count *and* pretraining-corpus size *and* vocabulary *and* (per the model cards) training schedule. The only clean parameter-scaling axis available to us is V2-104M → V2-316M (3× parameters at fixed corpus and vocab). We use both axes in this paper and label every comparison explicitly.

**Cell embeddings:** mean of last-hidden-state over non-padding tokens (Reimers & Gurevych 2019). **Gene embeddings:** the input token embedding lookup table itself, *not* contextualized aggregations. We chose the latter because it isolates what pretraining writes into the parameter space without per-batch confounds and because L2/L3 probe questions are properly gene-relational.

Tokenization: Geneformer's published rank-value scheme (median-scaled CP10K, descending sort, top-2048 tokens, dynamic-truncated per batch). Inference: V1 on Apple Silicon MPS (~3h on a 24 GB MacBook Air); V2-104M and V2-316M on AWS EC2 g5.xlarge (NVIDIA A10G 24 GB) in 1h17m and 3h24m respectively, batch size 64 (V2-104M) and 32 (V2-316M).

**Controls.** (a) **Random-initialization** baselines for V1 and V2-104M: identical architecture, untrained weights, same extraction pipeline. (b) **Bag-of-genes**: cell embedding = vocab-restricted log1p-CP10K vector (sparse, 22,029-dim); gene embedding = TruncatedSVD on the gene × cell matrix to 256 components.

### 1.3 Probes

All probes are scikit-learn `LogisticRegression` (`C=1.0`), chosen as the conservative probe class per Hewitt & Liang 2019 and Belinkov 2022. Sensitivity to probe capacity (1-hidden 64-unit MLP) reported in Section 2.5.

- **L1 cell-type identity** (PROBE-13). Multinomial; cell types with ≥200 cells (39 of 71 types, 98,434 cells); donor-stratified 5-fold; macro-F1.
- **L2 TF→target** (PROBE-14). Per-TF binary on gene embeddings, TFs with ≥30 in-vocab targets (375 of 1,194). Negatives: expression-matched (10 quantile bins on log-mean expression in `cells.h5ad`). 80/20 gene split per TF; AUC.
- **L3 hub identity** (PROBE-15). Top-decile out-degree (≥143 targets) defines hubs (120 of 1,194 TFs); 5-fold across TFs; AUC.
- **RSA** (Kriegeskorte et al. 2008): Spearman correlation between gene-embedding pairwise cosine similarity and the binary symmetric DoRothEA TF-TF adjacency, computed over the 1,194-TF subset.

### 1.4 Selectivity and sensitivity

**Selectivity (Hewitt & Liang 2019, PROBE-19):** rerun L2 with TF-target labels permuted within each TF's positive gene set (negative-sampling procedure unchanged). A well-behaved probe should fail.

**Sensitivity sweep:** to forestall reviewer attacks on probe and methodology choices, we rerun the V1, V2-104M, and V2-316M probes under three perturbations: (a) MLP probe (1 hidden layer, 64 units, sklearn `MLPClassifier`); (b) random-pair Layer 2 negatives (vs the matched default); (c) DoRothEA-only and CollecTRI-only restricted regulons (vs the union default).

### 1.5 Statistics

Paired Wilcoxon signed-rank (Wilcoxon 1945) on per-fold or per-TF differences. Bootstrap 95% CIs (Efron 1979; 5,000 resamples). For paired tests with n=5 folds the smallest two-sided Wilcoxon p attainable is 0.0625; we report effect-size CIs as the primary inference.

## 2. Results

### 2.1 Hub-TF identity is encoded at the smallest variant and is real model knowledge (Layer 3)

Geneformer-V1-10M's gene token embeddings — the input lookup table, before any cell flows through the model — separate top-decile-out-degree TFs from peripheral TFs at AUC **0.72** (mean over 5 folds; CI [+0.13, +0.29] vs bag-of-genes). This number does not move across 30× more parameters or ~3.5× more pretraining data:

| Variant | Hub identity AUC (mean) | Δ vs bag-of-genes (CI; p) |
|---|---:|---:|
| V1-10M | 0.72 | +0.219 [+0.132, +0.292]; p=0.063¹ |
| V2-104M | 0.71 | +0.210 [+0.138, +0.273]; p=0.063¹ |
| V2-316M | 0.70 | +0.199 [+0.133, +0.265]; p=0.063¹ |
| Bag-of-genes | 0.50 | — |
| V1-10M random-init | 0.42 | — |
| V2-104M random-init | 0.53 | — |

¹ The 0.063 p-floor is the smallest two-sided Wilcoxon attainable at n=5 paired observations.

The signal is not geometric capacity. **Random-initialized V1-10M at the same 256-dim embedding produces hub-AUC 0.42** — worse than chance, consistent with an uninformative embedding under the imbalanced class distribution. The 0.30 AUC gap between trained and untrained V1 at fixed dimensionality rules out the alternative that hubs are simply linearly separable in any 256-d space.

The signal is not the bag-of-genes ceiling. Bag-of-genes hub-AUC is 0.50 (chance). Geneformer separates hubs from non-hubs at AUC 0.20 above this baseline at every scale.

**Interpretation.** Hub-TF topology is encoded in the input gene embedding lookup table after as little pretraining as Genecorpus-30M with a 10M-parameter model. It is not refined by ~3.5× more pretraining data (V1→V2) nor by 3× more parameters at fixed data (V2-104M→V2-316M). Hub identity is a low-complexity property of transcriptomes — recoverable from co-expression statistics at modest model capacity — that does not require scale to extract. The paper's strongest single finding.

### 2.2 TF→target structure: V1→V2 gains are not separable from corpus expansion (Layer 2)

| Comparison | n | Δ AUC (mean; CI; p) |
|---|---:|---:|
| V2-104M − V1 *(confounded: data + params + vocab)* | 375 | **+0.054** [+0.041, +0.068]; **p=2×10⁻¹⁶** |
| V2-316M − V1 *(confounded)* | 375 | +0.051 [+0.038, +0.066]; p=3×10⁻¹⁵ |
| **V2-316M − V2-104M** *(parameter-only)* | 375 | **−0.003** [−0.014, +0.008]; **p=0.45** |
| V2-316M − bag_of_genes | 375 | +0.167 [+0.151, +0.182]; p=4×10⁻⁵⁰ |
| V2-316M − permuted (selectivity) | 375 | +0.203 [+0.181, +0.214]; p=5×10⁻⁵² |

**The +0.054 AUC gain across V1→V2 is large and significant but is confounded** with the Genecorpus-30M → 104M-cell pretraining-corpus expansion and the gc30M → gc104M vocabulary change. We cannot attribute it to parameter count alone.

**The clean parameter-scaling axis (V2-104M → V2-316M, ~3×) shows no measurable improvement** (Δ=−0.003, CI excludes anything larger than +0.008). At the same time, **V2-316M still substantially exceeds bag-of-genes by +0.17 AUC at p=10⁻⁵⁰** and selectivity passes (permuted-label AUC = 0.498). The L2 signal is real; what is not real is the per-parameter return at this scale.

DoRothEA-only V2-316M reaches median AUC **0.83** (n=214 high-confidence-edge TFs), the highest single value in our experiments and consistent with the interpretation that the V2 family encodes high-confidence regulatory edges substantially better than V1 — but again, the "V2 family" here means *both* more pretraining data and more parameters and a new vocab, not parameters alone.

### 2.3 Cell-type identity: no measurable parameter-scaling improvement; tie with raw expression at the largest scale (Layer 1)

| Variant | Macro-F1 (mean over 5 folds) | Macro-F1 (median) |
|---|---:|---:|
| V1-10M | 0.354 | 0.381 |
| V2-104M | 0.430 | 0.419 |
| V2-316M | 0.424 | 0.451 |
| Bag-of-genes | 0.404 | 0.442 |
| V1-10M random-init | 0.219 | 0.234 |
| V2-104M random-init | 0.191 | 0.183 |

| Comparison | Δ macro-F1 (mean; CI; p) |
|---|---:|
| V2-104M − V1 *(confounded)* | +0.077 [+0.040, +0.117]; p=0.063 |
| V2-316M − V1 *(confounded)* | +0.070 [+0.047, +0.093]; p=0.063 |
| **V2-316M − V2-104M** *(parameter-only)* | **−0.007** [−0.034, +0.017]; **p=1.00** |
| V2-316M − bag_of_genes | +0.019 [−0.001, +0.043]; p=0.31 |

**The L1 within-V2 null is power-limited.** With n=5 paired-fold observations, our 95% bootstrap CI on the V2-104M → V2-316M difference rules out improvements larger than ~+0.02 macro-F1; smaller effects we cannot detect. We report this honestly: parameter scaling within V2 does not produce a *measurable* gain on L1 at our n.

**V2-316M does not exceed the bag-of-genes baseline.** Δ=+0.019 [−0.001, +0.043], p=0.31. We did not perform formal equivalence testing (TOST with a pre-specified equivalence margin), so we report this as "fails to reject the null that they are equal" — not "is equivalent to" or "matches." The honest statement is that **at the largest publicly available Geneformer variant, frozen mean-pool cell embeddings do not measurably exceed log-normalized raw expression for cell-type identity probing**.

The full V1→V2 gain on L1 (Δ=+0.07, n=5, p=0.063) is the largest gain Geneformer shows on cell identity in our data — and it is the one comparison we cannot decompose into parameter-vs-data scaling.

### 2.4 RSA is null at every scale

Spearman correlation between gene-embedding pairwise cosine similarity and DoRothEA's binary symmetric TF-TF adjacency, restricted to the 1,194-TF subset: V1 ρ=−0.028, V2-104M ρ=−0.004, V2-316M ρ=+0.022. Slight positive trend with scale; all near null. The Layer 3 hub probe successfully separates hub from non-hub TFs, but the *fine-grained pairwise* TF-TF regulatory geometry is not captured by raw cosine similarity at any scale. Consistent with linear probing finding "is_hub" via the embedding's first few principal directions while higher-order structure stays unaligned.

### 2.5 Sensitivity: probe class, negative sampling, regulon source

| Sensitivity | Effect on V2-316M numbers |
|---|---|
| MLP probe (vs linear) | L1 +0.04; L2 −0.05; L3 −0.13. MLP overfits on small probe training sets; linear is the conservative correct choice. |
| Random-pair L2 negatives (vs matched) | L2 +0.026 (random easier as expected). V2 > V1 holds under both regimes. |
| DoRothEA-only L2 (vs union) | V2-316M L2 climbs to **0.83** (n=214). Highest single result. |
| CollecTRI-only L2 (vs union) | V2-316M L2 = 0.75. Stable. |
| DoRothEA-only L3 | V2-316M = 0.57 (vs V1 = 0.64). Slight *decrease* with scale on the canonical regulon. |
| CollecTRI-only L3 | V2-316M = 0.74 (vs V1 = 0.73). Stable. |

The L3 saturation finding is regulon-sensitive: on DoRothEA-defined hubs alone, the trend is slightly negative. We report this honestly. The cross-regulon-stable interpretation is that V1 already encodes a superset of what V2 picks up for hub identity; the regulon-restricted views differ in which subset each captures.

## 3. Discussion

### 3.1 What we can conclude

Restricting attention to the only clean parameter-scaling axis available to us (V2-104M → V2-316M, ~3× parameters at fixed pretraining corpus), **none of our three probes show measurable improvement**. Cell-type identity Δ=−0.007 [−0.034, +0.017]; TF→target Δ=−0.003 [−0.014, +0.008]; hub identity Δ=−0.011 [−0.028, +0.007]. The combined V1→V2 transition delivers measurable gains on L2 (large, p=10⁻¹⁶) and on L1 (modest, n=5-limited), but those gains cannot be attributed to parameter count alone — they confound parameter scaling, pretraining-corpus expansion, vocabulary change, and presumably training-schedule differences. Hub topology is encoded at the smallest variant and remains stable across both axes; this is the cleanest positive finding in the paper.

### 3.2 What we cannot conclude

- **We cannot conclude that parameter scaling does not help scFMs in general.** We tested one model family (Geneformer) and one parameter-scaling step (V2-104M → V2-316M, 3×). A null at this specific step does not generalize to scFMs writ large.
- **We cannot decompose the V1→V2 gain into parameter-vs-data contributions.** Doing so would require checkpoints that hold all-but-one factor constant (e.g., a 10M-parameter Geneformer trained on the V2 corpus, or a 104M-parameter Geneformer trained on Genecorpus-30M). Neither exists publicly to our knowledge.
- **We cannot conclude that Geneformer "fails to encode cell-type identity."** The L1 V1→V2 Δ is +0.07, a real if power-limited gain. What we can say is that at the largest publicly available variant, frozen mean-pool embeddings reach a statistical tie with the bag-of-genes baseline.
- **We cannot conclude that hub identity is generally a low-complexity signal across single-cell foundation models.** Our claim is specifically that *Geneformer* encodes hub topology at all tested scales above the bag-of-genes baseline. Generalization across model families requires testing scGPT, UCE, scFoundation, etc. on the same probes.
- **L1 within-V2 null is power-limited.** With n=5 fold-level paired observations, we can rule out improvements larger than ~+0.02 macro-F1 at 95% confidence; smaller effects we cannot detect.
- **The L1 V2-316M vs bag-of-genes comparison is "fails to reject," not "is equivalent to."** Formal equivalence testing was not performed.
- **L3 saturation is regulon-sensitive.** On DoRothEA-only hubs, V2 actually trends slightly below V1.

### 3.3 Connection to prior work

Our finding that scFM embeddings do not decisively beat raw-expression baselines on cell-type identity tasks at scale is consistent with the field's emerging consensus (Boiarsky et al. 2025; Kedzierska et al. 2025; Ahlmann-Eltze, Huber & Anders 2025). What we add: a controlled three-scale parameter sweep with selectivity-passing probes for two distinct levels of regulatory graph structure, and an explicit decomposition that separates parameter-only scaling from corpus expansion. The Theodoris-group scaling paper (Venkatesh et al. 2026, *Nature Computational Science*) reports monotonic accuracy improvements with Geneformer scale on fine-tuned classifiers; our finding that frozen-embedding linear probes show no within-V2 parameter-scaling benefit is not in conflict with that — fine-tuning lets the model rewrite its representation and frozen probing does not — but it implies the gains in Venkatesh et al. are largely due to fine-tuning, not to representations available off-the-shelf.

Civale et al. (2026, *arXiv*) recently showed that intermediate layers can substantially outperform final-layer pooling for scFM trajectory and perturbation tasks. We use final-layer mean-pool throughout and acknowledge that layer-wise probing of Geneformer for our regulatory-hierarchy probes is open.

### 3.4 What this implies for the field

The most actionable observation: **for Geneformer specifically, most of the apparent gain "from scale" attributed in V1-vs-V2 comparisons is gain from corpus expansion and vocabulary update, not from parameter count.** Practitioners using these embeddings as features for downstream models should not expect substantial off-the-shelf improvement from picking V2-316M over V2-104M; the cheaper variant is roughly equivalent on every probe we tested. The value of scaling in Geneformer's published checkpoints, to the extent it exists in our experiments, is concentrated in the V1→V2 transition where corpus and parameters scale jointly.

## 4. Outstanding tickets

- **PROBE-26** internal review pass (after 2–3 day cooling).
- **PROBE-28** repo polish + reproducible README.
- **PROBE-29** Zenodo snapshot + DOI.
- **PROBE-30/32** bioRxiv + journal submission.

## 5. Reproducibility

Pinned `uv.lock`; deterministic seed=0 throughout. Census release 2025-11-08; decoupler 2.1.6; transformers 4.57; PyTorch 2.5.1+cu124 (cloud) / 2.11+mps (laptop). Pipeline reruns end-to-end on a 24 GB MacBook Air for V1 and via two AWS EC2 g5.xlarge sessions (~$5 each) for V2-104M and V2-316M.
