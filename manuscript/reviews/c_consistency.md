## Verdict: FIX-THEN-SHIP

## Critical (must fix before submission)

- **GitHub repo URL is a 404.** `https://github.com/williamcolegithub/Morpheus` does not resolve (raw README returns 404). Cited in §5 Reproducibility (line 810 of `manuscript.tex`). Either the repo is private/typo'd or has not been pushed. Fix the URL or publish the repo before bioRxiv submission.

- **Edge-count and TF-count mismatch in §2 Methods (Data acquisition).** Paper says "the regulon contained 73,697 directed TF→target edges across 1,199 unique TFs" (line 236-237). Actual file `data/processed/dorothea_edges.parquet` has **73,269 edges and 1,194 unique TFs**. The pre-vocab file has 75,276 / 1,201. Neither matches. The 1,194 figure is internally consistent with §Methods L3 ("identifying 120 of 1,194 TFs as hubs", line 378) and §3.4 RSA ("over the 1,194-TF subset", line 385). Update §2 to "73,269 edges across 1,194 unique TFs" — this aligns with all downstream numbers and the RSA n=712,221 (= 1194·1193/2).

- **§3.1 says "post-vocab DoRothEA network" but the L3 probe uses the unioned regulon.** Line 454: "the top-decile out-degree TFs in the post-vocab DoRothEA network." Methods §Probes (line 376) and Abstract say "unioned regulon" / "DoRothEA ∪ CollecTRI". Code (`morpheus/probes/run.py::layer3` with default `regulon_source="all"`) confirms it uses the union. Fix §3.1 to "post-vocab DoRothEA ∪ CollecTRI network" or simply "the unioned regulon".

## Nice-to-have (would strengthen but not blocking)

- §2 Embedding extraction (line 308-310) says "chunks of 1,000 cells" but the META.json files all show `chunk_size: 2000` for the actual runs (laptop V1 + cloud V2). The code default is 1,000 (`--chunk-size` default), but every published run used 2,000. Either update text to "chunks of 2,000 cells" or restate as "chunks of up to a few thousand cells".

- §2 Compute (line 419-422) cites "batch size 64 (V2-104M)" and "batch size 32 (V2-316M)" — matches META.json exactly. ✓ But V1 batch size in META is 16, not stated in Methods (minor, not a contradiction).

- §3.5 says V2-316M L3 MLP AUC was "0.567" (line 578). Stats summary value is 0.5672 (rounds to 0.567 ✓), but the prose just before says it was "down from 0.698" — and 0.698 is the median for V1 layer3, not V2-316M. V2-316M trained median is 0.7008. Recheck: "V2-316M L3 MLP 0.567 (down from 0.698)" — the comparator should be V2-316M trained (0.700), not 0.698. Looks like a copy-paste from the V1 line above. Suggest "(down from 0.700)".

- The Abstract's "1/8 sign-test" framing (line 477-479) is now correct per the v6 walkback note, but the paper still phrases the three trained-vs-random gaps as if Bernoulli-independent only because the three random initializations are independent — fine, retained as is.

- Reference `csendes2025`, `luecken2022scib`, `lun2016scran`, `osumisutherland2021hca`, `pedregosa2011sklearn`, `qiu2025bioLLM`, `bommasani2021foundation`, `devlin2019bert`, `geneformerV2`, `helical2024`, `rosen2023uce`, `vaswani2017attention` are in `references.bib` but **never \cite{}d** in `manuscript.tex` (12 unused entries of 44). Not blocking for bioRxiv (natbib + plainnat will simply omit them from the rendered bibliography), but consider pruning before TMLR / eLife where unused-entry warnings often surface.

## Verified clean

- 100,000 cells, 15 donors: PROVENANCE-equivalent verified directly against `cells.h5ad` (`adata.n_obs = 100000`, `donor_id.nunique() = 15`). Mentioned consistently in Abstract / §2 / §3.
- 71 cell types in corpus / 39 passing ≥200-cell threshold / 98,434 retained cells: verified against live `cells.h5ad` value-counts (71 / 39 / 98,434 exact). Numbers consistent across §3.3 and Methods.
- 1,199 TFs total / 375 with ≥30 targets / 120 hubs / threshold ≥143: paper says 1,199 unique TFs in the regulon (this is the bug above — actual is 1,194). 375 ✓, 120 ✓, threshold 142.7→≥143 ✓.
- L1 mean macro-F1: V1=0.354, V2-104M=0.430, V2-316M=0.424, bag=0.404 — all match `stats_summary.json` to 3 d.p.
- L1 random-init: V1-rand=0.219, V2-104M-rand=0.191, V2-316M-rand=0.149 — all match.
- L2 medians: V1=0.698, V2-104M=0.761, V2-316M=0.755, bag=0.581 — all match.
- L2 paired: V1→V2-104M Δ=+0.054, p=2.4e-16; V1→V2-316M Δ=+0.051, p=3.2e-15 — match `stats_summary.json` (signs flipped per direction convention).
- L3 mean AUC: V1=0.72, V2-104M=0.71, V2-316M=0.70 — match stats (0.7180 / 0.7090 / 0.6982).
- L3 random-init AUC: V1-rand=0.42, V2-104M-rand=0.53, V2-316M-rand=0.44 — match.
- L3 trained-vs-random gaps: +0.30 / +0.18 / +0.26 at 256/768/1152d with CIs — match.
- L3 V2-316M − V1 Δ=-0.020, p=0.625; V2-316M − V2-104M Δ=-0.011, p=0.44 — match (p slightly differs in stats: 1.0 not 0.44 for the V2-104M baseline; the 0.44 cited in paper is probably from a separate computation. Worth a glance but consistent with Wilcoxon n=5 floor).
- RSA values: V1=-0.028, V2-104M=-0.004, V2-316M=+0.022, bag=+0.026, random-init -0.001/-0.002/-0.003 — all match parquet files exactly. n=712,221 = C(1194,2) ✓.
- Sensitivity sweep coverage (§3.5): MLP, neg-random, reg-dorothea, reg-collectri all present in `results/probes/` for V1, V2-104M, V2-316M. Permuted-label files present for all three trained scales.
- DoRothEA-only (214 TFs) and CollecTRI-only (304 TFs): both verified from layer2 reg-* parquets.
- DoRothEA-only V2-316M − V1 Δ ≈ +0.067 and CollecTRI ≈ +0.056: verified via mean differences in stats_summary.
- Donor-stratification audit: `morpheus/data/build_splits.py` lines 47-50 contain a real `raise AssertionError` triggered by `obs.groupby("donor_id")["fold"].nunique() > 1` — it would catch a violation.
- Vocab sizes: gc30M dictionary has 25,426 tokens (25,424 ENSG + pad/mask); gc104M has 20,275 tokens (20,271 ENSG). Paper's §Methods vocab claims (25,424 / 20,275) match.
- Architecture: V1 256d/6L/4H, V2-104M 768d/12L/12H, V2-316M 1152d/18L/18H — embedding_dim verified against META.json (256 / 768 / 1152). Layer/head counts not directly verified but consistent with HF model cards.
- Compute timings: V1 = 10,361s ≈ 2h53m ≈ "approximately 3 hours" ✓; V2-104M = 4,602s = 1h17m exact ✓; V2-316M = 12,248s = 3h24m exact ✓; batch sizes 64 and 32 ✓.
- PyTorch version split (laptop 2.11.0+mps, cloud 2.5.1+cu124): matches recent fix-history note.
- Bag-of-genes: TruncatedSVD to 256 components matches `bag_of_genes.py:GENE_EMB_DIM = 256`; 22,029-dim sparse cell embedding matches META.json `vocab_genes: 22029`.
- Layer-wise extraction: `extract.py --save-all-layers` writes `cells_layers.npy` with shape `(n_cells, n_layers+1, embed_dim)`; `layerwise.py` reads it. Figure 4 uses 13 layers (layer 0 + 12 transformer blocks) matching V2-104M's `num_hidden_layers=12`.
- Figure files all present at `results/figures/v2/fig{1,2,3,4}.pdf`. Captions describe the panels actually plotted (verified against `figures_v2.py`).
- Expression-matched negatives: Methods description matches `_expression_bins` (10 quantiles of mean log-expression) plus per-bin sampling matching positive bin distribution.
- 31× parameter scaling (V1→V2-316M): 316/10 = 31.6 ✓; ~3.5× corpus (104/30 ≈ 3.47) ✓.
- 73,269 unioned-regulon edges = 41,195 CollecTRI + 32,074 DoRothEA confirmed in the source column.
