# Morpheus

Code and analysis pipeline accompanying the preprint
**"What scale buys for single-cell foundation models: a controlled probe of Geneformer V1, V2-104M, and V2-316M"** (Cole, 2026).

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20046043.svg)](https://doi.org/10.5281/zenodo.20046043)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)

The badge above resolves to the **concept DOI** (always points to the latest version). The specific version accompanying the preprint is `v1.0-preprint` ([10.5281/zenodo.20046044](https://doi.org/10.5281/zenodo.20046044)).

The manuscript PDF is in [`manuscript/manuscript.pdf`](./manuscript/manuscript.pdf); the LaTeX source and BibTeX bibliography are alongside it.

## TL;DR

We probe whether Geneformer's frozen single-cell representations encode known regulatory structure, on 100,000 healthy human blood/immune cells (CELLxGENE Census), across three pretrained variants (V1-10M, V2-104M, V2-316M). Three probes (cell-type identity, TF→target, hub-TF identity) and four control families (random-init at matched dimensionality, bag-of-genes, Hewitt-Liang label permutation, donor-leakage audit).

**Three findings:**

1. **Hub-TF identity is encoded at the smallest variant.** V1-10M reaches AUC 0.72, beats bag-of-genes by ≈ +0.20 AUC at every scale, and the gap survives a random-init control at matched dimensionality (256/768/1152-d).
2. **TF→target gains across V1→V2 (+0.054 AUC, p = 10⁻¹⁶) do not survive within V2.** The only clean parameter-scaling axis we have (V2-104M → V2-316M, 3× parameters at fixed corpus) shows no measurable improvement (Δ = −0.003, p = 0.45).
3. **Cell-type identity does not exceed bag-of-genes at any scale.** V2-316M's macro-F1 is 0.42 vs bag-of-genes 0.40 (Δ = +0.019 [−0.001, +0.043], p = 0.31, fails to reject).

The methodological contribution is the parameter-vs-corpus decomposition: scaling claims for single-cell foundation models that compare across model-family transitions confound parameter count with pretraining-corpus and vocabulary changes; the within-family parameter-scaling axis is the one that reveals what scaling parameters alone actually buys.

## Reproducing the results

```bash
# 1. environment (Python 3.11 + pinned deps)
uv sync --extra dev

# 2. data acquisition (CELLxGENE Census + DoRothEA + CollecTRI)
uv run python -m morpheus.data.pull_dorothea
uv run python -m morpheus.data.pull_census
uv run python -m morpheus.data.build_gene_map --variant Geneformer-V2-104M
uv run python -m morpheus.data.build_splits --k 5 --seed 0

# 3. embedding extraction
#    V1 fits on Apple Silicon MPS (~3 hr)
uv run python -m morpheus.embed.extract --weights Geneformer-V1-10M --out geneformer
uv run python -m morpheus.embed.extract --weights random_init --variant Geneformer-V1-10M --out geneformer_random_init
uv run python -m morpheus.embed.bag_of_genes
#    V2-104M and V2-316M need a CUDA GPU; one AWS g5.xlarge session each (~1.3 hr / ~3.4 hr)
uv run python -m morpheus.embed.extract --weights Geneformer-V2-104M --out Geneformer-V2-104M
uv run python -m morpheus.embed.extract --weights Geneformer-V2-316M --out Geneformer-V2-316M
# random-init controls at each scale, plus layer-wise V2-104M:
uv run python -m morpheus.embed.extract --weights random_init --variant Geneformer-V2-104M --out Geneformer-V2-104M_random_init
uv run python -m morpheus.embed.extract --weights random_init --variant Geneformer-V2-316M --out Geneformer-V2-316M_random_init
uv run python -m morpheus.embed.extract --weights Geneformer-V2-104M --save-all-layers --out Geneformer-V2-104M_layerwise

# 4. probes (CPU-bound, runs on the laptop)
for model in geneformer Geneformer-V2-104M Geneformer-V2-316M geneformer_random_init Geneformer-V2-104M_random_init Geneformer-V2-316M_random_init bag_of_genes; do
  for probe in layer1 layer2 layer3 rsa; do
    uv run python -m morpheus.probes.run --probe $probe --model $model
  done
done
# selectivity controls and sensitivity sweep:
bash scripts/run_sensitivity.sh
uv run python -m morpheus.probes.layerwise --model Geneformer-V2-104M_layerwise

# 5. analysis (paired Wilcoxon + bootstrap CIs + figures)
uv run python -m morpheus.analysis.stats
uv run python -m morpheus.analysis.figures_v2
```

Total wall time on a 24 GB MacBook Air + 2 × AWS `g5.xlarge` sessions: ~7 hours of compute, ~$15 of cloud GPU.

The pipeline is fully deterministic (`seed = 0` threaded through every random component); rerunning produces bit-identical parquet outputs.

## Repository layout

- `morpheus/` — Python package
  - `data/` — Census pull, DoRothEA + CollecTRI, gene-ID map, donor-stratified splits
  - `embed/` — Geneformer extraction (cells, gene-token lookup, attention sample, layer-wise variant), bag-of-genes baseline
  - `probes/` — Layer 1/2/3 + RSA + per-layer probe runners
  - `analysis/` — paired stats, figures
- `scripts/` — bash entry points (`run_sensitivity.sh`, etc.)
- `data/` — gitignored; raw downloads, processed caches, embeddings (regenerable)
- `results/` — probe-result parquets (one per probe × model × sensitivity variant), RSA, stats summary, figures (PDF + PNG)
- `manuscript/` — `manuscript.md`, `manuscript.tex`, `manuscript.pdf`, `references.bib`, section drafts under `sections/`, reviewer reports under `reviews/`
- `shared/` — cross-component data contracts and the project-history log

## Data sources

- **CELLxGENE Census**, release `2025-11-08` ([CZI Cell Science Program et al. 2025, *NAR*](https://academic.oup.com/nar/article/53/D1/D886/7912032)) — single-cell RNA-seq corpus
- **DoRothEA** confidence levels A/B/C ([Garcia-Alonso et al. 2019, *Genome Research*](https://genome.cshlp.org/content/29/8/1363)) — TF→target regulons
- **CollecTRI** ([Müller-Dott et al. 2023, *NAR*](https://academic.oup.com/nar/article/51/20/10934/7318114)) — TF→target regulons
- **Geneformer V1-10M** ([Theodoris et al. 2023, *Nature*](https://www.nature.com/articles/s41586-023-06139-9)) — pretrained on Genecorpus-30M
- **Geneformer V2-104M and V2-316M** ([model card](https://huggingface.co/ctheodoris/Geneformer)) — pretrained on the V2 corpus (~104M cells)

## Citation

If you use this code, please cite both the preprint and the Zenodo archive:

```bibtex
@misc{cole2026morpheus,
  author = {Cole, William},
  title = {What scale buys for single-cell foundation models: a controlled probe of Geneformer V1, V2-104M, and V2-316M},
  year = {2026},
  note = {bioRxiv preprint},
  url = {https://github.com/williamcolegithub/Morpheus}
}

@software{cole2026morpheus_zenodo,
  author = {Cole, William},
  title = {williamcolegithub/Morpheus},
  year = {2026},
  publisher = {Zenodo},
  doi = {10.5281/zenodo.20046043},
  url = {https://doi.org/10.5281/zenodo.20046043},
  note = {Concept DOI; always resolves to the latest version. The v1.0-preprint snapshot is at DOI 10.5281/zenodo.20046044.}
}
```

## License

MIT — see [LICENSE](./LICENSE).

The CELLxGENE Census data, DoRothEA / CollecTRI regulons, and Geneformer model weights are distributed under their own respective licenses; this repository does not redistribute any of them, but provides scripts to fetch them from their canonical sources.

## Contact

Issues and pull requests are welcome on [github.com/williamcolegithub/Morpheus](https://github.com/williamcolegithub/Morpheus). For substantive scientific questions, the bioRxiv preprint comments thread is the right place once the preprint is live.
