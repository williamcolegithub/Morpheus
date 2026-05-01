# Shared Contracts — Morpheus (single-cell foundation-model probing)

This file is the source of truth for cross-component data shapes, on-disk cache layouts, and the history of how they change. Every contract change MUST be reflected here AND announced in the mailbox.

Update format: keep each artifact in its own H3 section, with a version in the heading (e.g. `### cells.h5ad  v1`). When you bump a version, leave the old one in place for one migration cycle before deleting.

Scope: this is a single-repo research codebase, not a service mesh. "Contracts" here means the things one agent produces and another agent consumes — cached datasets, embedding matrices, edge lists, probe result tables. If you only consume an artifact in your own module, it does not belong here.

---

## File-system layout (canonical)

All paths are relative to repo root. Agents MUST read/write through these paths; do not invent parallel locations.

```
data/
  raw/                        # immutable downloads, never edited in place
    cellxgene/
    dorothea/
    geneformer_weights/
  processed/
    cells.h5ad                # see "cells.h5ad v1"
    gene_id_map.parquet       # see "gene_id_map.parquet v1"
    dorothea_edges.parquet    # see "dorothea_edges.parquet v1"
  embeddings/
    geneformer/
      cells.npy               # see "cell_embeddings v1"
      cell_index.parquet      # row -> cell_id mapping
      genes.npy               # see "gene_embeddings v1"
      gene_index.parquet      # row -> ensembl_id mapping
      attention_sample.npz    # see "attention_sample v1"
    geneformer_random_init/   # same layout, untrained-weights control
    bag_of_genes/             # mean log-expression "embeddings", same layout
results/
  probes/                     # see "probe_result v1"
  rsa/                        # see "rsa_result v1"
  figures/
```

---

## Data Schemas

### `cells.h5ad`  v1
- Producer: data
- Consumers: embed, controls (bag-of-genes), probe
- Format: AnnData (h5ad), one row per cell.
- `.X`: log1p-normalized expression on protein-coding genes only. Dense or CSR; consumers must handle both.
- `.var.index`: Ensembl gene ID (e.g. `ENSG00000141510`). No version suffix.
- `.var["hgnc_symbol"]`: HGNC symbol or empty string when unmapped.
- `.obs["cell_id"]`: globally unique string, stable across reruns. **Primary key.**
- `.obs["donor_id"]`: required, used for donor-stratified CV splits. NEVER cross-reference donors across splits.
- `.obs["cell_type"]`: CL ontology label (e.g. `CL:0000236`).
- `.obs["cell_type_name"]`: human-readable label.
- `.obs["dataset_id"]`: CELLxGENE Census dataset UUID.
- Filters applied: protein-coding genes only; healthy donors only (no disease tag); blood/immune subset for v1.
- Notes: callers MUST NOT re-normalize. If you need raw counts, pull from `raw/cellxgene/` and document why.

### `gene_id_map.parquet`  v1
- Producer: data
- Consumers: embed, probe
- Columns: `ensembl_id: str` (PK), `hgnc_symbol: str`, `in_geneformer_vocab: bool`, `in_dorothea: bool`.
- One row per protein-coding gene present in `cells.h5ad`. Rows where both `in_geneformer_vocab` and `in_dorothea` are false MAY be omitted.
- AC for upstream producer: `≥95%` of cells.h5ad genes resolve to an HGNC symbol.

### `dorothea_edges.parquet`  v2
- Producer: data
- Consumers: probe (Layer 2, Layer 3), controls (selectivity)
- Columns: `tf_ensembl: str`, `tf_symbol: str`, `target_ensembl: str`, `target_symbol: str`, `confidence: enum("A","B","C")`, `mor: int (-1|0|1)`, `source: enum("dorothea","collectri")`.
- Confidence levels D/E excluded for v1. Bumping to include them is a future change.
- Source: DoRothEA human regulons via `decoupler.op.dorothea(organism='human', levels=['A','B','C'])` ∪ CollecTRI via `decoupler.op.collectri(organism='human')`. Commit the package version + access date to `data/raw/dorothea/PROVENANCE.md`.
- **No model-vocab filter applied here.** Endpoints are resolved to Ensembl IDs via a *union* of every Geneformer `gene_name_id_dict_*.pkl` present on disk (gc30M for V1, gc104M for V2; gc104M wins on collision). Per-model vocab restriction happens at probe time via each model's `gene_index.parquet`, so V1 and V2 probes draw from the same canonical edge set. v2 schema is identical to v1; only the filter contract changed.
- v1 (DEPRECATED): vocab-filtered to a single variant. Replaced because it forced V1 vs V2 probes to read different edge files.

### `gene_id_map_<variant>.parquet`  v1
- Producer: data
- Consumers: provenance only — not consumed by probes/controls.
- One file per `--variant` (e.g. `gene_id_map_Geneformer-V1-10M.parquet`, `gene_id_map_Geneformer-V2-104M.parquet`).
- Columns: `ensembl_id: str` (PK in this file), `hgnc_symbol: str`, `in_geneformer_vocab: bool` (membership in *that variant's* vocab), `in_dorothea: bool`.
- AC: `≥95%` of cells.h5ad genes resolve to an HGNC symbol.

### `cell_embeddings`  v1  (`embeddings/<model>/cells.npy` + `cell_index.parquet`)
- Producer: embed
- Consumers: probe (Layer 1), controls
- `cells.npy`: float32 array, shape `(N_cells, D)`. `D` is model-dependent; record in `embeddings/<model>/META.json`.
- `cell_index.parquet`: columns `row: int64`, `cell_id: str`. `row` is the row index into `cells.npy`. Order MUST match.
- Coverage: every cell in `cells.h5ad` MUST have a row. Consumers may assert this.
- **Streaming requirement**: producer MUST write via `np.memmap` (or chunked HDF5) and append cell-by-cell or in chunks of ≤5k cells. Do NOT hold raw expression data and full embedding array in memory simultaneously — 24 GB MacBook Air is the target environment. Densifying `cells.h5ad.X` to a single dense array is forbidden; slice sparse rows then `.toarray()` per chunk.

### `gene_embeddings`  v1  (`embeddings/<model>/genes.npy` + `gene_index.parquet`)
- Producer: embed
- Consumers: probe (Layer 2, Layer 3), RSA
- `genes.npy`: float32, shape `(N_genes, D)`. `D` matches `cells.npy`.
- `gene_index.parquet`: columns `row: int64`, `ensembl_id: str`, `hgnc_symbol: str`.
- Coverage: every protein-coding gene in the model vocab. Genes outside `cells.h5ad` but in vocab are still included; probe code filters.

### `attention_sample`  v1  (`embeddings/<model>/attention_sample.npz`)
- Producer: embed
- Consumers: writeup (figures), probe (optional)
- `.npz` keys: `cell_ids: (S,) str`, `attention: (S, L, H, T, T) float16`, where `S = sample size (~1000)`, `L = layers`, `H = heads`, `T = tokens (rank-value sequence length)`.
- Sampling strategy: stratified by `cell_type`, documented in `embeddings/<model>/META.json` under `attention_sampling`.
- Float16 to keep file size sane; consumers cast to float32 before any aggregation.

### `probe_result`  v1  (`results/probes/<probe_name>.parquet`)
- Producer: probe, controls
- Consumers: writeup
- One row per (probe target, CV fold).
- Columns: `probe_name: str`, `target: str` (e.g. cell-type label, TF symbol, "is_hub"), `model: str` (`geneformer` | `geneformer_random_init` | `bag_of_genes`), `fold: int`, `metric: enum("auc","macro_f1","accuracy")`, `value: float`, `n_train: int`, `n_test: int`, `seed: int`.
- Selectivity controls (PROBE-19) emit rows with `model = "<model>__permuted"`.
- Consumers compute medians/IQRs/paired tests from this table; never inline into figure code.

### `rsa_result`  v1  (`results/rsa/<name>.parquet`)
- Producer: probe (PROBE-15)
- Consumers: writeup
- Columns: `lhs: str` (e.g. `gene_embeddings_geneformer`), `rhs: str` (e.g. `dorothea_adjacency`), `correlation_metric: enum("spearman","pearson","kendall")`, `value: float`, `n: int`, `seed: int`.

### `splits`  v1  (`data/processed/splits.parquet`)
- Producer: data
- Consumers: probe, controls
- Columns: `cell_id: str`, `fold: int (0..K-1)`.
- Built by donor-stratified K-fold (K=5 for v1) — every cell from a given donor sits in exactly one fold. PROBE-20 audits this.

---

## Code Contracts

### Probe runner CLI  v1
- Producer: probe
- Consumers: controls (re-uses the same CLI for baselines)
- Entry point: `python -m morpheus.probes.run --probe {layer1|layer2|layer3} --model {geneformer|geneformer_random_init|bag_of_genes} [--permute-labels]`
- Reads from canonical paths above. Writes to `results/probes/`.
- MUST be deterministic given a fixed `--seed` (default 0).

### META.json (per embedding directory)  v1
- Required keys: `model`, `weights_source`, `weights_sha256`, `embedding_dim`, `extracted_at` (ISO-8601 UTC), `commit_sha`.
- Optional: `attention_sampling: {strategy, n, seed}`.

---

---

## Operating environment (informational, not a contract)

Primary dev machine: MacBook Air, 24 GB RAM, Apple Silicon. All code MUST run end-to-end here at the v1 scale (100k cells). The 500k–1M-cell scale-up is an optional one-shot cloud GPU run; nothing in the codebase may assume CUDA or >24 GB RAM.

- Default scale for v1: **100k cells** (not 500k). Validate the entire pipeline at this scale before considering a scale-up. Embed inference at 100k on MPS is 15–30 min vs. 1–3 hours at 500k.
- Torch device selection: `"mps" if torch.backends.mps.is_available() else "cpu"`. Centralize this in one helper; do not scatter `torch.device(...)` literals.
- Cache rule: once an embedding artifact is written, treat it as immutable. Re-running PROBE-9/10/11 must produce a new directory (e.g. `embeddings/geneformer_v2_100k/`) rather than overwriting v1.
- Fine-tuning is out of scope for v1 — frozen-embedding probing only. Any fine-tuning (planarian extension, scGPT/UCE comparison) is a separate epic and a separate paper.

---

## Change log

- 2026-04-30 — Initial contract scaffold (lead). All v1 sections seeded; no producers have shipped yet, so shapes may still shift before first real publish. Bump versions once anything is materialized and consumed.
- 2026-04-30 — Added CollecTRI fallback + early vocab restriction to `dorothea_edges.parquet` (lead). Added streaming/memmap requirement to `cell_embeddings`. Added "Operating environment" section pinning the 24 GB MacBook Air target and the 100k-cell v1 scale.
- 2026-05-01 — `dorothea_edges.parquet` bumped to v2: removed the per-variant vocab filter so V1/V2 probes share one canonical edge file. Endpoint resolution now uses a union of every Geneformer name→ensembl dictionary on disk. `gene_id_map.parquet` split into per-variant `gene_id_map_<variant>.parquet` files (provenance only — not consumed downstream).
