# Morpheus

## Project overview
Probing project: do single-cell foundation models implicitly encode known regulatory hierarchies? Three layered probes against pretrained Geneformer (frozen embeddings only in v1) — cell-type identity, TF→target edges, and master-regulator/hub identity — each scored against the random-init, bag-of-genes, and label-permutation controls that most probing papers skip. Output is a bioRxiv preprint + open-source repo + Zenodo snapshot.

Ticket breakdown and current status live in `shared/tickets-status.md`. All cross-component data shapes and on-disk cache layouts live in `shared/contracts.md`.

## Stack
- **Language:** Python 3.11+, env pinned via `uv` (or conda) — see PROBE-2.
- **Core libs:** `torch` (MPS backend), `transformers`, `scanpy`, `anndata`, `cellxgene-census`, `decoupler`, `scikit-learn`, `pyarrow`, `numpy`, `scipy`.
- **Model:** Geneformer (HuggingFace `ctheodoris/Geneformer`), frozen weights. scGPT/UCE are out of scope for v1.
- **Data sources:** CELLxGENE Census (human blood/immune, healthy donors); DoRothEA via `decoupler.get_dorothea`; CollecTRI fallback via `decoupler.get_collectri` when post-vocab DoRothEA is too sparse.
- **Probes:** `scikit-learn` linear / logistic-regression probes. No torch on the probe side.
- **Hardware:** 24 GB MacBook Air, Apple Silicon. Frozen-embedding inference must complete here at v1 scale (100k cells). Cloud GPU is reserved for an optional post-v1 scale-up; do not assume CUDA in any module.

## Repo layout (owners in parentheses)
- `/morpheus/**` — Python package (all roles; module ownership documented per-file once created)
  - `morpheus/data/` — Census pull, DoRothEA pull, gene-ID map, splits (data)
  - `morpheus/embed/` — Geneformer load, batched extraction, attention sampling (embed)
  - `morpheus/probes/` — Probe runner CLI, layer 1/2/3 probes, RSA (probe)
  - `morpheus/controls/` — Random-init, bag-of-genes, selectivity, leakage audit (controls)
  - `morpheus/analysis/` — Stats, figures (writeup)
- `/data/` — gitignored; all inputs and processed caches (see contracts §file-system layout)
- `/results/` — gitignored except small parquet result tables and final figures
- `/notebooks/` — exploratory; nothing here is on the critical path
- `/shared/` — Cross-role coordination surface
  - `contracts.md` — Cache shapes, embedding layouts, probe-result schemas, change log (source of truth)
  - `tickets-status.md` — Append-only ticket progress log + current backlog
  - `mailbox/<role>.md` — Inter-role notices, one file per agent role
- `CLAUDE.md` — This file

## Agent roles
Seven roles, one mailbox each under `shared/mailbox/`:
- `lead` — coordination, contract arbitration, epic sign-off
- `data` — env, Census pull, DoRothEA, gene-ID map, splits (PROBE-1 epic)
- `embed` — Geneformer load + extraction + attention (PROBE-8 epic)
- `probe` — Layer 1/2/3 probes + RSA (PROBE-12 epic)
- `controls` — Random-init, bag-of-genes, selectivity, leakage audit (PROBE-16 epic)
- `writeup` — Stats, figures, manuscript, Zenodo, bioRxiv, journal (PROBE-21 + PROBE-27 epics)
- `reviewer` — Adversarial internal-review pass (PROBE-26)

## How to run
The pipeline is staged. Each command depends on the previous stage's outputs being on disk.

### Environment
- Install: `uv sync` (lockfile committed at repo root). Pinned in PROBE-2.
- Tests: `uv run pytest`
- Lint: `uv run ruff check . && uv run mypy morpheus`

### Data stage (PROBE-2 → PROBE-7)
- `uv run python -m morpheus.data.pull_census` → writes `data/processed/cells.h5ad`
- `uv run python -m morpheus.data.pull_dorothea` → writes `data/processed/dorothea_edges.parquet`
- `uv run python -m morpheus.data.build_gene_map` → writes `data/processed/gene_id_map.parquet`, applies vocab filter to dorothea edges, auto-unions CollecTRI when needed
- `uv run python -m morpheus.data.build_splits --k 5` → writes `data/processed/splits.parquet`

### Embedding stage (PROBE-9 → PROBE-11)
- `uv run python -m morpheus.embed.extract --weights ctheodoris/Geneformer --out embeddings/geneformer/`
- `uv run python -m morpheus.embed.extract --weights random_init --out embeddings/geneformer_random_init/`
- `uv run python -m morpheus.embed.attention_sample --n 1000 --out embeddings/geneformer/attention_sample.npz`

### Probe stage (PROBE-13 → PROBE-15, PROBE-17 → PROBE-19)
- `uv run python -m morpheus.probes.run --probe layer1 --model geneformer`
- `uv run python -m morpheus.probes.run --probe layer2 --model geneformer`
- `uv run python -m morpheus.probes.run --probe layer3 --model geneformer`
- Repeat with `--model geneformer_random_init` and `--model bag_of_genes` for controls; add `--permute-labels` for selectivity.

### Analysis stage (PROBE-22 → PROBE-23)
- `uv run python -m morpheus.analysis.stats` → paired tests + effect sizes
- `uv run python -m morpheus.analysis.figures` → PDFs + PNGs into `results/figures/`

## Coordination protocol
1. Shared contracts live in `/shared/contracts.md`. Any role changing a cache shape, embedding layout, or result schema MUST update `contracts.md` AND broadcast via mailbox: `CONTRACT CHANGE: <what> <impact>`. Bumping the version in the H3 heading is part of the change.
2. Before starting any task: read latest `/shared/contracts.md` and `/shared/mailbox/<your-role>.md`.
3. After completing a task: append a row to `/shared/tickets-status.md` with ticket id, owner, status, summary, files touched. Move it out of the backlog list when it transitions to `done`.
4. If blocked on another role's artifact, post in their mailbox referencing the ticket — do not guess at the file shape.
5. If a file or folder referenced in `contracts.md` doesn't exist yet, the first role to need it creates it and announces via mailbox.

## Hard rules (reviewer will enforce)
1. **Donor-stratified splits, no exceptions.** Every cell from a given donor lives in exactly one fold. PROBE-20 audits.
2. **All four controls run before writeup ships.** Random-init + bag-of-genes + selectivity + leakage audit. No figure or claim ships without them.
3. **No fine-tuning in v1.** Frozen embeddings only.
4. **Streaming embedding extraction.** Never densify `cells.h5ad.X` in full; extract via memmap chunks. The 24 GB target is non-negotiable for v1.
5. **Vocab filter happens at PROBE-7.** Geneformer vocab is fixed; surface the DoRothEA edge loss before PROBE-14, not during it. CollecTRI auto-unions when post-filter edges <50k.
6. **Cached artifacts are immutable.** A re-run goes to a new directory (`embeddings/geneformer_v2_*/`), not an overwrite.

## References
- Backlog + status: [`shared/tickets-status.md`](./shared/tickets-status.md)
- Shared contracts: [`shared/contracts.md`](./shared/contracts.md)
- Mailboxes: [`shared/mailbox/`](./shared/mailbox/)
