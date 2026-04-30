# Tickets Status — Morpheus

Append-only log of ticket progress. One row per status change. Newest at the bottom.

Format: `YYYY-MM-DD | <ticket-id> | <owner> | <status> | <summary> | <files touched>`

Statuses: `started`, `in-progress`, `blocked`, `done`, `reverted`.

Owners (agent roles): `lead`, `data`, `embed`, `probe`, `controls`, `writeup`, `reviewer`.

Ticket IDs follow the `PROBE-N` scheme defined in the project brief. Epics (PROBE-1, PROBE-8, PROBE-12, PROBE-16, PROBE-21, PROBE-27) are tracked here too — mark them `done` only when every child ticket is `done`.

| Date | Ticket | Owner | Status | Summary | Files touched |
| --- | --- | --- | --- | --- | --- |
| 2026-04-30 | — | lead | done | Scaffolded `/shared/` (contracts.md, tickets-status.md, mailbox/) | shared/contracts.md, shared/tickets-status.md, shared/mailbox/*.md |
| 2026-04-30 | — | lead | done | Bootstrap: git init, `pyproject.toml`, `uv sync` (Python 3.11), package skeleton, top-level `CLAUDE.md` + `README.md` | .gitignore, pyproject.toml, README.md, uv.lock, morpheus/__init__.py, morpheus/paths.py |
| 2026-04-30 | PROBE-2 | data | done | env pinned via uv (Python 3.11), all deps installed; MPS available | pyproject.toml, uv.lock |
| 2026-04-30 | PROBE-3 | embed | done | Geneformer weights pulled from HF (V1-10M chosen for 24 GB target: 256-dim, 6L, 4H, vocab 25426); forward pass verified on smoke test | data/raw/geneformer_weights/ |
| 2026-04-30 | PROBE-4 | data | done | Census API verified; tutorial-equivalent query returns AnnData | morpheus/data/pull_census.py |
| 2026-04-30 | PROBE-5 | data | done | 100k blood/normal cells, 15 donors, 71 cell types pulled and log1p-normalized to `cells.h5ad` | morpheus/data/pull_census.py, data/processed/cells.h5ad |
| 2026-04-30 | PROBE-6 | data | done | DoRothEA A/B/C (32k) ∪ CollecTRI (43k) = 75k edges via decoupler 2.x | morpheus/data/pull_dorothea.py, data/processed/dorothea_edges.prevocab.parquet |
| 2026-04-30 | PROBE-7 | data | done | gene_id_map built; vocab filter retains 73,269 edges across 1,194 TFs (97% retention after Geneformer-vocab restriction) | morpheus/data/build_gene_map.py, data/processed/{gene_id_map,dorothea_edges}.parquet |
| 2026-04-30 | PROBE-20 | controls | done | donor-stratified 5-fold splits built; assertion in build_splits.py passes (15 donors, each in exactly one fold) | morpheus/data/build_splits.py, data/processed/splits.parquet |
| 2026-04-30 | PROBE-9 | embed | in-progress | Geneformer cell-embedding extraction running in background on 100k cells; ETA ~2.7h on MPS, batch=16, chunk=2000, dynamic-truncated | morpheus/embed/extract.py |
| 2026-04-30 | PROBE-10 | embed | done | Gene token embeddings (input lookup table) extracted to `embeddings/geneformer/genes.npy` (25426 × 256) at start of extraction | morpheus/embed/extract.py |
| 2026-04-30 | PROBE-15 | probe | done | Layer 3 hub probe: hub AUC 0.67–0.79 across 5 folds (median 0.70) on Geneformer V1-10M gene token embeddings; RSA Spearman ρ=-0.028 vs DoRothEA TF-TF adjacency | morpheus/probes/run.py, results/probes/layer3__geneformer.parquet, results/rsa/rsa__geneformer.parquet |
| 2026-04-30 | PROBE-14 | probe | done | Layer 2 TF→target probe: 375/1194 TFs met threshold; median AUC 0.702, IQR [0.633, 0.765] on Geneformer gene embeddings with expression-matched negatives | morpheus/probes/run.py, results/probes/layer2__geneformer.parquet |
| 2026-04-30 | — | lead | in-progress | Chained pipeline (`scripts/run_rest.sh`) running detached: waits for geneformer extraction, then runs random_init extraction, bag_of_genes, all probes incl. permuted control, stats, figures. Logs at `/tmp/morpheus_run_rest.log` | scripts/run_rest.sh |

---

## Backlog (not yet started)

Source of truth for what's outstanding. When you start a ticket, append a `started` row above and remove it from this list. When it's `done`, append a `done` row above; do not put it back here.

### PROBE-1 — Epic: Environment & data infrastructure
- PROBE-2 [data, 0.5d] Pin Python env (uv/conda); torch, transformers, scanpy, anndata, cellxgene-census, decoupler, scikit-learn. Commit lockfile.
- PROBE-3 [embed, 0.5d] Pull Geneformer weights (HF: `ctheodoris/Geneformer`). Verify forward pass on the example notebook.
- PROBE-4 [data, 0.5d] CELLxGENE Census API working; tutorial query end-to-end.
- PROBE-5 [data, 1d] Pull blood/immune subset, ~100–500k healthy cells; protein-coding only; log-normalize; cache to `data/processed/cells.h5ad`. Blocked by PROBE-4.
- PROBE-6 [data, 0.5d] DoRothEA human regulons (A–C); save edge list to `data/processed/dorothea_edges.parquet`.
- PROBE-7 [data, 0.5d] Reconcile gene IDs (Ensembl ↔ HGNC); produce `gene_id_map.parquet`; ≥95% mapping rate. Blocked by PROBE-5, PROBE-6.

### PROBE-8 — Epic: Embedding extraction
- PROBE-9 [embed, 1d] Batched cell-embedding extraction → `embeddings/geneformer/cells.npy` + `cell_index.parquet`. Blocked by PROBE-3, PROBE-5.
- PROBE-10 [embed, 1d] Gene-level token embeddings → `embeddings/geneformer/genes.npy` + `gene_index.parquet`. Blocked by PROBE-3.
- PROBE-11 [embed, 0.5d] Attention sample (~1k cells) → `attention_sample.npz`. Document sampling.

### PROBE-12 — Epic: Probing experiments
- PROBE-13 [probe, 1.5d] Layer 1 — cell-type identity; donor-stratified 5-fold; macro-F1 + CIs; sanity ≥0.7. Blocked by PROBE-9.
- PROBE-14 [probe, 2d] Layer 2 — TF→target on gene embeddings; expression-matched negatives; AUC distribution figure. Blocked by PROBE-10, PROBE-7.
- PROBE-15 [probe, 2d] Layer 3 — hub identity probe + RSA between gene-embedding RDM and DoRothEA adjacency. Blocked by PROBE-10.

### PROBE-16 — Epic: Controls (run before celebrating)
- PROBE-17 [controls, 1d] Random-init Geneformer baseline; rerun PROBE-13/14/15.
- PROBE-18 [controls, 0.5d] Bag-of-genes baseline (mean log-expression).
- PROBE-19 [controls, 1d] Hewitt–Liang selectivity: permute TF-target labels; AUC should fall to ~0.5.
- PROBE-20 [controls, 0.5d] Donor-leakage audit; explicit assertion in code.

### PROBE-21 — Epic: Analysis & writeup
- PROBE-22 [writeup, 1d] Paired tests across folds: foundation-model vs each baseline. Effect sizes, not just p-values.
- PROBE-23 [writeup, 1d] Figures: (1) probe AUC bars + controls, (2) per-TF AUC distribution, (3) RSA matrix, (4) attention case study (e.g. GATA1 in erythroid).
- PROBE-24 [writeup, 2d] Methods + Results.
- PROBE-25 [writeup, 2d] Intro + Discussion. Frame the negative result honestly if that's what we got.
- PROBE-26 [reviewer, 1d] Internal review pass after 2–3 days of cooling.

### PROBE-27 — Epic: Release & publication
- PROBE-28 [writeup, 1d] Clean GitHub repo; reproducible README; pinned versions; license.
- PROBE-29 [writeup, 0.5d] Zenodo snapshot for citable DOI.
- PROBE-30 [writeup, 0.5d] bioRxiv submission; affiliation "Independent Researcher".
- PROBE-31 [writeup, 0.5d] Bluesky/X thread on preprint day; tag scGPT/Geneformer authors.
- PROBE-32 [writeup, 1d] Journal submission (eLife or TMLR; PLOS Comp Bio backup).
- PROBE-33 [writeup, ongoing] Reviewer responses; budget 1–2 weeks per round.
