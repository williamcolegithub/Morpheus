# Mailbox — embed

Append-only. Newest at the bottom.

Owns: Geneformer weight load, cell-embedding extraction, gene-token-embedding extraction, attention sampling. Produces everything under `data/embeddings/geneformer/`. Also owns the `geneformer_random_init/` control directory (same code path, untrained weights) — coordinate with controls.

Format:
```
### YYYY-MM-DD HH:MM | from <role> | re: <ticket-id or topic>
<body>
```

---

### 2026-04-30 | from lead | re: PROBE-8 epic
PROBE-3 (weight load + forward-pass smoke test) is unblocked — start now. PROBE-10 (gene-token embeddings) only needs the model, not `cells.h5ad`, so you can land it before data finishes PROBE-5. PROBE-9 (cell embeddings) blocks on `cells.h5ad` v1 — watch the data mailbox.

Two things that have bitten prior probing efforts: (1) make sure the rank-value tokenization matches Geneformer's published preprocessing exactly, including the gene-median normalization step — don't reinvent it from the paper; use the official tokenizer in the repo. (2) When you write `cells.npy`, the row order MUST match `cell_index.parquet`. Add an assertion that round-trips a few cell_ids before declaring done. Post here once `META.json` is written — include `embedding_dim` so downstream probes can size their classifiers.

### 2026-04-30 | from lead | re: hardware constraints (MacBook Air, 24 GB)
Target is a 24 GB MacBook Air on Apple Silicon. Frozen-embedding inference only — no fine-tuning in v1. Hard rules:

1. **Streaming is mandatory.** Iterate over cells in chunks of 1k–5k. Open `cells.h5ad` in backed mode (`anndata.read_h5ad(path, backed='r')`), slice sparse rows for the chunk, `.toarray()` only that chunk, tokenize, run forward pass, write embeddings into a pre-allocated `np.memmap` on `cells.npy`. Never hold the full expression matrix densified in RAM. Never hold the full embedding array in RAM during extraction.
2. **Device:** `device = "mps" if torch.backends.mps.is_available() else "cpu"`. Put the helper in one place and import it. Expect 50–200 cells/sec on MPS depending on Geneformer variant. If you hit a known MPS attention bug, fall back to CPU for that op rather than chasing torch versions — note the fallback in `META.json`.
3. **Scale:** 100k cells for v1. Wall-clock budget at this scale: ~15–30 min on MPS. If your run is taking 2+ hours, something is wrong (likely batch size 1, or you're re-tokenizing per cell instead of caching tokens). Post here before doing a 500k run — that's a deliberate decision, possibly on cloud, not a default.
4. **Immutability:** once `cells.npy` is written and consumed by probe, do not overwrite. New extraction → new directory (`embeddings/geneformer_v2_*/`).
5. **Random-init control reuses your code.** Parameterize the extraction on `weights_source` (HF id OR `"random_init"`); controls will call your CLI with the latter to populate `embeddings/geneformer_random_init/`. Don't fork the file.
