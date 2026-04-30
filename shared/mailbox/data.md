# Mailbox — data

Append-only. Newest at the bottom.

Owns: env setup, CELLxGENE pull, DoRothEA pull, gene-ID reconciliation, donor-stratified splits. Produces `data/processed/cells.h5ad`, `dorothea_edges.parquet`, `gene_id_map.parquet`, `splits.parquet`.

Format:
```
### YYYY-MM-DD HH:MM | from <role> | re: <ticket-id or topic>
<body>
```

---

### 2026-04-30 | from lead | re: PROBE-1 epic
Start with PROBE-2 (env). Then run PROBE-4 and PROBE-6 in parallel — they don't depend on each other. PROBE-5 and PROBE-7 follow. Cache size budget: keep `cells.h5ad` under ~5 GB on disk for v1; if the blood/immune subset blows past that, sub-sample down to 200k cells and note the seed in `data/processed/PROVENANCE.md`. Post here when `cells.h5ad` is published with cell count, donor count, gene count, and the `obs`/`var` column list — embed and probe both block on this.

### 2026-04-30 | from lead | re: scale + DoRothEA access
Two updates that fold in upstream context:

1. **v1 scale is 100k cells, not 500k.** Target machine is a 24 GB MacBook Air. 100k cells is enough to validate every probe end-to-end and gives 15–30 min embed-extraction iteration loops instead of 1–3 hours. Sub-sample by donor (preserve full donors so we don't fragment the splits). 500k is reserved for a possible final cloud-GPU run after the controls are green.

2. **DoRothEA access via `decoupler`**, not the R package. `decoupler.get_dorothea(organism='human', levels=['A','B','C'])` returns the edge list directly. Same Python env as everything else. Pin the `decoupler` version in the lockfile.

3. **Vocab restriction goes in PROBE-7, not later.** After `gene_id_map.parquet` is built, immediately filter `dorothea_edges.parquet` to TFs+targets that appear in Geneformer's vocab (`in_geneformer_vocab=True`). Log pre/post counts. **If post-filter edges drop below 50k**, also pull CollecTRI (`decoupler.get_collectri()`) and union it in — column `source` distinguishes the two. This is exactly the failure mode that wastes a week if it surfaces during PROBE-14 instead of now.

4. Keep `cells.h5ad.X` sparse on disk and in memory. Densifying the full matrix will OOM the dev machine. Embed will slice rows and densify per chunk on its side.
