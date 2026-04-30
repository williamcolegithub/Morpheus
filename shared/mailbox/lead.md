# Mailbox — lead

Append-only. Newest at the bottom. One message per H3 block.

Format:
```
### YYYY-MM-DD HH:MM | from <role> | re: <ticket-id or topic>
<body>
```

The lead is responsible for: ticket assignment, contract arbitration when two roles disagree, unblocking, and final sign-off on each epic.

---

### 2026-04-30 | from lead | re: kickoff
Project scaffolded. `shared/contracts.md` seeded with v1 schemas for the cell/gene caches, embedding artifacts, probe results, and the splits table — read it before producing any artifact, and post here if a shape needs to change. Backlog is in `shared/tickets-status.md`. Suggested first wave to run in parallel: PROBE-2 (data), PROBE-3 (embed), PROBE-6 (data). PROBE-4/5/7 follow data; PROBE-9/10/11 follow embed.

Ground rules:
1. Any contract change → edit `contracts.md` AND post in this mailbox before downstream consumers run.
2. Donor-stratified splits are non-negotiable (PROBE-20 will audit).
3. Controls (PROBE-16 epic) run before any results are written up — random-init + bag-of-genes + selectivity. No exceptions.

### 2026-04-30 | from lead | re: environment + data-source decisions
Folding in upstream context everyone needs to know up front:

- **Hardware target: 24 GB MacBook Air, Apple Silicon.** v1 runs end-to-end here. No assumption of CUDA or >24 GB RAM in any module. Cloud GPU is reserved for an optional final scale-up; not part of v1.
- **v1 scale: 100k cells.** Not 500k. Validate the whole pipeline at 100k first; scaling is cheap once the science works.
- **DoRothEA via `decoupler`** (`get_dorothea` / `get_collectri`) — pure Python, no R dep. CollecTRI is the auto-fallback when post-vocab-filter DoRothEA edges drop below 50k. See contracts §`dorothea_edges.parquet v1`.
- **Vocab filter happens in PROBE-7**, not lazily later. Geneformer's vocab is fixed (~25k human protein-coding); any DoRothEA TF/target outside it is dead weight. Surface the loss now, not during PROBE-14.
- **Frozen embeddings only in v1.** No fine-tuning. Planarian extension and multi-model (scGPT/UCE) comparison are separate epics, separate papers.
- **Streaming + memmap for embedding extraction** is a hard requirement, not a perf optimization. See contracts §`cell_embeddings v1`.
