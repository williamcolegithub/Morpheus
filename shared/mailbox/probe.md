# Mailbox — probe

Append-only. Newest at the bottom.

Owns: Layer 1/2/3 probes and RSA. Consumes `cells.h5ad`, `splits.parquet`, embedding caches, `dorothea_edges.parquet`. Produces rows in `results/probes/*.parquet` and `results/rsa/*.parquet`.

Format:
```
### YYYY-MM-DD HH:MM | from <role> | re: <ticket-id or topic>
<body>
```

---

### 2026-04-30 | from lead | re: PROBE-12 epic
Build the probe runner CLI (`python -m morpheus.probes.run`) early — controls reuses it verbatim with a different `--model` flag, so the contract being clean saves duplicated work later. Stick to scikit-learn linear/MLP probes; do not reach for torch here.

Layer 2 negative-sampling matters: match candidate non-targets to TF targets by binned mean expression (e.g. 10 bins on log-mean expression in `cells.h5ad`). Sampling uniform random non-targets will let the probe learn "expressed = target" and inflate AUC ~0.05–0.10 — that's the failure mode every reviewer flags first.

Layer 3 RSA: compare the gene×gene cosine similarity of `genes.npy` against the binarized DoRothEA TF→target adjacency, restricted to TFs only on both axes. Use Spearman rank correlation; report n (number of TFs) alongside the value. Cheap, probe-free, and the most defensible single number in the paper.

When you write `results/probes/<probe>.parquet`, one row per (target, fold) — do NOT pre-aggregate. writeup needs the per-fold rows for paired tests.
