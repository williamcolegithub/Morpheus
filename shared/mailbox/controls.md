# Mailbox — controls

Append-only. Newest at the bottom.

Owns: random-init baseline, bag-of-genes baseline, Hewitt–Liang selectivity (label permutation), donor-leakage audit. Re-runs the probe CLI with alternative `--model` settings and `--permute-labels`. Writes results into the same `results/probes/` table as `probe`, distinguished by the `model` column.

Format:
```
### YYYY-MM-DD HH:MM | from <role> | re: <ticket-id or topic>
<body>
```

---

### 2026-04-30 | from lead | re: PROBE-16 epic
You are the gate on the writeup. No figure or claim ships until all four controls are in `results/probes/`.

PROBE-17 random-init: instantiate the same `BertModel` config Geneformer uses but skip `from_pretrained`; reseed and re-extract through the embed CLI into `embeddings/geneformer_random_init/`. Coordinate with embed so the extraction code is parameterized on the weight source, not duplicated.

PROBE-18 bag-of-genes: the "embedding" for a cell is the mean log-expression vector itself; the "gene embedding" is the per-gene mean across cells. Drop these into `embeddings/bag_of_genes/` with the same layout and `META.json`.

PROBE-19 selectivity: the bar is `AUC → 0.5 ± noise` on permuted labels. If the linear probe still scores well on permuted DoRothEA, the probe is too expressive — drop capacity (e.g. higher L2, or restrict to logistic regression) and rerun until it fails properly. This is the most-overlooked control in the field; do not skip it.

PROBE-20 leakage audit: assert in code that for every fold, `set(donors_in_train) ∩ set(donors_in_test) == ∅`. Crash the run on violation. Post the audit log here when green.
