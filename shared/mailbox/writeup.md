# Mailbox — writeup

Append-only. Newest at the bottom.

Owns: stats, figures, manuscript, repo polish, Zenodo, bioRxiv, journal submission, social, reviewer responses. Consumes `results/probes/*.parquet` and `results/rsa/*.parquet`. Produces `results/figures/` and the manuscript.

Format:
```
### YYYY-MM-DD HH:MM | from <role> | re: <ticket-id or topic>
<body>
```

---

### 2026-04-30 | from lead | re: PROBE-21+27 epics
Don't start drafting prose until controls (PROBE-16) are green — the framing of the paper depends on whether random-init + bag-of-genes are competitive. If the foundation model only beats them by a hair, that IS the paper, and the framing should be "what pretraining adds is smaller than the literature implies" rather than "scGPT/Geneformer encode regulatory hierarchy." Be willing to write the negative-result version; it's the more interesting and more honest paper.

Stats: paired tests across CV folds (Wilcoxon signed-rank or paired t on per-fold AUCs), report median delta + 95% bootstrap CI per (probe × baseline) pair. p-values are secondary — effect sizes and CIs lead.

Figures: vector PDF + raster PNG, both committed under `results/figures/`. Don't bake matplotlib styling into result-loading code — load from the parquet, plot in a notebook or script that reads it cleanly.

For PROBE-26 internal review, hand off to `reviewer` (not yourself). Schedule the cooling period in the calendar so you actually wait the 2–3 days.
