#!/usr/bin/env bash
# Resume from probe stage — extractions already done.
set -euo pipefail
cd "$(dirname "$0")/.."

LOG=/tmp/morpheus_run_probes.log
: > "$LOG"

for model in geneformer geneformer_random_init bag_of_genes; do
  for probe in layer1 layer2 layer3 rsa; do
    echo "[$(date -u +%FT%TZ)] probe=$probe model=$model" | tee -a "$LOG"
    uv run python -m morpheus.probes.run --probe "$probe" --model "$model" 2>&1 | tee -a "$LOG"
  done
done

for probe in layer2; do
  echo "[$(date -u +%FT%TZ)] permuted probe=$probe model=geneformer" | tee -a "$LOG"
  uv run python -m morpheus.probes.run --probe "$probe" --model geneformer --permute-labels 2>&1 | tee -a "$LOG"
done

echo "[$(date -u +%FT%TZ)] running analysis ..." | tee -a "$LOG"
uv run python -m morpheus.analysis.stats 2>&1 | tee -a "$LOG"
uv run python -m morpheus.analysis.figures 2>&1 | tee -a "$LOG"
echo "[$(date -u +%FT%TZ)] PROBES_PIPELINE_COMPLETE" | tee -a "$LOG"
