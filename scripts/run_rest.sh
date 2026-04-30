#!/usr/bin/env bash
# Runs the remaining pipeline AFTER the primary geneformer extraction completes.
# Sequence: wait → random_init → bag_of_genes → probes → controls (permuted) → stats → figures.
set -euo pipefail
cd "$(dirname "$0")/.."

LOG=/tmp/morpheus_run_rest.log
echo "[$(date -u +%FT%TZ)] waiting for primary geneformer extraction to finish ..." | tee -a "$LOG"
while pgrep -f "morpheus.embed.extract.*--out geneformer$" > /dev/null || pgrep -f "morpheus.embed.extract --batch-size 16 --chunk-size 2000 --attention-sample 200 --out geneformer" > /dev/null; do
  sleep 30
done
echo "[$(date -u +%FT%TZ)] primary done; starting random_init extraction" | tee -a "$LOG"

if [ ! -f data/embeddings/geneformer/META.json ]; then
  echo "ERROR: geneformer/META.json missing; primary extraction may have failed" | tee -a "$LOG"
  exit 1
fi

uv run python -m morpheus.embed.extract --weights random_init --batch-size 16 --chunk-size 2000 --attention-sample 0 --out geneformer_random_init 2>&1 | tee -a "$LOG"
echo "[$(date -u +%FT%TZ)] random_init done" | tee -a "$LOG"

uv run python -m morpheus.embed.bag_of_genes 2>&1 | tee -a "$LOG"
echo "[$(date -u +%FT%TZ)] bag_of_genes done" | tee -a "$LOG"

for model in geneformer geneformer_random_init bag_of_genes; do
  for probe in layer1 layer2 layer3 rsa; do
    echo "[$(date -u +%FT%TZ)] probe=$probe model=$model" | tee -a "$LOG"
    uv run python -m morpheus.probes.run --probe "$probe" --model "$model" 2>&1 | tee -a "$LOG"
  done
done

# Selectivity (PROBE-19): permuted Layer 2 + Layer 3 on geneformer.
for probe in layer2; do
  echo "[$(date -u +%FT%TZ)] permuted probe=$probe model=geneformer" | tee -a "$LOG"
  uv run python -m morpheus.probes.run --probe "$probe" --model geneformer --permute-labels 2>&1 | tee -a "$LOG"
done

echo "[$(date -u +%FT%TZ)] running analysis ..." | tee -a "$LOG"
uv run python -m morpheus.analysis.stats 2>&1 | tee -a "$LOG"
uv run python -m morpheus.analysis.figures 2>&1 | tee -a "$LOG"
echo "[$(date -u +%FT%TZ)] PIPELINE COMPLETE" | tee -a "$LOG"
