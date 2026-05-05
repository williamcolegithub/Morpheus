#!/usr/bin/env bash
# Run all probes + sensitivity sweep for Geneformer-V2-316M locally.
set -euo pipefail
cd "$(dirname "$0")/.."
source .venv/bin/activate

LOG=/tmp/v2_316m_probes.log
: > "$LOG"
log() { echo "[$(date -u +%FT%TZ)] $*" | tee -a "$LOG"; }

MODEL=Geneformer-V2-316M

log "=== primary probes ==="
for probe in layer1 layer2 layer3 rsa; do
  log "  probe=$probe"
  uv run python -m morpheus.probes.run --probe "$probe" --model "$MODEL" 2>&1 | tee -a "$LOG" | tail -3
done

log "=== layer2 permuted (selectivity) ==="
uv run python -m morpheus.probes.run --probe layer2 --model "$MODEL" --permute-labels 2>&1 | tee -a "$LOG" | tail -3

log "=== sensitivity: MLP ==="
for probe in layer1 layer2 layer3; do
  log "  mlp probe=$probe"
  uv run python -m morpheus.probes.run --probe "$probe" --model "$MODEL" --probe-type mlp 2>&1 | tee -a "$LOG" | tail -3
done

log "=== sensitivity: random-neg ==="
uv run python -m morpheus.probes.run --probe layer2 --model "$MODEL" --neg-mode random 2>&1 | tee -a "$LOG" | tail -3

log "=== sensitivity: regulon source ==="
for src in dorothea collectri; do
  for probe in layer2 layer3; do
    log "  reg=$src probe=$probe"
    uv run python -m morpheus.probes.run --probe "$probe" --model "$MODEL" --regulon-source "$src" 2>&1 | tee -a "$LOG" | tail -3
  done
done

log "=== regenerate stats + figures ==="
uv run python -m morpheus.analysis.stats 2>&1 | tee -a "$LOG" | tail -10
uv run python -m morpheus.analysis.figures 2>&1 | tee -a "$LOG" | tail -5

log "V2_316M_PROBES_COMPLETE"
