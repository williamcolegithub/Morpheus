#!/usr/bin/env bash
# Sensitivity analyses (audit hardening):
#   A. MLP probe across all probe layers
#   B. Random-pair vs expression-matched negatives for Layer 2
#   F. DoRothEA-only vs CollecTRI-only for Layer 2 + Layer 3
set -euo pipefail
cd "$(dirname "$0")/.."

LOG=/tmp/morpheus_sensitivity.log
: > "$LOG"
log() { echo "[$(date -u +%FT%TZ)] $*" | tee -a "$LOG"; }

MODELS=(geneformer Geneformer-V2-104M bag_of_genes)

log "=== A. MLP probe sensitivity (linear is the headline; MLP is the robustness check) ==="
for model in "${MODELS[@]}"; do
  for probe in layer1 layer2 layer3; do
    log "  mlp probe=$probe model=$model"
    uv run python -m morpheus.probes.run --probe "$probe" --model "$model" --probe-type mlp 2>&1 | tee -a "$LOG" | tail -3
  done
done

log "=== B. Random-pair negatives for Layer 2 (vs default matched) ==="
for model in geneformer Geneformer-V2-104M; do
  log "  random-neg model=$model"
  uv run python -m morpheus.probes.run --probe layer2 --model "$model" --neg-mode random 2>&1 | tee -a "$LOG" | tail -3
done

log "=== F. DoRothEA-only and CollecTRI-only for Layer 2 + Layer 3 ==="
for model in geneformer Geneformer-V2-104M; do
  for src in dorothea collectri; do
    for probe in layer2 layer3; do
      log "  reg=$src probe=$probe model=$model"
      uv run python -m morpheus.probes.run --probe "$probe" --model "$model" --regulon-source "$src" 2>&1 | tee -a "$LOG" | tail -3
    done
  done
done

log "=== SENSITIVITY_COMPLETE ==="
