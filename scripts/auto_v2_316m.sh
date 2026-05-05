#!/usr/bin/env bash
# Watchdog: poll EC2 for V2-316M extraction completion, then auto-pull embeddings,
# run all probes locally (incl. sensitivity sweep), and regenerate stats/figures.
#
# Run with: nohup bash scripts/auto_v2_316m.sh > /tmp/auto_v2_316m_stdout.log 2>&1 &
# Tail with: tail -f /tmp/auto_v2_316m.log
set -euo pipefail
cd "$(dirname "$0")/.."

EC2_HOST="ec2-user@54.224.93.14"
KEY="$HOME/.ssh/morpheus-key.pem"
REMOTE_LOG="/home/ec2-user/Morpheus/v2_316m.log"
REMOTE_DIR="/home/ec2-user/Morpheus/data/embeddings/Geneformer-V2-316M"
LOCAL_DIR="data/embeddings/Geneformer-V2-316M"

LOG=/tmp/auto_v2_316m.log
: > "$LOG"
log() { echo "[$(date -u +%FT%TZ)] $*" | tee -a "$LOG"; }

ssh_opts=(-o ConnectTimeout=20 -o ServerAliveInterval=30 -o StrictHostKeyChecking=accept-new -i "$KEY")

log "watchdog starting; polling $EC2_HOST every 5 min for V2_316M_REAL_COMPLETE"

# Poll loop. Tolerate transient SSH failures (network blips, IP-rotation, brief overload).
while true; do
  if ssh "${ssh_opts[@]}" "$EC2_HOST" "grep -q V2_316M_REAL_COMPLETE $REMOTE_LOG 2>/dev/null"; then
    log "completion marker detected on remote"
    break
  fi
  sleep 300
done

# Quick sanity: META.json + cells.npy must exist remotely.
log "verifying remote outputs"
if ! ssh "${ssh_opts[@]}" "$EC2_HOST" "test -f $REMOTE_DIR/META.json && test -f $REMOTE_DIR/cells.npy"; then
  log "ERROR: remote outputs missing despite completion marker. Aborting."
  exit 1
fi

# scp the whole directory down.
log "pulling embeddings to laptop"
mkdir -p "$LOCAL_DIR"
scp -i "$KEY" -r "$EC2_HOST:$REMOTE_DIR/" "data/embeddings/" 2>&1 | tail -5 | tee -a "$LOG"

# Show what landed.
log "local files:"
ls -la "$LOCAL_DIR" 2>&1 | tee -a "$LOG"

# Run all primary probes on V2-316M.
log "running primary probes on Geneformer-V2-316M"
source .venv/bin/activate
for probe in layer1 layer2 layer3 rsa; do
  log "  primary probe=$probe"
  uv run python -m morpheus.probes.run --probe "$probe" --model Geneformer-V2-316M 2>&1 | tee -a "$LOG" | tail -3
done

# Selectivity control on layer2.
log "  permuted layer2"
uv run python -m morpheus.probes.run --probe layer2 --model Geneformer-V2-316M --permute-labels 2>&1 | tee -a "$LOG" | tail -3

# Sensitivity sweep.
log "running sensitivity sweep on V2-316M (MLP, random-neg, regulon-source)"
for probe in layer1 layer2 layer3; do
  log "  mlp probe=$probe"
  uv run python -m morpheus.probes.run --probe "$probe" --model Geneformer-V2-316M --probe-type mlp 2>&1 | tee -a "$LOG" | tail -3
done
log "  layer2 random-neg"
uv run python -m morpheus.probes.run --probe layer2 --model Geneformer-V2-316M --neg-mode random 2>&1 | tee -a "$LOG" | tail -3
for src in dorothea collectri; do
  for probe in layer2 layer3; do
    log "  reg=$src probe=$probe"
    uv run python -m morpheus.probes.run --probe "$probe" --model Geneformer-V2-316M --regulon-source "$src" 2>&1 | tee -a "$LOG" | tail -3
  done
done

log "regenerating stats + figures"
uv run python -m morpheus.analysis.stats 2>&1 | tee -a "$LOG" | tail -5
uv run python -m morpheus.analysis.figures 2>&1 | tee -a "$LOG" | tail -5

log "AUTO_V2_316M_COMPLETE"
