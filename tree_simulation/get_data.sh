#!/usr/bin/env bash
set -euo pipefail

# ---- adjust these ----
REMOTE="aimolay@euler.ethz.ch"
REMOTE_DIR="/cluster/home/aimolay/scratch/thesis/tree_simulation/data"   # trailing slash = copy contents
LOCAL_DIR="./data_cluster/"
# ----------------------

mkdir -p "$LOCAL_DIR"
rsync -az --partial --progress --human-readable "$REMOTE:$REMOTE_DIR" "$LOCAL_DIR"
echo "Done."