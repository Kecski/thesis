#!/usr/bin/env bash
set -euo pipefail

# ---- adjust these ----
REMOTE="euler"
REMOTE_DIR="/cluster/home/aimolay/scratch/thesis/tree_simulation/data/"   # trailing slash = copy contents
LOCAL_DIR="./data_cluster/simple_parameter_test_run"
# ----------------------

mkdir -p "$LOCAL_DIR"
rsync -az --partial --progress --human-readable "$REMOTE:$REMOTE_DIR" "$LOCAL_DIR"
echo "Done."