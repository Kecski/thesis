#!/usr/bin/env bash
set -euo pipefail

# ---- adjust these ----
REMOTE="euler"
REMOTE_DIR="/cluster/home/aimolay/scratch/thesis/tree_simulation/data/"   # trailing slash = copy contents
LOCAL_DIR="./data_cluster_multiple_chains/big_rt_in_the_middle_section"
# ----------------------

mkdir -p "$LOCAL_DIR"
rsync -az --partial --progress --human-readable "$REMOTE:$REMOTE_DIR" "$LOCAL_DIR"
echo "Done."