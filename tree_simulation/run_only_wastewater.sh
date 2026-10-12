#!/bin/zsh

TREE_NAME="sim"
MAX_TIME=0.4
N_SIMS=1



WW_OUTPUT_DIR=./data/wastewater
TREE_FILE=./data/remaster_trees/logs/$TREE_NAME.full.trees

N_DAYS=$(echo "$MAX_TIME * 365.25" | bc -l)
N_DAYS=${N_DAYS%.*}
WW_DATA_DIR="data"
SEED=42
SEQUENCING_DELAY=14

# Shedding Parameters
TOTAL_LOAD_PER_INFECTION=1e9
SHED_MAX_DAYS=21
SHED_MEAN=4
SHED_SD=3
TRANSPORT_MEAN=0.5
TRANSPORT_SD=0.5
FLOW_MEAN=1e8
MEAS_CV=0.3
LOD=1

uv run wastewater_data.py \
    --tree $TREE_FILE \
    --n-days $N_DAYS \
    --output_dir $WW_OUTPUT_DIR \
    --total-load-per-infection $TOTAL_LOAD_PER_INFECTION \
    --shed-max-days $SHED_MAX_DAYS \
    --shed-mean $SHED_MEAN \
    --shed-sd $SHED_SD \
    --transport-mean $TRANSPORT_MEAN \
    --transport-sd $TRANSPORT_SD \
    --flow-mean $FLOW_MEAN \
    --meas-cv $MEAS_CV \
    --lod $LOD \
    --sequencing-delay $SEQUENCING_DELAY \
    --seed $SEED


# Rscript episewer.R \
#     --data-dir $WW_OUTPUT_DIR \
#     --n-sims $N_SIMS \
#     --seed $SEED