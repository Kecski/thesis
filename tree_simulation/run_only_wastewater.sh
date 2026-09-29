#!/bin/zsh

N_DAYS=110
SAMPLE_DAYS_EVERY=1
N_REPLICATES=1
OUTPUT_DIR="episewer/test"
EPISEWER_PREFIX="episewer"

TOTAL_LOAD_PER_INFECTION=1e9
SHED_MAX_DAYS=21
SHED_MEAN=4
SHED_SD=3
TRANSPORT_MEAN=7
TRANSPORT_SD=2
FLOW_MEAN=1e8
MEAS_CV=0.3
LOD=1
SEED=42

uv run wastewater_data.py remaster_sims/logs/sim.full.trees \
    --n-days $N_DAYS \
    --sample-days-every $SAMPLE_DAYS_EVERY \
    --n-replicates $N_REPLICATES \
    --output_dir $OUTPUT_DIR \
    --episewer-prefix $EPISEWER_PREFIX \
    --total-load-per-infection $TOTAL_LOAD_PER_INFECTION \
    --shed-max-days $SHED_MAX_DAYS \
    --shed-mean $SHED_MEAN \
    --shed-sd $SHED_SD \
    --transport-mean $TRANSPORT_MEAN \
    --transport-sd $TRANSPORT_SD \
    --flow-mean $FLOW_MEAN \
    --meas-cv $MEAS_CV \
    --lod $LOD \
    --seed $SEED


cd episewer

Rscript simulate_wastewater_metrics.R