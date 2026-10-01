#!/bin/zsh

TREE_NAME="sim"

TREE_FILE=./data/remaster_trees/logs/$TREE_NAME.sampled.trees
ROOT_SEQUENCE_FILE=./data/root_sequences/test_root_sequence.txt
ROOT_FILE_TYPE=txt
OUTPUT_DIR=./data/simulated_sequences

# Evolutionary & Substitution Model Settings
TRANSITION_RATES=(0.5 3 0.5 0.5 6 1)   # Relative exchangeability rates: [AC, AG, AT, CG, CT, GT]
FREQUENCIES=(0.30 0.18 0.20 0.32)      # Base equilibrium frequencies: [A, C, G, T]
CLOCK_RATE=1e-3                        # Molecular clock rate (substitutions/site/unit time)
GAMMA_SHAPE=0.5                        # Shape parameter (alpha) for gamma site rate variation
GAMMA_CATEGORIES=4                     # Number of discrete rate categories
ROOT_DATE=2002                         # Temporal calendar date assigned to root node

uv run sequence_mutation.py \
    --tree_file $TREE_FILE \
    --root_sequence_file $ROOT_SEQUENCE_FILE \
    --root_file_type $ROOT_FILE_TYPE \
    --transition_rates $TRANSITION_RATES \
    --frequencies $FREQUENCIES \
    --clock_rate $CLOCK_RATE \
    --gamma_shape $GAMMA_SHAPE \
    --gamma_categories $GAMMA_CATEGORIES \
    --root_date $ROOT_DATE \
    --output_dir $OUTPUT_DIR

