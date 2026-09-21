#!/bin/zsh

TREE_FILE=./remaster_sims/logs/sim.sampled.trees
# ROOT_SEQUENCE_FILE=test_root_sequence.txt
ROOT_SEQUENCE_FILE=focal_sequences_strain.fasta
ROOT_FILE_TYPE=fasta
TRANSITION_RATES=(0.5 3 0.5 0.5 6 1)   # Relative exchangeability rates: [AC, AG, AT, CG, CT, GT]
FREQUENCIES=(0.30 0.18 0.20 0.32)
CLOCK_RATE=1.0
GAMMA_SHAPE=0.5
GAMMA_CATEGORIES=4
ROOT_DATE=2002
OUTPUT_DIR=sequences

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


