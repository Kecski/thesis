#!/bin/zsh

# PARAMETRIZATION:
# 
# Re = effective reproduction number
# delta = become uninfectious rate (rate at with each individual is removed from the infectious pool)
# s = sampling proportion (the probability with which infected individuals are sampled through time)
# 
# lambda = birth rate (per individual per unit time)
# mu = death rate
# psi = sampling rate
#
# 
# FORMULAS FOR CONVERSION
# 
# psi = s * delta
# mu = delta - s * delta
# lambda = Re * (delta)
#
# Re = lambda/(mu + psi)
# delta = (mu + psi)
# s = psi / (mu + psi)

# PARAMETERS FOR FIRST PART
RE=1.2
DELTA=$(echo "365 / 7" | bc -l)
S=$(echo "5e-3" | bc -l)

MAX_TIME=0.5

PSI=$(echo "$S * $DELTA" | bc -l)
MU=$(echo "$DELTA - ($S * $DELTA)" | bc -l)
LAMBDA=$(echo "$RE * $DELTA" | bc -l)


# PARAMETERS FOR SECOND PART
TREE_FILE=./remaster_sims/logs/sim.sampled.trees
ROOT_SEQUENCE_FILE=test_root_sequence.txt
ROOT_FILE_TYPE=txt
TRANSITION_RATES=(1 1 1 1 1 1)
FREQUENCIES=(0.25 0.25 0.25 0.25)
CLOCK_RATE=1.0
GAMMA_SHAPE=0.5
GAMMA_CATEGORIES=4
ROOT_DATE=2002
OUTPUT_DIR=sequences


# python3 simulation.py \
#     --output_name sim \
#     --birth_rate 1.6 1.4 \
#     --change_times_birth 3 \
#     --death_rate 1 0.8 \
#     --change_times_death 3 \
#     --sampling_rate 0.1 \
#     --max_time 10 \
#     --init_pop 1 \
#     --traj_type StochasticTrajectory \
#     --n_sims 1 \
#     --min_tips 1

# echo "Valameddig eljut"

uv run simulation.py \
    --output_name sim \
    --birth_rate "$LAMBDA" \
    --death_rate "$MU" \
    --sampling_rate "$PSI" \
    --max_time "$MAX_TIME" \
    --init_pop 1 \
    --traj_type StochasticTrajectory \
    --n_sims 1 \
    --min_tips 2

cd remaster_sims
beast sim.xml

cd ..
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
