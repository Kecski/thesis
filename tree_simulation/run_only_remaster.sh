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
# Become uninfectious (delta) = (mu + psi)
# Sampling proportion (s) = psi / (mu + psi)

RE=(1.5 1 0.8)            # Target effective reproduction number
CHANGE_TIMES_BIRTH=(0.1 0.2)
DELTA=$(echo "365 / 7" | bc -l)        # Removal rate (~52.14/year = 7-day infectious period)
S=$(echo "5e-2" | bc -l)               # Sampling proportion (5% sampled)

MAX_TIME=0.3                           # Maximum simulation duration (0.5 time units / years)
MIN_SAMPLED=10                         # Minimum number of sampled genomes
OUTPUT_DIR="data/remaster_trees"
TREE_NAME="sim"
N_SIMS=2

# Derived Rate Calculations (evaluating floating-point arithmetic via 'bc')
PSI=$(echo "$S * $DELTA" | bc -l)      # Sampling rate (psi)
MU=$(echo "$DELTA - ($S * $DELTA)" | bc -l) # Death/un-sampled removal rate (mu)
LAMBDA=()
for n in "${RE[@]}"; do
    LAMBDA+=("$(echo "$n * $DELTA" | bc -l)")
done
#LAMBDA


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


uv run create_xml.py \
    --output_dir $OUTPUT_DIR \
    --output_name $TREE_NAME \
    --birth_rate "${LAMBDA[@]}" \
    --change_times_birth "${CHANGE_TIMES_BIRTH[@]}" \
    --death_rate "$MU" \
    --sampling_rate "$PSI" \
    --max_time "$MAX_TIME" \
    --init_pop 1 \
    --traj_type StochasticTrajectory \
    --n_sims $N_SIMS \
    --min_tips $MIN_SAMPLED

cd data/remaster_trees
beast $TREE_NAME.xml
cd -