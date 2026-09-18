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

RE=1.1
DELTA=$(echo "365 / 7" | bc -l)
S=$(echo "1e-4" | bc -l)

MAX_TIME=0.2

PSI=$(echo "$S * $DELTA" | bc -l)
MU=$(echo "$DELTA - ($S * $DELTA)" | bc -l)
LAMBDA=$(echo "$RE * $DELTA" | bc -l)
echo "psi = $PSI"
echo "mu = $MU"
echo "lambda = $LAMBDA"


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
    --min_tips 1

cd remaster_sims
beast sim.xml
