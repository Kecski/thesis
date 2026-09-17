#!/bin/zsh

python3 simulation.py \
    --output_name sim \
    --birth_rate 1.6 1.4 \
    --change_times_birth 3.7 \
    --death_rate 0.5 0.6 \
    --change_times_death 3.7 \
    --sampling_rate 0.5 \
    --max_time 10.0 \
    --init_pop 1 \
    --traj_type StochasticTrajectory \
    --n_sims 1

cd remaster_sims
./run_beast.sh