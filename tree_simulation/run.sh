#!/bin/zsh

python3 simulation.py \
    --output_name sim \
    --birth_rate 1.4 1.2 \
    --change_times_birth 3.7 \
    --death_rate 1.0 0.8 \
    --change_times_death 3.7 \
    --sampling_rate 1.0 \
    --max_time 10.0 \
    --init_pop 1 \
    --n_sims 1

cd remaster_sims
./run_beast.sh