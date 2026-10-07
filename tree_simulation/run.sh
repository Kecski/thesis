#!/bin/zsh
#
# ==============================================================================
# Pipeline: Epidemiological Trajectory, Tree & Sequence Mutation Simulator
# ==============================================================================
# Workflow Description:
#   1. Parameter Conversion: Translates high-level epidemiological parameters 
#      (Re, delta, s) into birth, death, and sampling rates (lambda, mu, psi).
#   2. Trajectory Simulation: Runs `simulation.py` to generate a stochastic
#      birth-death-sampling trajectory.
#   3. Tree Synthesis: Executes BEAST (`beast sim.xml`) inside `remaster_sims/`
#      to build phylogenetic trees from the simulated trajectory.
#   4. Sequence Mutation: Runs `sequence_mutation.py` to evolve nucleotide 
#      sequences down the generated phylogenetic tree.
#
# Dependencies:
#   - zsh shell
#   - bc (command-line calculator)
#   - uv (Python project & tool runner)
#   - BEAST (Phylogenetic tree simulator)
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. EPIDEMIOLOGICAL PARAMETER CONVERSION
# ------------------------------------------------------------------------------
# Biological Definitions:
#   Re    = Effective reproduction number
#   delta = Becoming uninfectious rate (1 / duration of infectiousness)
#   s     = Sampling proportion (probability an infected individual is sampled)
#
# Conversion Formulas:
#   psi    = s * delta           (sampling rate)
#   mu     = delta - (s * delta) (un-sampled clearance/death rate)
#   lambda = Re * delta          (transmission/birth rate)
# ------------------------------------------------------------------------------

# Input Biological Parameters
RE=(2.7 1.8 0.7 1.3)            # Target effective reproduction number
CHANGE_TIMES_BIRTH=(0.1 0.2 0.3)
DELTA=$(echo "365 / 12" | bc -l)         # Removal rate (~30.42/year = 7-day infectious period)
S=$(echo "2.6e-2" | bc -l)               # Sampling proportion (5% sampled)

MAX_TIME=0.4                           # Maximum duration of the sampling. The simulation will run an extra 0.1 years without sampling to gather wastewater data
MIN_SAMPLED=10                         # Minimum number of sampled genomes
REMASTER_OUTPUT_DIR=data/remaster_trees
TREE_NAME=sim
N_SIMS=1

# Derived Rate Calculations (evaluating floating-point arithmetic via 'bc')
PSI=$(echo "$S * $DELTA" | bc -l)      # Sampling rate (psi)
MU=$(echo "$DELTA - ($S * $DELTA)" | bc -l) # Death/un-sampled removal rate (mu)
LAMBDA=()
for n in "${RE[@]}"; do
    LAMBDA+=("$(echo "$n * $DELTA" | bc -l)")
done
#LAMBDA
#LAMBDA=$(echo "$RE * $DELTA" | bc -l)  # Birth/transmission rate (lambda)

# ------------------------------------------------------------------------------
# 2. SEQUENCE & EVOLUTIONARY SIMULATION PARAMETERS
# ------------------------------------------------------------------------------
SAMPLED_TREE_FILE=./data/remaster_trees/logs/$TREE_NAME.sampled.trees
ROOT_SEQUENCE_FILE=./data/root_sequences/test_root_sequence_1000.txt
# ROOT_FILE_TYPE=txt
SEQ_OUTPUT_DIR=./data/simulated_sequences

# Evolutionary & Substitution Model Settings
# TRANSITION_RATES=(0.5 3 0.5 0.5 6 1)   # Relative exchangeability rates: [AC, AG, AT, CG, CT, GT]
SUB_MODEL="HKY"
KAPPA=3.77
FREQUENCIES=(0.30 0.18 0.20 0.32)      # Base equilibrium frequencies: [A, C, G, T]
CLOCK_RATE=1.12e-3                        # Molecular clock rate (substitutions/site/unit time)
GAMMA_SHAPE=0.5                        # Shape parameter (alpha) for gamma site rate variation
GAMMA_CATEGORIES=6                     # Number of discrete rate categories
ROOT_DATE=2002                         # Temporal calendar date assigned to root node

# ------------------------------------------------------------------------------
# 3. WASTEWATER SIMULATION PARAMETERS
# ------------------------------------------------------------------------------
WW_OUTPUT_DIR=./data/wastewater
FULL_TREE_FILE=./data/remaster_trees/logs/$TREE_NAME.full.trees

N_DAYS=$(echo "$MAX_TIME * 365.25" | bc -l)
N_DAYS=${N_DAYS%.*}
WW_DATA_DIR="data"
SEED=42
SEQUENCING_DELAY=14

# Shedding Parameters
TOTAL_LOAD_PER_INFECTION=1.44e12
SHED_MAX_DAYS=12
SHED_MEAN=4
SHED_SD=2
TRANSPORT_MEAN=0.5
TRANSPORT_SD=0.5
FLOW_MEAN=1e8
MEAS_CV=0.3
LOD=1

# Episewer params
NUM_SEGMENTS=4

# ------------------------------------------------------------------------------
# 4. PIPELINE EXECUTION
# ------------------------------------------------------------------------------

# --- Step 3.1: Generate Epidemiological Trajectory ---
# Simulates stochastic population dynamics using calculated birth, death, and sampling rates.
uv run create_xml.py \
    --output_dir $REMASTER_OUTPUT_DIR \
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

# --- Step 3.2: Simulate Phylogenetic Tree with BEAST ---
# Navigates to working subfolder to run the XML tree simulation script.
cd data/remaster_trees
beast $TREE_NAME.xml
cd -

# --- Step 3.3: Simulate Nucleotide Mutations along the Tree ---
# Returns to root directory and mutates root sequence down the tree leaves.
uv run sequence_mutation.py \
    --tree_file $SAMPLED_TREE_FILE \
    --root_sequence_file $ROOT_SEQUENCE_FILE \
    --sub_model $SUB_MODEL\
    --frequencies $FREQUENCIES \
    --kappa $KAPPA \
    --clock_rate $CLOCK_RATE \
    --gamma_shape $GAMMA_SHAPE \
    --gamma_categories $GAMMA_CATEGORIES \
    --root_date $ROOT_DATE \
    --output_dir $SEQ_OUTPUT_DIR


uv run wastewater_data.py \
    --tree $FULL_TREE_FILE \
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


Rscript episewer.R \
    --data-dir $WW_OUTPUT_DIR \
    --n-sims $N_SIMS \
    --num-segments $NUM_SEGMENTS \
    --seed $SEED

