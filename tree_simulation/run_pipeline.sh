#!/usr/bin/env zsh
set -euo pipefail
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
#
# Configuration:
#   All tunable parameters live in pipeline.config (same directory by default).
#   Override the config path: ./run_pipeline.sh path/to/other.config
# ==============================================================================

SCRIPT_DIR="${0:A:h}"
CONFIG_FILE="${1:-$SCRIPT_DIR/pipeline.config}"

if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "Error: config file not found: $CONFIG_FILE" >&2
    exit 1
fi

source "$CONFIG_FILE"

# ------------------------------------------------------------------------------
# 1. EPIDEMIOLOGICAL PARAMETER CONVERSION
# ------------------------------------------------------------------------------
# Conversion Formulas:
#   psi    = s * delta           (sampling rate)
#   mu     = delta - (s * delta) (un-sampled clearance/death rate)
#   lambda = Re * delta          (transmission/birth rate)
# ------------------------------------------------------------------------------

DELTA=$(echo "$DELTA_EXPR" | bc -l)
S=$(echo "$S_EXPR" | bc -l)

PSI=$(echo "$S * $DELTA" | bc -l)
MU=$(echo "$DELTA - ($S * $DELTA)" | bc -l)
LAMBDA=$(echo "$RE * $DELTA" | bc -l)

echo "Derived rates -> lambda: $LAMBDA, mu: $MU, psi: $PSI"

# ------------------------------------------------------------------------------
# 2. PIPELINE EXECUTION
# ------------------------------------------------------------------------------

# --- Step 2.1: Generate Epidemiological Trajectory ---
uv run simulation.py \
    --output_name "$SIM_OUTPUT_NAME" \
    --birth_rate "$LAMBDA" \
    --death_rate "$MU" \
    --sampling_rate "$PSI" \
    --max_time "$MAX_TIME" \
    --init_pop "$SIM_INIT_POP" \
    --traj_type "$SIM_TRAJ_TYPE" \
    --n_sims "$SIM_N_SIMS" \
    --min_tips "$SIM_MIN_TIPS"

# --- Step 2.2: Simulate Phylogenetic Tree with BEAST ---
cd "$BEAST_WORKDIR"
beast "$BEAST_XML"
cd -

# --- Step 2.3: Simulate Nucleotide Mutations along the Tree ---
uv run sequence_mutation.py \
    --tree_file "$TREE_FILE" \
    --root_sequence_file "$ROOT_SEQUENCE_FILE" \
    --root_file_type "$ROOT_FILE_TYPE" \
    --transition_rates $TRANSITION_RATES \
    --frequencies $FREQUENCIES \
    --clock_rate "$CLOCK_RATE" \
    --gamma_shape "$GAMMA_SHAPE" \
    --gamma_categories "$GAMMA_CATEGORIES" \
    --root_date "$ROOT_DATE" \
    --output_dir "$OUTPUT_DIR"

echo "Pipeline complete."
