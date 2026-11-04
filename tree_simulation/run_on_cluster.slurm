#!/bin/bash
#SBATCH --job-name=test_simulation
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=20G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

PROJECT=$SCRATCH/thesis
cd "$PROJECT"

module purge
module load eth_proxy
if [ -n "${https_proxy:-}" ]; then
  p=${https_proxy#*://}; p=${p%%/*}
  export JAVA_TOOL_OPTIONS="-Dhttps.proxyHost=${p%%:*} -Dhttps.proxyPort=${p##*:} -Dhttp.proxyHost=${p%%:*} -Dhttp.proxyPort=${p##*:}"
fi

# ---- BEAST packages ----
module load stack/2024-06 gcc/12.2.0 openjdk beast2 libbeagle
# packagemanager -add BDMM-Prime
# packagemanager -add feast
# packagemanager -add remaster
# packagemanager -add SA 2.1.1

# ---- R environment ----
module load gcc/12.2.0 r/4.4.0 
# module load gdal
# mkdir -p "$SCRATCH/thesis/Rlibs" "$SCRATCH/thesis/renv"
# export R_LIBS_USER=$SCRATCH/Rlibs
# export RENV_PATHS_ROOT=$SCRATCH/thesis/renv
# Rscript -e 'if (!requireNamespace("renv", quietly=TRUE)) install.packages("renv", repos="https://cloud.r-project.org")'
# Rscript -e 'renv::restore(prompt = FALSE)'

# ---- Python environment ----
export UV_NO_CONFIG=1
# uv sync --frozen

cd tree_simulation
mkdir -p logs

# uv run snakemake --unlock
uv run snakemake --cores "$SLURM_CPUS_PER_TASK"
