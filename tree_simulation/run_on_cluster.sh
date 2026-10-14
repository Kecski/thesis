#!/bin/bash
#SBATCH --job-name=test_simulation
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

module purge
module load stack/2024-06 gcc/12.2.0 openjdk beast2 libbeagle
module load r/4.4.0
packagemanager -add bdmm-prime
packagemanager -add feast

export UV_NO_CONFIG=1
uv sync --frozen

mkdir -p logs

uv run snakemake