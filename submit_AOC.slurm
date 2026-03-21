#!/bin/bash
#SBATCH --job-name=AOC
#SBATCH --output=AOC_%j.out
#SBATCH --error=AOC_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=48:00:00

set -euo pipefail

# Go to the directory you submitted from
cd "${SLURM_SUBMIT_DIR}"

# --- Environment setup (edit to match your cluster) ---
# Option A: modules + conda
#module purge || true
#module load miniconda || true

# Activate your env (edit env name)
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate AOC

# --- Run AOC (samples.csv-driven) ---
# Uses the Snakemake 9-compatible wrapper in run_AOC.sh.
# Adjust samples path if needed.
SAMPLES_CSV="${1:-samples.csv}"

# Use allocated CPUs
CORES="${SLURM_CPUS_PER_TASK:-16}"

# Jobs = how many rules Snakemake can have in flight (can be > cores if many IO-bound steps)
JOBS="${2:-20}"

bash run_AOC.sh \
  --samples "${SAMPLES_CSV}" \
  --cores "${CORES}" \
  --jobs "${JOBS}" \
  --latency-wait 60
